##### SET UP #####

## load packages
library(plyr)
library(tidyverse)
library(lubridate)
library(sandwich)
library(car)
library(plyr)
library(MASS)
library(nlme)
library(ggplot2)+theme_set(theme_bw()+
                             theme(panel.grid = element_blank(),
                                   strip.background = element_rect(fill="white"),
                                   axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
                                   axis.text.x = element_text(colour = "black", size = 8),
                                   legend.text = element_text(size = 8, face ="bold", colour ="black"),
                                   legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
                                   axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
                                   legend.title = element_text(size = 14, colour = "black", face = "bold"),
                                   legend.key=element_blank()))
  
## read in data
tempdata =read.csv("2022_temp_logger_data.csv")
logger_elevation = read.csv("logger_elevation.csv")
tide_obs = read.csv("gov_canada_tide_observations.csv")

## make color vector for site
site.colors = c("#39B54A", "darkgreen", "plum", "#A30262")


##### FORMAT THE TEMPERATURE LOGGER DATA #####
## remove bad data (before deployment temps and first 12h outside)
tempdata = subset(tempdata, tempdata$keep=="yes")

## set variable format
tempdata$datetime.utc = as.POSIXct(tempdata$datetime.utc, format="%Y-%m-%d %H:%M")
tempdata$datetime.pdt = as.POSIXct(tempdata$datetime.pdt, format="%Y-%m-%d %H:%M")
tempdata$temp = as.numeric(tempdata$temp)

## split columns
tempdata = separate(data = tempdata, 
                         col = site, 
                         into = c("site","position"), 
                         sep = "-")

## rename sites
site_id <- c(gw="GW", cp="CP", tb="TB", scp="SCP")
tempdata$site <- as.character(site_id[tempdata$site])

## order the sites
tempdata$site = factor(tempdata$site, levels=c("TB","SCP","GW","CP"))



##### FORMAT THE TIDE DATA #####
## rename observations columns
names(tide_obs)[names(tide_obs)=="observations(m)"]<-"tide_height_m"

## split the dates column
tide_obs = tide_obs %>%
  separate(col=Date,
           sep=" ",
           into=c("date", "time", "timezone"))

## get date and time together to make datetime colums
tide_obs$datetime.pdt = paste(tide_obs$date,tide_obs$time)

## format datetime 
tide_obs$datetime.pdt = as.POSIXct(tide_obs$datetime.pdt, format="%Y-%m-%d %H:%M")

##### PLOT THE TEMPERATURE DATA #####
## temp data all separate
ggplot(tempdata, aes(x=datetime.pdt, y=temp, color=site))+
  geom_line(cex=0.5)+
  ylab("Reccorded temperature (C)")+
  xlab("Date and time")+
  facet_grid(site~position)+
  scale_color_manual(values=site.colors)



##### MERGE THE DATAFRAMES TOGETHER AND ADD UNDERWATER VARIABLE  ######
## merge
full_tide = full_join(tempdata, tide_obs)
                  
full_tide = full_join(full_tide, logger_elevation)


## add underwater variable
full_tide$underwater <- if_else(full_tide$elevation_m < full_tide$tide_height_m, "yes", "no")

##### PLOT THE LOGGER TEMPERATURE BY UNDER/ABOVE WATER #####
full_tide=na.omit(full_tide)

## temp data all separate
ggplot(full_tide, aes(x=datetime.pdt, y=temp, color=underwater))+
  geom_point(cex=0.3)+
  facet_grid(site+position~.)+
  labs(x="Date", y="Reccorded Temperature (C)", color="Underwater?")


##### ARE THE SITES DIFFERENT TEMPERATURES? #####

lm.temp = lm(temp~site*position*datetime.pdt*underwater, data=full_tide)
summary(lm.temp) ## GW and CP water colder than SCP and TB

##### SELECT THE TOP 5 WARMEST AND COLDEST MEASURE FOR EVERY SITE AND DAY #####

## make groups
full_tide$daygroups = paste0(full_tide$site, full_tide$position, full_tide$date)
#write.csv(full_tide, paste0(path,"temp_profile_dataframe.csv"))

## get list of groups
groups = c(unique(full_tide$daygroups))

## get temp extremes
## sort data by relative abundance
sorted = full_tide[order(-full_tide$temp),]

## make empty dataframe to store data
temp.df = NULL


## start loop
for(i in groups) {
  for(j in i) {
    ## subset dataframe by samples
    sample = subset(sorted, sorted$daygroup %in% c(j))
    
    ## get top 10 warmest temps
    top = sample[c(1:5),]
    
    ## get top 10 coldest temps
    coldest = nrow(sample)
    start.b5 = coldest - 4
    
    bottom = sample[c(start.b5:coldest),]
    
    ## save list of top  abundance taxa
    t.tmp <- top
    b.tmp <- bottom
    
    ## add top of bottom indicator
    t.tmp$place = "top5"
    b.tmp$place = "bottom5"
    
    ## combine both top and bottom dfs
    new.temps = merge(t.tmp, b.tmp, all=TRUE)
    
    ## add group indicator
    new.temps$group = j
    
    ## bind with temporary dataframe
    temp.df <- rbind.fill(temp.df, new.temps)
    
    ## close loop 
  }
}


##### calculate the mean +/- SD for each group
temps.grouped = ddply(temp.df, c("site","position","date", "place"), 
                        summarise,
                        N.temp = length(temp), ## sample size
                        mean.temp = mean(temp), ## mean
                        sd.temp = sd(temp))%>% ## standard deviation) %>%
  mutate(se.temp = sd.temp/sqrt(N.temp), ## standard error
         lci.temp = mean.temp - sd.temp, ## lower bound of 95% confidence interval
         uci.temp = mean.temp + sd.temp) ## upper bound of 95% confidence interval



##### MAKE DOT PLOT WITH TOP AND BOTTOM TEMPS BY DAY BY SITE #####
## plot mean water temps by date
ggplot(temps.grouped, aes(x=as.Date(date), y=mean.temp, colour=site))+ 
  geom_errorbar(aes(ymin=lci.temp, ymax=uci.temp), width=1) +
  geom_point(cex=2)+
  geom_line()+
  facet_grid(position~place, scales = "free")+
  scale_color_manual(values=c("#A30262", "plum", "darkgreen","#39B54A"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x="Date", y="Mean of the 5 temperature extemes")

## only plot lower temps
bottoms = subset(temps.grouped, temps.grouped$place=="bottom5")
ggplot(bottoms, aes(x=date, y=mean.temp, colour=site, group = site))+ 
  geom_errorbar(aes(ymin=lci.temp, ymax=uci.temp), width=0.4) +
  geom_point(cex=2)+
  geom_line(linewidth=1)+
  facet_grid(paste0(position, "-intertidal")~.)+
  scale_color_manual(values=c("#A30262", "plum", "darkgreen","#39B54A"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(y="Mean of the lowest 5 temperatures (C) by day")
ggsave(filename = ("mean_bottom5_temps_site.pdf"), width=17, height=8, units="in", path=path)

tops = subset(temps.grouped, temps.grouped$place=="top5")
ggplot(tops, aes(x=date, y=mean.temp, colour=site, group = site))+ 
  geom_errorbar(aes(ymin=lci.temp, ymax=uci.temp), width=0.4) +
  geom_point(cex=2)+
  geom_line(linewidth=1)+
  facet_grid(paste0(position, " intertidal")~.)+
  scale_color_manual(values=c("#A30262", "plum", "darkgreen","#39B54A"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(y="Mean of the highest 5 temperatures (C) by day")
ggsave(filename = ("mean_top5_temps_site.pdf"), width=17, height=8, units="in", path=path)



##### LINEAR REGRESSION MODEL #####
temp.df = temp.df %>% separate(col=date, into=c("year","month","day"), sep="-")

top5 = subset(temp.df, temp.df$place=="top5")

mod2 = lm(as.numeric(temp)~position*site, data=top5)
summary(mod2)
capture.output(summary(mod2),file=paste0(path,"lm_top5_temps_site.txt"))



bottom5 =subset(temp.df, temp.df$place=="bottom5")

mod2 = lm(as.numeric(temp)~position*site, data=bottom5)
summary(mod2)
capture.output(summary(mod2),file=paste0(path,"lm_bottom5_temps_site.txt"))


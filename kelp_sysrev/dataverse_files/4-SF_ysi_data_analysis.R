##### WORKSPACE SETUP #####
## load packages
library(plyr)
library(tidyverse)
library(ggpubr)
library(ggplot2)+theme_set(theme_bw()+
                             theme(strip.background = element_rect(fill="white"),
                                   axis.text.y = element_text(colour = "black", size = 12, face = "bold"),
                                   axis.text.x = element_text(colour = "black", face = "bold", size = 12),
                                   legend.text = element_text(size = 8, face ="bold", colour ="black"),
                                   legend.position = "right", axis.title.y = element_text(face = "bold", size = 14),
                                   axis.title.x = element_text(face = "bold", size = 14, colour = "black"),
                                   legend.title = element_text(size = 14, colour = "black", face = "bold"),
                                   legend.key=element_blank(),
                                   #axis.ticks = element_blank()
                                   ))


## set the position dodge for ggplot
pd = position_dodge(0.2) # move them .05 to the left and right

## make color vector for site
site.colors = c("seagreen1", "#39B54A", "darkgreen", "plum", "#A30262")

## load in data
abiotic = read.csv("ysi_data_both_years.csv")

lmnotp = subset(abiotic, abiotic$water_type !="tidepool")


## group the data by site, sampling date, and water type
abiotic.grouped = ddply(abiotic, c("site_id", "water_type", "day","month","year"), 
                       summarise,
                       N.conduct = length(salinity), ## sample size
                       mean.conduct = mean(salinity), ## mean
                       sd.conduct = sd(salinity),
                      N.temp = length(water_temp_C), ## sample size
                      mean.temp = mean(water_temp_C), ## mean
                      sd.temp = sd(water_temp_C))%>% ## standard deviation) %>%
  mutate(se.conduct = sd.conduct/sqrt(N.conduct), ## standard error
         lci.conduct = mean.conduct - sd.conduct, ## lower bound of 95% confidence interval
         uci.conduct = mean.conduct + sd.conduct,
         se.temp = sd.temp/sqrt(N.temp), ## standard error
         lci.temp = mean.temp - (sd.temp), ## lower bound of 95% confidence interval
         uci.temp = mean.temp + (sd.temp)) ## upper bound of 95% confidence interval

#write.csv(abiotic.grouped, "summarised_abiotic_2021_and_2022.csv")


##### ADD JULIAN DAY TO DATAFRAME TO HAVE SEASONALITY INFO LATER #####
abiotic.grouped$jul.date = paste0(as.numeric(abiotic.grouped$year), "-", 
                                  as.numeric(abiotic.grouped$month), "-", 
                                  as.numeric(abiotic.grouped$day))
abiotic.grouped$jul.date = as.Date(abiotic.grouped$jul.date)
abiotic.grouped$jul.date = as.POSIXlt(abiotic.grouped$jul.date)$yday

## order the sites
abiotic.grouped$site_id = factor(abiotic.grouped$site_id, levels=c("lhp", "tb","scp","gw","cp"))

## remove tidepools
abiotic.notp = subset(abiotic.grouped, abiotic.grouped$water_type !="tidepool")


##### GRAPH - SALINITY #####

## plot mean water salinity by date
salprof=ggplot(abiotic.notp, aes(x=jul.date, y=mean.conduct, colour=site_id))+ 
  geom_errorbar(aes(ymin=lci.conduct, ymax=uci.conduct), 
                width=1) +
  geom_point(cex=3)+
  geom_line(linewidth=1.3)+
  ylab("Mean Conductance")+
  xlab(" ")+
  labs(colour="Site")+
  facet_wrap(.~year)+
  scale_colour_manual(values = site.colors)+
  theme(legend.background = element_rect(fill="white",size=0.5, linetype="solid",colour ="grey"))+
  theme(axis.text.x=element_text(angle=90, hjust=1))
salprof


##### GRAPH TEMPERATURE ####
## plot mean water salinity by date
tempprof=ggplot(abiotic.notp, aes(x=jul.date, y=mean.temp, colour=site_id))+ 
  geom_errorbar(aes(ymin=lci.temp, ymax=uci.temp), 
                width=1) +
  geom_point(cex=3)+
  geom_line(linewidth=1.3)+
  ylab("Mean Temperature (C)")+
  xlab(" ")+
  labs(colour="Site")+
  facet_wrap(.~year)+
  scale_colour_manual(values = site.colors)+
  theme(legend.background = element_rect(fill="white",size=0.5, linetype="solid",colour ="grey"))+
  theme(axis.text.x=element_text(angle=90, hjust=1))
tempprof


##### arrange graphs #####
ggarrange(salprof, tempprof, ncol=1, nrow=2, labels=c("A", "B"))
ggsave(filename = ("ysi_data.pdf"), width=9.4, height=9.4, units="in", path=path)

##### CORRELATE TEMPERATURE AND SALINITY ######

## run corelation test
cor.test(abiotic.grouped$mean.conduct, abiotic.grouped$mean.temp)


## plot temperature vs salinity
ts=ggplot(abiotic.grouped, aes(x=mean.conduct, y=mean.temp))+
  geom_point()+
  geom_smooth(method='lm', formula= y~x, se=F)+
  labs(y="Mean Temperature (C)", x="Mean Conductance")+
  annotate(geom="text", x=15, y =11.5, 
           label="r = -0.5102358, p < 0.001")


##### CORELATE TEMPERATURE AND SALINITY WITH JULIAN DAY #####
## temperature
cor.test(abiotic.grouped$jul.date, abiotic.grouped$mean.temp)
jt=ggplot(abiotic.grouped, aes(y=mean.temp, x=jul.date))+
  geom_point()+
  geom_smooth(method='lm', formula= y~x, se=F)+
  labs(y="Mean Temperature (C)", x="Julian day of the year")+
  annotate(geom="text", x=185, y =11.5, label="r = 0.7195482 , p < 0.001")

## salinity
cor.test(abiotic.grouped$jul.date, abiotic.grouped$mean.conduct)
js=ggplot(abiotic.grouped, aes(y=mean.conduct, x=jul.date))+
  geom_point()+
  geom_smooth(method='lm', formula= y~x, se=F)+
  labs(y="Mean Conductance", x="Julian day of the year")+
  annotate(geom="text", x=185, y =11.5, label="r = -0.4664766 , p < 0.001")


##### SAVE THE CORELATION PLOTS #####
ggarrange(js, jt, ts, ncol=3, nrow=1)


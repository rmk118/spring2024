##### WORKSPACE SETUP  ####
## load packages
library(tidyverse)
library(plyr)
library(stats)
library(phyloseq)
library(ggpubr)
library(ggh4x)
library(car)
library(rstatix)
library(ggplot2); theme_set(theme_bw()+
                              theme(panel.grid = element_blank(),
                                    strip.background = element_rect(fill="white"),
                                    axis.text = element_text(size = 12, colour = "black", face="bold"),
                                    axis.title = element_text(size=15),
                                    strip.text = element_text(color="black", size=12),
                                    legend.text=element_text(size=12),
                                    axis.line = element_line(colour = "black")))

dephyloseq = function(phylo_obj){
  
  ## get the metadata
  meta = as.data.frame(as.matrix(phylo_obj@sam_data))
  
  ## how many metadata columns you have 
  metacols = ncol(meta)+1
  
  ## get out the otu table 
  otu = as.data.frame(t(as.matrix(phylo_obj@otu_table)))
  
  ## merge the metadata and otu table by the rownames (sample ids from the Illumina sequencing data)
  mo = merge(meta, otu, by=0)
  
  ## get out the taxonomy file 
  tax = as.data.frame(phylo_obj@tax_table)
  
  ## get the ASV ID out. This the matches the placeholder ASV ID in the OTU table
  tax = tax %>% rownames_to_column(var="asv_id")
  
  ## pivot longer to be able to match the ASVs in the OTU table to the taxonomy table 
  mo = mo %>% pivot_longer(cols = -c(1:metacols), names_to = "asv_id", values_to="asv_abundance")
  
  ## Join the metadata and otu table with the taoxnomy table 
  mot = full_join(mo, tax)
  
  ## Specify the output for the dephyloseq funciton 
  output = mot
}

## set color vector
salinity.colors = c("cyan","dodgerblue2","navy")

## read in data
wetweight = read.csv("lab_condition_metadata.csv")

##### FORMAT DATA FOR ANALYSIS #####

## fix variable type
wetweight$exp_day = as.numeric(wetweight$exp_day)
wetweight$culture_round = as.numeric(wetweight$culture_round)

## add tracking variable 
wetweight$meristem_track = paste0(wetweight$meristem_id, "-",wetweight$culture_round)
wetweight$sample_id = paste0(wetweight$salinity_percent, "-", wetweight$culture_round, "-",wetweight$meristem_id)


##### DOES THE DAMMAGE REMAIN ONCE IT OCCURS? #####
## only keep days 1, 2, and 5
relv.days = subset(wetweight, wetweight$exp_day %in% c(1, 2, 5))

## remove one with blisterign at day 1
relv.days = subset(relv.days, relv.days$meristem_track!="4-2")

sal.labs <- c("10", "20", "full strength")
names(sal.labs) <- c("10", "20", "30")

ct <- c("0", "1", "4")

blistering=ggplot(relv.days, aes(x=sample_id, y=as.character(exp_day), fill=as.factor(blistering_lab)))+
  geom_tile()+
  scale_fill_manual(labels=c("no", "yes"), values=c("grey90", "orange"))+
  facet_nested(.~aquarium_salinity+culture_round, labeller = labeller(aquarium_salinity=sal.labs), scales = "free") +
  theme(panel.spacing = unit(0,"line")) +
  theme(axis.text.x = element_blank(), panel.spacing = unit(0,"lines"))+
  labs(x=" ", y="Time in lab (days)", fill="Presence of Blistering")+
  theme(legend.position = "none")+ 
  scale_y_discrete(labels= ct)


blistering
ggsave(filename = ("meristem_dammage_blistering.pdf"), width=11, height=3.8, units="in", path=path)


##### FISHER TEST - DOES LOWER SALINITY LEAD TO MORE DAMMAGE #####

## get dataframe for fisher

fisher.df = relv.days[,c("aquarium_salinity", "exp_day", "blistering_lab")]

fisher.df = subset(fisher.df, fisher.df$exp_day=="5")

fisher.df = ddply(fisher.df, c("aquarium_salinity"),
                  summarise,
                  sum = sum(blistering_lab),
                  all = length(blistering_lab))

## calculate yes/no for matrix
fisher.df$blistered = fisher.df$sum
fisher.df$notblistered = fisher.df$all - fisher.df$sum

## remove unnecessary columns 
fisher.df = fisher.df[,c("aquarium_salinity", "blistered", "notblistered")]

## format dtaa
fisher.df = fisher.df %>% 
  column_to_rownames(var="aquarium_salinity")

## caulcuate stat
fisher.test(fisher.df)

## post-hoc
pairwise_fisher_test(as.matrix(fisher.df), p.adjust.method = "BH")


##### ANOVA - did the 10 psu take on more water than 20 and 30? #####
## group by day 1, d2, and d5 for the same kelp, by same salinity and culture round 
d5 = subset(wetweight, wetweight$exp_day=="5")

## get ration of dry weiht over wet weight 
d5$dw.to.ww = (d5$post_dry_weight_g - d5$petri_dish_weight_g)/(d5$pre_dry_weight_g - d5$petri_dish_weight_g)

## remove datapoints with no dry weight over wet weight
d5 = subset(d5, d5$dw.to.ww>0)

## ANOVA for dry weight over wet weight ration by group
a1 = aov(dw.to.ww~as.factor(aquarium_salinity), data=d5)
a1
summary(a1)
capture.output(summary(a1),file=paste0(path,"anova_dw_to_ww_d5.txt"))
qplot(residuals(a1))
lt=leveneTest(dw.to.ww~as.factor(culture_round)/as.factor(aquarium_salinity), data=d5)
lt
capture.output(lt,file=paste0(path,"levene_dw_to_ww_d5.txt"))

## The responses for each factor level have a normal population distribution.
## These distributions have the same variance. (yes, Levenes)
## The data are independent. (between salinities, yes)

hsd=TukeyHSD(a1, conf.level = 0.95)
hsd 
capture.output(hsd,file=paste0(path,"tukey_dw_to_ww_d5.txt"))


##### PLOTS - did the 10 psu take on more water than 20 and 30? #####

salinity = c("10", "20", "full strength")

## plot ratio in dry weight to wet weight  
dw=ggplot(d5, aes(x=as.factor(aquarium_salinity), y=dw.to.ww, fill=as.factor(aquarium_salinity), color=as.factor(aquarium_salinity)))+
  geom_boxplot(linewidth=0.5, outlier.size = 5, alpha = 0.3)+
  scale_fill_manual(values=c(salinity.colors))+
  scale_color_manual(values=c(salinity.colors))+
  labs(x=" ", y="Ratio of dry weight to wet 
weight on day four", fill="Salinity")+
  guides(colour = "none")+
  ylim(0,0.175)+
  #20-10 0.012603587 0.004610629 0.02059655 0.0009605
  annotate("text", x = 1.5, y = 0.165, label = "p = 0.006", size=4, color="black")+
  annotate("segment", x = 1.1, xend = 1.9, y = 0.16, yend = 0.16, colour = "black") +
  #30-20 0.009043122 0.001050163 0.01703608 0.0227933
  annotate("text", x = 2.5, y = 0.165, label = "p = 0.066", size=4, color="black")+
  annotate("segment", x = 2.1, xend = 2.9, y = 0.16, yend = 0.16, colour = "black") +
  #30-10 0.021646709 0.013653750 0.02963967 0.0000000
  annotate("text", x = 2, y = 0.175, label = "p <0.001", size=4, color="black")+
  annotate("segment", x = 1.1, xend = 2.9, y = 0.169, yend = 0.169, colour = "black")+
  theme(legend.position = "none")+
  scale_x_discrete(labels= salinity)

dw

ggsave(filename = ("day5_dw_to_ww_boxplot.pdf"), width=4.8, height=5, units="in", path=path)



##### How do different salinity affect growth? ####

d1 = subset(wetweight, wetweight$exp_day=="1")

## get day 1 weight accross entire datafram e
d1$d1w = (d1$wet_weight_g)
d1 = d1[,c("d1w", "meristem_track")]

relv.days = full_join(d1, wetweight)

## normalize as % of starting weight (d1 weight)
relv.days$norm.ww = (relv.days$wet_weight_g)/(relv.days$d1w)

relv.days = subset(relv.days, relv.days$exp_day %in% c("1","5") & relv.days$culture_round %in% c("3","4","5","6","7","8"))


## only keep day 5 and make different plot 
d5 = subset(relv.days, relv.days$exp_day=="5")

## ANOVA for dry weight over wet weight ration by group
a2 = aov(norm.ww~as.factor(aquarium_salinity), data=d5)
a2
summary(a2)
capture.output(summary(a2),file=paste0(path,"anova_growth_d5vsd1.txt"))
qplot(residuals(a2))
lt=leveneTest(norm.ww~as.factor(aquarium_salinity)*as.factor(culture_round), data=d5)
capture.output(lt,file=paste0(path,"levene_growth_d5vsd1.txt"))
lt

## The responses for each factor level have a normal population distribution.
## These distributions have the same variance. (yes, Levenes)
## The data are independent. (between salinities, yes)

tsd=TukeyHSD(a2, conf.level = 0.95)
tsd
capture.output(tsd,file=paste0(path,"tukey_growth_d5vsd1.txt"))

## plot ratio in dry weight to wet weight  
ww=ggplot(d5, aes(x=as.factor(aquarium_salinity), y=norm.ww, fill=as.factor(aquarium_salinity), color=as.factor(aquarium_salinity)))+
  geom_boxplot(linewidth=0.5, outlier.size = 5, alpha = 0.3)+
  scale_fill_manual(values=c(salinity.colors))+
  scale_color_manual(values=c(salinity.colors))+
  labs(x=" ", y="Ratio of wet weight on day 
four vs day zero", fill="Salinity")+
  guides(colour = "none", fill = "none")+
  ylim(0.9,1.51)+
  scale_x_discrete(labels= salinity)

ggsave(filename = ("ww_d1_vd_d5_boxplot.pdf"), width=4.8, height=5, units="in", path=path)
ww




##### ANOVA- HOW DOES SALINTIY AFFECT PAM (Y) #####
## only look at day 5

## ANOVA for dry weight over wet weight ration by group
a3 = aov(as.numeric(Y)~as.factor(aquarium_salinity), data=d5)
a3
summary(a3)
capture.output(summary(a3),file=paste0(path,"anova_pam_d5.txt"))
qplot(residuals(a3))
lt=leveneTest(as.numeric(Y)~as.factor(aquarium_salinity)*as.factor(culture_round), data=d5)
capture.output(lt,file=paste0(path,"levene_test_pam_d5.txt"))
lt

tsd=TukeyHSD(a3, conf.level = 0.95)
tsd
capture.output(tsd,file=paste0(path,"tukey_pam_d5.txt"))


## plot pam  
pam=ggplot(d5, aes(x=as.factor(aquarium_salinity), y=as.numeric(Y), fill=as.factor(aquarium_salinity), color=as.factor(aquarium_salinity)))+
  geom_boxplot(linewidth=0.5, outlier.size = 5, alpha = 0.3)+
  scale_fill_manual(values=c(salinity.colors))+
  scale_color_manual(values=c(salinity.colors))+
  labs(x=" ", y="Estimated effective quantum 
yeild (Fv/Fm) on day four", fill="Salinity")+
  guides(colour = "none")+
  theme(legend.position = "none")+
  scale_x_discrete(labels= salinity)

ggsave(filename = ("day5_pam_boxplot.pdf"), width=4.8, height=5, units="in", path=path)

##### group salinity plots together #####
p1= ggarrange(dw, ww, pam, ncol=3, nrow=1, labels=c("B", "C","D"))

ggarrange(blistering, p1,
          ncol=1, nrow=2,
          heights=c(0.8, 1),
          labels=c("A", " "))

ggsave(filename = ("meristem_dammage_boxplots.pdf"), 
       width=14.18, height=7.78, units="in", path=path)

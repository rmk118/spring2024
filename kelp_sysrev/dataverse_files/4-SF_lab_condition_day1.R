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
## only keep 2022 freshet samples
#wetweight = subset(wetweight, wetweight$X2022_freshet_sample=="yes" & wetweight$pam_notes!="demerestiamaybeinaquaria")

## fix variable type
#wetweight$exp_day = as.numeric(wetweight$exp_day)
wetweight$culture_round = as.numeric(wetweight$culture_round)

## add tracking variable 
wetweight$meristem_track = paste0(wetweight$meristem_id, "-",wetweight$culture_round)
wetweight$sample_id = paste0(wetweight$salinity_percent, "-", wetweight$culture_round, "-",wetweight$meristem_id)


##### DOES THE DAMMAGE REMAIN ONCE IT OCCURS? #####
## remove one with blisterign at day 1
relv.days = subset(wetweight, wetweight$meristem_track!="4-2")

sal.labs <- c("10", "20", "full strength")
names(sal.labs) <- c("10", "20", "30")

#ct <- c("0", "1", "4")

blistering=ggplot(relv.days, aes(x=sample_id, 
                                 y=as.character(exp_day), 
                                 fill=as.factor(blistering_lab)))+
  geom_tile()+
  scale_fill_manual(labels=c("no", "yes"), values=c("grey90", "orange"))+
  facet_nested(.~aquarium_salinity+culture_round, labeller = labeller(aquarium_salinity=sal.labs), scales = "free") +
  theme(panel.spacing = unit(0,"line")) +
  theme(axis.text.x = element_blank(), panel.spacing = unit(0,"lines"))+
  labs(x=" ", y="Time in lab (days)", fill="Presence of Blistering")+
  theme(legend.position = "none")#+scale_y_discrete(labels= ct)


blistering
ggsave(filename = ("meristem_dammage_blistering_alldays.pdf"), width=11, height=3.8, units="in", path=path)


##### FISHER TEST - DOES LOWER SALINITY LEAD TO MORE DAMMAGE #####

## get dataframe for fisher

fisher.df = relv.days[,c("aquarium_salinity", "exp_day", "blistering_lab")]

fisher.df = subset(fisher.df, fisher.df$exp_day=="2")

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
pairwise_fisher_test(as.matrix(fisher.df), p.adjust.method = "fdr")


##### How do different salinity affect growth? ####

d1 = subset(wetweight, wetweight$exp_day=="1")

## get day 1 weight accross entire datafram e
d1$d1w = (d1$wet_weight_g)
d1 = d1[,c("d1w", "meristem_track")]

relv.days = full_join(d1, wetweight)

## normalize as % of startign weight (d1 weight)
relv.days$norm.ww = (relv.days$wet_weight_g)/(relv.days$d1w)

## only keep day 5 and make different plot 
d2 = subset(relv.days, relv.days$exp_day=="2")

## ANOVA for dry weight over wet weight ration by group
a2 = aov(norm.ww~as.factor(aquarium_salinity), data=d2)
a2
summary(a2)
capture.output(summary(a2),file=paste0(path,"anova_growth_d2vsd1.txt"))
qplot(residuals(a2))
lt=leveneTest(norm.ww~as.factor(aquarium_salinity)*as.factor(culture_round), data=d2)
capture.output(lt,file=paste0(path,"levene_growth_d2vsd1.txt"))
lt

## The responses for each factor level have a normal population distribution.
## These distributions have the same variance. (yes, Levenes)
## The data are independent. (between salinities, yes)

tsd=TukeyHSD(a2, conf.level = 0.95)
tsd
capture.output(tsd,file=paste0(path,"tukey_growth_d2vsd1.txt"))

## plot ratio in dry weight to wet weight  
ww=ggplot(d2, aes(x=as.factor(aquarium_salinity), y=norm.ww, 
                  fill=as.factor(aquarium_salinity), color=as.factor(aquarium_salinity)))+
  geom_boxplot(linewidth=0.5, outlier.size = 5, alpha = 0.3)+
  scale_fill_manual(values=c(salinity.colors))+
  scale_color_manual(values=c(salinity.colors))+
  labs(x=" ", y="Ratio of wet weight on day 
one vs day zero", fill="Salinity")+
  guides(colour = "none", fill = "none")+
  ylim(0.9,1.51)

ggsave(filename = ("ww_d1_vd_d2_boxplot.pdf"), width=4.8, height=5, units="in", path=path)
ww




##### ANOVA- HOW DOES SALINTIY AFFECT PAM (Y) #####
## only look at day 5

## ANOVA for dry weight over wet weight ration by group
a3 = aov(as.numeric(Y)~as.factor(aquarium_salinity), data=d2)
a3
summary(a3)
capture.output(summary(a3),file=paste0(path,"anova_pam_d2.txt"))
qplot(residuals(a3))
lt=leveneTest(as.numeric(Y)~as.factor(aquarium_salinity)*as.factor(culture_round), data=d2)
capture.output(lt,file=paste0(path,"levene_test_pam_d2.txt"))
lt

tsd=TukeyHSD(a3, conf.level = 0.95)
tsd
capture.output(tsd,file=paste0(path,"tukey_pam_d2.txt"))


## plot pam  
pam=ggplot(d2, aes(x=as.factor(aquarium_salinity), y=as.numeric(Y), fill=as.factor(aquarium_salinity), color=as.factor(aquarium_salinity)))+
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

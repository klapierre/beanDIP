####setup####
library(tidyverse)
setwd("~/Dropbox/bean_dip_2018-2024/data")

treatment_key <- c("67" = "untreated","27" = "treated")
site_key <- c("Clarksville" = "C","Keedysville" = "K","Poplar Hill"= "PH","Wye"="W")

####import data####
harvest<-read.csv("clean_data/clean_harvest_2019.csv")
harvest$site<-factor(harvest$site, levels = c("K", "C", "W", "PH"))
harvest$variety<-as.factor(harvest$variety)

#for comparing only treated v untreated
# KM<-harvest[(harvest$variety %in% c("27","67")),]
# KM$treatment<-recode(KM$variety, !!!treatment_key)
# KM$treatment<-droplevels(KM$treatment)
# KM$treatment<-factor(KM$treatment, levels = c("untreated", "treated"))

damage<-read.csv("clean_data/clean_damage_2019.csv")
damage$site<-factor(damage$site, levels = c("K", "C", "W", "PH"))
damage$variety<-as.factor(damage$variety)
damage$treatment<-recode(damage$variety, !!!treatment_key)
damage$treatment<-factor(damage$treatment, levels = c("untreated", "treated"))

traits<-read.csv("clean_data/clean_traits_2019.csv")
traits$site<-fct_relevel(traits$site, "K", "C", "W", "PH")
traits$variety<-as.factor(traits$variety)
traits$treatment<-recode(traits$variety, !!!treatment_key)
traits$treatment<-droplevels(traits$treatment)

####graphing####
#boxplots
boxes<-ggplot(traits, aes(x = variety, y= Phi2))
boxes+geom_boxplot()
boxes+geom_boxplot()+facet_wrap(~site,nrow=1)

ggsave(filename = "figures/2019 phi2.tiff",width = 6,height = 4,units = "in")

#line graphs
avgInsects<-damage %>%
  group_by(date,treatment,site) %>%
  summarise(avg_insect_count=mean(insects),se_insect_count=sd(insects)/sqrt(n()))

lines<-ggplot(avgInsects, aes(x=date,y=avg_insect_count,colour=treatment,group=treatment))+facet_wrap(~site,nrow=1)
lines+geom_point()+geom_line()+geom_errorbar(aes(ymin=avg_insect_count-se_insect_count,ymax=avg_insect_count+se_insect_count))

avgChew<-damage %>%
  group_by(date,treatment,site) %>%
  summarise(avg_pct_chew=mean(chew_pct),se_pct_chew=sd(chew_pct)/sqrt(n()))

lines<-ggplot(avgChew, aes(x=date,y=avg_pct_chew,colour=treatment,group=treatment))+facet_wrap(~site,nrow=1)
lines+geom_point()+geom_line()+geom_errorbar(aes(ymin=avg_pct_chew-se_pct_chew,ymax=avg_pct_chew+se_pct_chew))

#can't figure out saving code for this style of figure- just export for now




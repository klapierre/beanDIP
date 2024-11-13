##############################################################################################
##
##  graph_non_damage_metrics.R: for any data except damage. can make line graphs with means  
##                              and SE, scatterplot matrices, and correlation heat maps
##
##  Author: Kelsey McGurrin
##
##############################################################################################


#### setup ####
library(tidyverse)
library(gridExtra)
library(ggcorrplot)

# working directory (add yours if different)
setwd("~/Dropbox/bean_dip_2018-2024/field trials/data")

# clean data file
data<-read_csv("clean_data/clean_all_years_long.csv")

#### all varieties, all years ####

# pick a column/metric to plot

siteAvgDf<-data %>%
  group_by(site,year,brandline,seed_treat) %>%
  summarise(site_avg=mean(beans_g,na.rm=T),
            se=sd(beans_g,na.rm = T)/sqrt(n()))

#not working now because of linetype?
# ggplot(siteAvgDf, aes(x=site,y=site_avg,group=seed_treat,linetype=seed_treat,color=brandline))+
#   geom_line(size=1.5)+
#   geom_point(size=2)+
#   geom_errorbar(aes(ymin=site_avg-se,
#                     ymax=site_avg+se),width=0.2,size=1.5)+
#   facet_wrap(~year,nrow=1)+
#   ylab("Site Average Beans (g)")


#### one variety across time ####
siteAvgDf<-data %>%
  filter(brandline_seedcoat=="AG38X8.untreated") %>%
  group_by(site,year,brandline_seedcoat) %>%
  summarise(site_avg=mean(beans_g,na.rm=T),
            se=sd(beans_g,na.rm = T)/sqrt(n()))

ggplot(siteAvgDf,aes(x=year,y=site_avg,group=site,color=site))+
  geom_point()+
  geom_line(size=1.5)+
  geom_errorbar(aes(ymin=site_avg-se,
                    ymax=site_avg+se),width=0.2,size=1.5)

#### treated/untreated pairs only ####
# AG38X8 2019-2021 
siteAvgDf<-data %>%
  filter(brandline=="AG38X8" & year== c(2019,2020,2021)) %>%
  group_by(site,year,brandline,seed_treat) %>%
  summarise(site_avg=mean(beans_g,na.rm=T),
            se=sd(beans_g,na.rm = T)/sqrt(n()))

ggplot(siteAvgDf,aes(x=site,y=site_avg,group=seed_treat,linetype=seed_treat,color=brandline))+
         geom_line(size=1.5)+
         geom_point()+
         geom_errorbar(aes(ymin=site_avg-se,
                           ymax=site_avg+se),width=0.2,size=1.5)+
         facet_wrap(~year,nrow=1)

# SH3814 2021
siteAvgDf<-data %>%
  filter(brandline=="SH3814 LL" & year== 2021) %>%
  group_by(site,year,brandline,seed_treat) %>%
  summarise(site_avg=mean(beans_g,na.rm=T),
            se=sd(beans_g,na.rm = T)/sqrt(n()))

ggplot(siteAvgDf,aes(x=site,y=site_avg,group=seed_treat,linetype=seed_treat,color=brandline))+
  geom_line(size=1.5)+
  geom_point()+
  geom_errorbar(aes(ymin=site_avg-se,
                    ymax=site_avg+se),width=0.2,size=1.5)+
  facet_wrap(~year,nrow=1)

####2021####
#import data
all<-read.csv("clean_data/clean_all_2021.csv")
all$variety<-as.factor(all$variety)
all$site<-factor(all$site, levels = c("K", "C", "W", "PH"))

#boxplots
boxes<-ggplot(all, aes(x = variety, y= beans_g))
boxes+geom_boxplot()+facet_wrap(~site,nrow = 1)

#scatterplots
sp<-ggplot(all, aes(x = plant_biomass_g, y= trial_yield))
sp+geom_point()+stat_smooth(method="lm")
#color code by variety
csp<-ggplot(all, aes(x = plant_biomass_g, y= trial_yield, color= site))
csp+geom_point()

#correlation of metrics
#select only metrics of interest
subsall<-select(all,
                plant_biomass_g,beans_g,healthy,
                root_biomass_g,
                nodule_count,stem_diameter_cm=root_collar_diameter_cm,
                avg_toughness,avg_trichomes,SLA,leaf_dry_matter,pressure_chamber,
                gH.,Phi2,PhiNPQ,Relative.Chlorophyll,Thickness,trial_yield
)

#heat map with correlation values
corrdata <- round(cor(subsall,use="complete.obs"), 1)
ggcorrplot(corrdata,type="lower",show.legend=F,lab=T,hc.order = T)

#scatterplot matrices
#trait data measured by hand
pairs(~SLA+leaf_dry_matter+avg_toughness+avg_trichomes+pressure_chamber,data=subsall)
#trait data from photosynq
pairs(~PhiNPQ+Phi2+Relative.Chlorophyll+gH.+Thickness,data=subsall)
#harvest data
pairs(~plant_biomass_g+beans_g+healthy+root_biomass_g+nodule_count+stem_diameter_cm,data=subsall)

#line graphs across environment
#calculate means and error bars
siteAvgDf<-all %>%
  group_by(site,variety) %>%
  summarise(site_avg_trichomes=mean(avg_trichomes,na.rm=T),
            se_avg_trichomes=sd(avg_trichomes,na.rm = T)/sqrt(n()),
            site_Relative.Chlorophyll=mean(Relative.Chlorophyll,na.rm=T),
            se_Relative.Chlorophyll=sd(Relative.Chlorophyll,na.rm=T)/sqrt(n()),
            site_Phi2=mean(Phi2,na.rm=T),
            se_Phi2=sd(Phi2,na.rm=T)/sqrt(n()),
            site_SLA=mean(SLA,na.rm=T),
            se_SLA=sd(SLA,na.rm = T)/sqrt(n()),
            site_beans_g=mean(beans_g,na.rm=T),
            se_beans_g=sd(beans_g,na.rm = T)/sqrt(n()),
            site_plant_biomass_g=mean(plant_biomass_g,na.rm=T),
            se_plant_biomass_g=sd(plant_biomass_g,na.rm = T)/sqrt(n()),
            site_root_biomass_g=mean(root_biomass_g,na.rm=T),
            se_root_biomass_g=sd(root_biomass_g,na.rm = T)/sqrt(n()),
            site_RtoS=mean(root_biomass_g/plant_biomass_g,na.rm=T),
            se_RtoS=sd((root_biomass_g/plant_biomass_g),na.rm=T/sqrt(n())),
            site_trial_yield=mean(trial_yield,na.rm=T),
            se_trial_yield=sd(trial_yield,na.rm = T)/sqrt(n())
  )
#graph with individual facets
facets<-ggplot(siteAvgDf, aes(x=site,y=site_trial_yield,colour=variety,group=variety))+facet_wrap(~variety,nrow=1)+geom_point()+geom_line()+
  geom_errorbar(aes(ymin=site_trial_yield-se_trial_yield,
                    ymax=site_trial_yield+se_trial_yield),width=0.2)
facets
#overlayed all together
overlay<-ggplot(siteAvgDf, aes(x=site,y=site_trial_yield,colour=variety,group=variety))+geom_line(size=1.5)+
  geom_errorbar(aes(ymin=site_trial_yield-se_trial_yield,
                    ymax=site_trial_yield+se_trial_yield),width=0.2,size=1.5)
overlay
#plot both together to make output easier
lay <- rbind(c(1,1,1,1),
             c(1,1,1,1),
             c(2,2,2,2))
grid.arrange(overlay, facets, layout_matrix=lay)
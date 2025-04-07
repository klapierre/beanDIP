##############################################################################################
##
##  graph_non_damage_metrics.R: for any data except damage. can make line graphs with means  
##                              and SE, scatterplot matrices, and correlation heat maps.
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

# open clean data file 
data<-read_csv("clean_data/clean_all_years_long.csv")


# filter so only varieties/years where we have a treated/untreated comparison
janitor::tabyl(data,year,seed_treat)

data %>% janitor::tabyl(year,brandline_seedcoat)

# should be 218 rows
data <- data %>%
  filter(brandline=="AG38X8" & year== c(2019,2020,2021) #122
         | brandline=="SH3814 LL" & year== 2021 #96
         ) 

# so sites display properly from west to east
data$site<-factor(data$site, levels = c("K", "C", "W", "PH"))

#### treated/untreated pairs only ####
# AG38X8 2019-2021 
siteAvgDf<-data %>%
  filter(brandline=="AG38X8") %>%
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
  filter(brandline=="SH3814 LL") %>%
  group_by(site,year,brandline,seed_treat) %>%
  summarise(site_avg=mean(beans_g,na.rm=T),
            se=sd(beans_g,na.rm = T)/sqrt(n()))

ggplot(siteAvgDf,aes(x=site,y=site_avg,group=seed_treat,linetype=seed_treat,color=brandline))+
  geom_line(size=1.5)+
  geom_point()+
  geom_errorbar(aes(ymin=site_avg-se,
                    ymax=site_avg+se),width=0.2,size=1.5)+
  facet_wrap(~year,nrow=1)

#### boxplots ####
boxes<-ggplot(data, aes(x = brandline_seedcoat, y= beans_g))
boxes+geom_boxplot()+facet_wrap(~site,nrow = 1)

#### scatterplots ####
sp<-ggplot(data, aes(x = plant_biomass_g, y= est_trial_yield))
sp+geom_point()+stat_smooth(method="lm")


#### correlation of metrics ####
#select only metrics of interest
subsdata<-select(data,
                plant_biomass_g,beans_g,healthy_count,
                root_biomass_g,
                nodule_count,stem_diameter_cm,
                leaf_toughness,leaf_trichomes_count,SLA,leaf_dry_matter,lwp_bar,
                Phi2_prop,PhiNPQ_prop,Relative.Chlorophyll,Thickness_mm,est_trial_yield
)

#heat map with correlation values
corrdata <- round(cor(subsdata,use="complete.obs"), 1)
ggcorrplot(corrdata,type="lower",show.legend=F,lab=T,hc.order = T)

#scatterplot matrices
#trait data measured by hand
pairs(~SLA+leaf_dry_matter+leaf_toughness+leaf_trichomes_count+lwp_bar,data=subsdata)
#trait data from photosynq
pairs(~PhiNPQ_prop+Phi2_prop+Relative.Chlorophyll+Thickness_mm,data=subsdata)
#harvest data
pairs(~plant_biomass_g+beans_g+healthy_count+root_biomass_g+nodule_count+stem_diameter_cm,data=subsdata)

#### line graphs across environment ####
#calculate means and error bars
siteAvgDf<-data %>%
  group_by(site,brandline_seedcoat) %>%
  summarise(site_avg_trichomes=mean(leaf_trichomes_count,na.rm=T),
            se_avg_trichomes=sd(leaf_trichomes_count,na.rm = T)/sqrt(n()),
            site_Relative.Chlorophyll=mean(Relative.Chlorophyll,na.rm=T),
            se_Relative.Chlorophyll=sd(Relative.Chlorophyll,na.rm=T)/sqrt(n()),
            site_Phi2=mean(Phi2_prop,na.rm=T),
            se_Phi2=sd(Phi2_prop,na.rm=T)/sqrt(n()),
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
            site_trial_yield=mean(est_trial_yield,na.rm=T),
            se_trial_yield=sd(est_trial_yield,na.rm = T)/sqrt(n())
  )
#graph with individual facets
facets<-ggplot(siteAvgDf, aes(x=site,y=site_trial_yield,colour=brandline_seedcoat,group=brandline_seedcoat))+
  facet_wrap(~brandline_seedcoat,nrow=1)+geom_point()+geom_line()+
  geom_errorbar(aes(ymin=site_trial_yield-se_trial_yield,
                    ymax=site_trial_yield+se_trial_yield),width=0.2)
facets
#overlayed all together
overlay<-ggplot(siteAvgDf, aes(x=site,y=site_trial_yield,colour=brandline_seedcoat,group=brandline_seedcoat))+geom_line(size=1.5)+
  geom_errorbar(aes(ymin=site_trial_yield-se_trial_yield,
                    ymax=site_trial_yield+se_trial_yield),width=0.2,size=1.5)
overlay
#plot both together to make output easier
lay <- rbind(c(1,1,1,1),
             c(1,1,1,1),
             c(2,2,2,2))
grid.arrange(overlay, facets, layout_matrix=lay)

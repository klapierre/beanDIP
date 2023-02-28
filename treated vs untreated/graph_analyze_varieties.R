####setup####
library(tidyverse)
library(lme4)
library(lmerTest)
library(gridExtra)
library(ggcorrplot)
setwd("~/Dropbox/bean_dip_2018-2024/field trials/data")

####2019####
#import data
all<-read.csv("clean_data/clean_all_2019.csv")
all$variety<-as.factor(all$variety)
all$site<-factor(all$site, levels = c("K", "C", "W", "PH"))

#compare overall trial yield to our plants
#per plot calculations
per_plot<-all %>%
  group_by(site,variety,plot_id,trial_yield) %>%
  summarise(plot_density=mean(density),
            plot_beans_g=mean(beans_g),
            plot_yield=((((plot_density/8)*plot_beans_g)/454)/60)*43560)
#graph both plot level yield metrics together
p1<-ggplot(per_plot,aes(x=variety,y=trial_yield))+geom_boxplot()+facet_wrap(~site,nrow=1)+coord_cartesian(ylim=c(25,100))
p2<-ggplot(per_plot,aes(x=variety,y=plot_yield))+geom_boxplot()+facet_wrap(~site,nrow=1)+coord_cartesian(ylim=c(25,100))
grid.arrange(p1,p2,nrow=2)

#boxplots
boxes<-ggplot(all, aes(x = variety, y= Thickness))
boxes+geom_boxplot()+facet_wrap(~site,nrow = 1)#+coord_cartesian(ylim=c(0,400))
#ggsave(filename = "figures/2019 phi2.tiff",width = 6,height = 4,units = "in")

#modeling
#probably best model from 1st meeting but rank deficient & singular
mod4<-lmer(avg_trichomes~variety*site+(1|site:plot:variety), data=all)
#dropping to one fixed effect solves singularity error
mod4<-lmer(avg_trichomes~variety+(1|site:plot:variety), data=all)
mod4<-lmer(avg_trichomes~site+(1|site:plot:variety), data=all)
#trying new organization still singularity
mod5<-lmer(avg_trichomes~variety*site+(1|plot_id), data=all)
#model results
anova(mod5)
summary(mod5)
ranef(mod5)
#check model assumptions
plot(mod5)
hist(resid(mod5))

#scatterplots
sp<-ggplot(all, aes(x = Thickness, y= avg_trichomes))
sp+geom_point()+stat_smooth(method="lm")
#color code by variety
csp<-ggplot(all, aes(x = beans_g, y= trial_yield, color= site))
csp+geom_point()

#correlation of metrics
#select only metrics of interest
subsall<-select(all,
                plant_biomass_g,beans_g,healthy,trial_yield,root_biomass_g,functional_nodules,stem_diameter_cm,
                avg_toughness,avg_trichomes,SLA,leaf_dry_matter,osmom_final,pressure_chamber,
                gH.,Phi2,PhiNPQ,Relative.Chlorophyll,Thickness,trial_yield
)

#heat map with correlation values
corrdata <- round(cor(subsall,use="complete.obs"), 1)
ggcorrplot(corrdata,type="lower",show.legend=F,lab=T,hc.order = T)

#scatterplot matrices
#all data is too much
#pairs(~.,data=subsall)
#trait data measured by hand
pairs(~SLA+leaf_dry_matter+avg_toughness+avg_trichomes+pressure_chamber,data=subsall)
#color by variety
pairs(~SLA+leaf_dry_matter+avg_toughness+avg_trichomes+pressure_chamber,data=subsall)
#trait data from photosynq
pairs(~PhiNPQ+Phi2+Relative.Chlorophyll+gH.,data=subsall)
#harvest data
pairs(~plant_biomass_g+beans_g+healthy+trial_yield+root_biomass_g+functional_nodules+stem_diameter_cm,data=subsall)
#alternate colored matrix
# library(car)
# scatterplotMatrix(~avg_trichomes+Phi2+Relative.Chlorophyll+beans_g|variety,data=all)

#line graphs across time
#calculate means and error bars
damage <- read_csv("clean_data/clean_damage_2019.csv")
avgInsects<-damage %>%
  group_by(date,variety,site) %>%
  summarise(avg_insect_count=mean(insects),se_insect_count=sd(insects)/sqrt(n()))
avgInsects$variety<-as.factor(avgInsects$variety)
#graph
lines<-ggplot(avgInsects, aes(x=date,y=avg_insect_count,colour=variety,group=variety))+facet_wrap(~site,nrow=1)
lines+geom_point()+geom_line()+geom_errorbar(aes(ymin=avg_insect_count-se_insect_count,ymax=avg_insect_count+se_insect_count))
#calculate means and error bars
avgChew<-damage %>%
  group_by(date,variety,site) %>%
  summarise(avg_pct_chew=mean(chew_pct),se_pct_chew=sd(chew_pct)/sqrt(n()))
avgChew$variety<-as.factor(avgChew$variety)
#graph
lines<-ggplot(avgChew, aes(x=date,y=avg_pct_chew,colour=variety,group=variety))+facet_wrap(~site,nrow=1)
lines+geom_point()+geom_line()+geom_errorbar(aes(ymin=avg_pct_chew-se_pct_chew,ymax=avg_pct_chew+se_pct_chew))

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
            site_trial_yield=mean(trial_yield,na.rm=T),
            se_trial_yield=sd(trial_yield,na.rm = T)/sqrt(n()))
#graph with individual facets
facets<-ggplot(siteAvgDf, aes(x=site,y=site_beans_g,colour=variety,group=variety))+facet_wrap(~variety,nrow=1)+geom_point()+geom_line()+
  geom_errorbar(aes(ymin=site_beans_g-se_beans_g,
                    ymax=site_beans_g+se_beans_g),width=0.2)
facets
#overlayed all together
overlay<-ggplot(siteAvgDf, aes(x=site,y=site_beans_g,colour=variety,group=variety))+geom_line(size=1.5)+
  geom_errorbar(aes(ymin=site_beans_g-se_beans_g,
                    ymax=site_beans_g+se_beans_g),width=0.2,size=1.5)
overlay
#plot both to make output easier
lay <- rbind(c(1,1,1,1),
             c(1,1,1,1),
             c(2,2,2,2))
grid.arrange(overlay, facets, layout_matrix=lay)

####2020####
#import data
all<-read.csv("clean_data/clean_all_2020.csv")
all$variety<-as.factor(all$variety)
all$site<-factor(all$site, levels = c("K", "C", "W", "PH"))

#boxplots
boxes<-ggplot(all, aes(x = variety, y= root_biomass_g/plant_biomass_g))
boxes+geom_boxplot()+facet_wrap(~site,nrow = 1)

#scatterplots
sp<-ggplot(all, aes(x = beans_g, y= trial_yield))
sp+geom_point()+stat_smooth(method="lm")
#color code by variety
csp<-ggplot(all, aes(x = beans_g, y= trial_yield, color= site))
csp+geom_point()

#correlation of metrics
#select only metrics of interest
subsall<-select(all,
                plant_biomass_g,beans_g,healthy,root_biomass_g,nodule_count,stem_diameter_cm,
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
                #root_biomass_g,
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
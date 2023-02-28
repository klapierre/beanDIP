#######################################################################################
##  combine_data.R: Takes individual metrics for each year, calculates more,
##                  removes outliers (when necessary),
##                  and saves the data into a single file for each year.
##                  Then combines files for each year into a single one
##                  with "long" format
##
##  Author: Kelsey McGurrin
##
#######################################################################################

####setup####
library(tidyverse)
site_key<-c("K", "C", "W", "PH")

# working directory path (add yours if different)
setwd("~/Dropbox/bean_dip_2018-2024/field trials/data")

####2019####
#import data
harvest<-read.csv("clean_data/clean_harvest_2019.csv")
damage<-read.csv("clean_data/clean_damage_2019.csv")
dens<-read_csv("clean_data/clean_density_2019.csv")
traits<-read.csv("clean_data/clean_traits_2019.csv")
trial<-read.csv("clean_data/clean_trial_harvest_2019.csv")

#remove outliers
#one nest of baby stinkbugs from 8/28 Wye
damage<-subset(damage,insects<55)
#one really low gH point from PH 67
traits<-subset(traits,gH. > -500)
#two really high leaf area points from W 
traits<-subset(traits,area_cm < 400)
#one really thick leaf from PH
traits<-subset(traits,Thickness < 2)

#create single table with trait and harvest data
all<-left_join(harvest,traits,by=c("site", "variety","plot","indiv"))
all<-left_join(all,dens,by=c("site", "variety","plot","indiv","row"))
all<-left_join(all,trial,by=c("site", "variety","plot"))
all$variety<-as.factor(all$variety)
all$site<-fct_relevel(all$site, site_key)

#compute more useful metrics
all<-mutate(all,
             avg_toughness=(toughness_1+toughness_2)/2,
             avg_trichomes=(trichomes_top+trichomes_bottom)/2,
             SLA=area_cm/dry_mass,
             leaf_dry_matter=dry_mass/wet_mass,
            nodules_per_g=functional_nodules/root_biomass_g,
            beans_per_g=beans_g/plant_biomass_g
             )

#create unique IDs to make it easier to group random variables
all$plot_id<-interaction(all$site,all$variety,all$plot)
all$block_id<-interaction(all$site,all$plot)

#save all data file for easier input with later analyses
write.csv(all,file="clean_data/clean_all_2019.csv",row.names=F)

####2020####
#import data
harvest<-read.csv("clean_data/clean_harvest_2020.csv")
damage<-read.csv("clean_data/clean_damage_2020.csv")
dens<-read.csv("clean_data/clean_density_2020.csv")
traits<-read.csv("clean_data/clean_traits_2020.csv")
trial<-read.csv("clean_data/clean_trial_harvest_2020.csv")

#remove outliers
#one really high gH value from PH 2-2-2
traits<-subset(traits,gH. < 1500)

#one really small plant from PH 2-3-2
harvest<-subset(harvest,plant_biomass_g > 1)

#create single table with trait and harvest data
all<-left_join(harvest,traits,by=c("site", "variety","plot","indiv"))
all<-left_join(all,dens,by=c("site", "variety","plot","row"))
all<-left_join(all,trial,by=c("site", "variety","plot"))
all$variety<-as.factor(all$variety)
all$site<-fct_relevel(all$site, site_key)

#compute more useful metrics
all<-mutate(all,
            avg_toughness=(toughness_1+toughness_2)/2,
            avg_trichomes=(trichomes_top+trichomes_bottom)/2,
            SLA=area_cm/dry_mass,
            leaf_dry_matter=dry_mass/wet_mass,
            nodules_per_g=nodule_count/root_biomass_g,
            beans_per_g=beans_g/plant_biomass_g,
            plants_per_foot=final/4)

#create unique IDs to make it easier to group random variables
all$plot_id<-interaction(all$site,all$variety,all$plot)
all$block_id<-interaction(all$site,all$plot)

#save all data file for easier input with later analyses
write.csv(all,file="clean_data/clean_all_2020.csv",row.names=F)

####2021####
#import data
harvest<-read.csv("clean_data/clean_harvest_2021.csv")
damage<-read.csv("clean_data/clean_damage_2021.csv")
dens<-read.csv("clean_data/clean_density_2021.csv")
traits<-read.csv("clean_data/clean_traits_2021.csv")
trial<-read.csv("clean_data/clean_trial_harvest_2021.csv")

#remove outliers (none known)

#create single table with trait and harvest data
all<-left_join(harvest,traits,by=c("site", "variety","plot","indiv"))
all<-left_join(all,dens,by=c("site", "variety","plot","row"))
all<-left_join(all,trial,by=c("site", "variety","plot"))
all$variety<-as.factor(all$variety)
all$site<-fct_relevel(all$site, site_key)

#compute more useful metrics
all<-mutate(all,
            area_cm=(leaf_area_cm_1+leaf_area_cm_2)/2,
            avg_toughness=(toughness_1+toughness_2)/2,
            avg_trichomes=(trichomes_top+trichomes_bottom)/2,
            SLA=area_cm/dry_mass,
            leaf_dry_matter=dry_mass/wet_mass,
            nodules_per_g=nodule_count/root_biomass_g,
            beans_per_g=beans_g/plant_biomass_g,
            plants_per_foot=final/4)

#create unique IDs to make it easier to group random variables
all$plot_id<-interaction(all$site,all$variety,all$plot)
all$block_id<-interaction(all$site,all$plot)

#save all data file for easier input with later analyses
write.csv(all,file="clean_data/clean_all_2021.csv",row.names=F)

####2022####
#import data
harvest<-read.csv("clean_data/clean_harvest_2022.csv")
damage<-read.csv("clean_data/clean_damage_2022.csv")
dens<-read.csv("clean_data/clean_density_2022.csv")
traits<-read.csv("clean_data/clean_traits_2022.csv")
trial<-read.csv("clean_data/clean_trial_harvest_2022.csv")

#remove outliers (none known)

#create single table with trait and harvest data
all<-left_join(harvest,traits,by=c("site", "variety","plot","indiv"))
all<-left_join(all,dens,by=c("site", "variety","plot","row"))
all<-left_join(all,trial,by=c("site", "variety","plot"))
all$variety<-as.factor(all$variety)
all$site<-fct_relevel(all$site, site_key)

#compute more useful metrics
all<-mutate(all,
            area_cm=(leaf_area_cm_1+leaf_area_cm_2)/2,
            avg_toughness=(toughness_1+toughness_2)/2,
            avg_trichomes=(trichomes_top+trichomes_bottom)/2,
            SLA=area_cm/dry_mass,
            leaf_dry_matter=dry_mass/wet_mass,
            nodules_per_g=nodule_count/root_biomass_g,
            beans_per_g=beans_g/plant_biomass_g,
            plants_per_foot=final/4)

#create unique IDs to make it easier to group random variables
all$plot_id<-interaction(all$site,all$variety,all$plot)
all$block_id<-interaction(all$site,all$plot)

#save all data file for easier input with later analyses
write.csv(all,file="clean_data/clean_all_2022.csv",row.names=F)

####combine years####
dat_19<-read.csv("clean_data/clean_all_2019.csv")
dat_20<-read.csv("clean_data/clean_all_2020.csv")
dat_21<-read.csv("clean_data/clean_all_2021.csv")
dat_22<-read.csv("clean_data/clean_all_2022.csv")

#add in brand names provided to SVT
dat_19$variety<-as.factor(dat_19$variety)
dat_19_brand_key <- c("27"="AG38X8","67"="AG38X8", "31"="S39-G2X", "83"="SH3814 LL",
                 "32"="S39XT68","55"="7390ET")
dat_19$brandline<-recode(dat_19$variety, !!!dat_19_brand_key)
#add in groupings for treated or untreated seed coats
dat_19_treat_key <- c("27"="treated", "67"="untreated", "31"= "treated", "83"="treated"
                 ,"32"="treated","55"="treated")
dat_19$seed_treat<-recode(dat_19$variety, !!!dat_19_treat_key)
dat_19$brandline_seedcoat<-interaction(dat_19$brandline,dat_19$seed_treat)

#clean up names and metrics
dat_19<-mutate(dat_19,
          nodule_count=functional_nodules+scenesced_nodules,
          lwp_bar=pressure_chamber*10)

dat_19<-select(dat_19,site,plot,indiv,
          seed_treat,brandline,brandline_seedcoat,
          row,pheno_stage,
          plant_biomass_g,beans_g,row_dens=density,plot_yield_buac=trial_yield,
          root_biomass_g,nodule_count,stem_diameter_cm,
          wafer_count=wafer,discolored_count=purple,
          discolored_wrinkled_count=purple_wrinkled,
          wrinkled_count=wrinkled,healthy_count=healthy,
          leaf_area_cm=area_cm,wet_leaf_g=wet_mass,lwp_bar,dry_leaf_g=dry_mass,
          leaf_toughness=avg_toughness,leaf_trichomes_count=avg_trichomes,
          SLA,leaf_dry_matter,nodules_per_g,beans_per_g,
          humid_pct=Ambient.Humidity,temp_C=Ambient.Temperature,
          leaf_angle_deg=Leaf.Angle,
          leaf_temp_diff_C=Leaf.Temperature.Differential,
          PAR_uE=Light.Intensity..PAR.,
          Phi2_prop=Phi2,PhiNPQ_prop=PhiNPQ,
          Relative.Chlorophyll,Thickness_mm=Thickness
)

#same details for 2020 data
dat_20$variety<-as.factor(dat_20$variety)
dat_20_brand_key <- c("2"="AG38X8","5"="AG38X8", "70"="S39-G2X", "50"="SH3814 LL")
dat_20$brandline<-recode(dat_20$variety, !!!dat_20_brand_key)
dat_20_treat_key <- c("5"="treated", "2"="untreated", "50"= "treated", "70"="treated")
dat_20$seed_treat<-recode(dat_20$variety, !!!dat_20_treat_key)
dat_20$brandline_seedcoat<-interaction(dat_20$brandline,dat_20$seed_treat)
#also add in plot sizes for 2020
dat_20_length_key<- c("PH"=18,"W"=18,"C"=18,"K"=13)
dat_20$plot_length<-recode(dat_20$site,!!!dat_20_length_key)
#clean up names and metrics
dat_20<-mutate(dat_20,planted_length=plot_length*2*3,
          est_num_plants=planted_length*plants_per_foot,
          est_trial_yield=est_num_plants*beans_g)
dat_20<-select(dat_20,site,plot,indiv,
          seed_treat,brandline,brandline_seedcoat,
          row,pheno_stage,
          plant_biomass_g,beans_g,row_dens=initial,
          root_biomass_g,nodule_count,stem_diameter_cm,
          wafer_count=wafer,discolored_count=purple,
          discolored_wrinkled_count=purple_wrinkled,
          wrinkled_count=wrinkled,healthy_count=healthy,
          leaf_area_cm=area_cm,wet_leaf_g=wet_mass,lwp_bar=pressure_chamber,
          dry_leaf_g=dry_mass,
          leaf_toughness=avg_toughness,leaf_trichomes_count=avg_trichomes,
          SLA,leaf_dry_matter,nodules_per_g,beans_per_g,
          humid_pct=Ambient.Humidity,temp_C=Ambient.Temperature,
          leaf_angle_deg=Leaf.Angle,
          leaf_temp_diff_C=Leaf.Temperature.Differential,
          PAR_uE=Light.Intensity..PAR.,
          Phi2_prop=Phi2,PhiNPQ_prop=PhiNPQ,
          Relative.Chlorophyll,Thickness_mm=Thickness,
          plot_yield_buac=trial_yield,est_trial_yield
)

#same details for 2021 data
dat_21$variety<-as.factor(dat_21$variety)
dat_21_brand_key <- c("19"="AG38X8","59"="AG38X8", "82"="S39-G2X", "57"="SH3814 LL", "58"="SH3814 LL")
dat_21$brandline<-recode(dat_21$variety, !!!dat_21_brand_key)
dat_21_treat_key <- c("59"="treated", "19"="untreated", "82"= "treated", "58"="treated", "57"="untreated")
dat_21$seed_treat<-recode(dat_21$variety, !!!dat_21_treat_key)
dat_21$brandline_seedcoat<-interaction(dat_21$brandline,dat_21$seed_treat)
#also add in plot sizes
dat_21_length_key<- c("PH"=18,"W"=13,"C"=13,"K"=13)
dat_21$plot_length<-recode(dat_21$site,!!!dat_21_length_key)
#clean up names and metrics
dat_21<-mutate(dat_21,planted_length=plot_length*2*3,
          est_num_plants=planted_length*plants_per_foot,
          est_trial_yield=est_num_plants*beans_g)
dat_21<-select(dat_21,site,plot,indiv,
          seed_treat,brandline,brandline_seedcoat,
          row,
          pheno_stage,
          plant_biomass_g,beans_g,row_dens=initial,
          root_biomass_g,
          nodule_count,stem_diameter_cm=root_collar_diameter_cm,
          wafer_count=wafer,discolored_count=purple,
          discolored_wrinkled_count=purple_wrinkled,
          wrinkled_count=wrinkled,healthy_count=healthy,
          leaf_area_cm=area_cm,wet_leaf_g=wet_mass,lwp_bar=pressure_chamber,
          dry_leaf_g=dry_mass,
          leaf_toughness=avg_toughness,leaf_trichomes_count=avg_trichomes,
          SLA,leaf_dry_matter,
          nodules_per_g,
          beans_per_g,
          humid_pct=Ambient.Humidity,temp_C=Ambient.Temperature,
          leaf_angle_deg=Leaf.Angle,
          leaf_temp_diff_C=Leaf.Temperature.Differential,
          PAR_uE=Light.Intensity..PAR.,
          Phi2_prop=Phi2,PhiNPQ_prop=PhiNPQ,
          Relative.Chlorophyll,Thickness_mm=Thickness,
          plot_yield_buac=trial_yield,est_trial_yield
)

#same details for 2022 data
dat_22$variety<-as.factor(dat_22$variety)
dat_22_brand_key <- c("73"="AG38X8")
dat_22$brandline<-recode(dat_22$variety, !!!dat_22_brand_key)
dat_22_treat_key <- c("73"="untreated")
dat_22$seed_treat<-recode(dat_22$variety, !!!dat_22_treat_key)
dat_22$brandline_seedcoat<-interaction(dat_22$brandline,dat_22$seed_treat)
#also add in plot sizes
dat_22_length_key<- c("PH"=18,"W"=13,"C"=13,"K"=13)
dat_22$plot_length<-recode(dat_22$site,!!!dat_22_length_key)
#clean up names and metrics
dat_22<-mutate(dat_22,planted_length=plot_length*2*3,
          est_num_plants=planted_length*plants_per_foot,
          est_trial_yield=est_num_plants*beans_g)
dat_22<-select(dat_22,site,plot,indiv,
          seed_treat,brandline,brandline_seedcoat,
          row,
          pheno_stage,
          plant_biomass_g,beans_g,row_dens=initial,
          root_biomass_g,
          nodule_count,stem_diameter_cm=root_collar_diameter_cm,
          wafer_count=wafer,discolored_count=purple,
          discolored_wrinkled_count=purple_wrinkled,
          wrinkled_count=wrinkled,healthy_count=healthy,
          leaf_area_cm=area_cm,wet_leaf_g=wet_mass,lwp_bar=pressure_chamber,
          dry_leaf_g=dry_mass,
          leaf_toughness=avg_toughness,leaf_trichomes_count=avg_trichomes,
          SLA,leaf_dry_matter,
          nodules_per_g,
          beans_per_g,
          humid_pct=Ambient.Humidity,temp_C=Ambient.Temperature,
          leaf_angle_deg=Leaf.Angle,
          leaf_temp_diff_C=Leaf.Temperature.Differential,
          PAR_uE=Light.Intensity..PAR.,
          Phi2_prop=Phi2,PhiNPQ_prop=PhiNPQ,
          Relative.Chlorophyll,Thickness_mm=Thickness,
          plot_yield_buac=trial_yield,est_trial_yield
)

#long join data
dat_19<-add_column(dat_19,.before="site",year=as.factor(2019))
dat_20<-add_column(dat_20,.before="site",year=as.factor(2020))
dat_21<-add_column(dat_21,.before="site",year=as.factor(2021))
dat_22<-add_column(dat_22,.before="site",year=as.factor(2022))
long<-bind_rows(dat_19,dat_20,dat_21,dat_22)
write.csv(long,file="clean_data/clean_all_years_long.csv",row.names=F)







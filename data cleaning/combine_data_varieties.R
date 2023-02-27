####setup####
library(tidyverse)
setwd("~/Dropbox/bean_dip_2018-2024/field trials/data")
site_key<-c("K", "C", "W", "PH")

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

#remove outliers

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
            #nodules_per_g=nodule_count/root_biomass_g,
            beans_per_g=beans_g/plant_biomass_g,
            plants_per_foot=final/4)

#create unique IDs to make it easier to group random variables
all$plot_id<-interaction(all$site,all$variety,all$plot)
all$block_id<-interaction(all$site,all$plot)

#save all data file for easier input with later analyses
write.csv(all,file="clean_data/clean_all_2021.csv",row.names=F)

####combine years####
a<-read.csv("clean_data/clean_all_2019.csv")
b<-read.csv("clean_data/clean_all_2020.csv")
c<-read.csv("clean_data/clean_all_2021.csv")

#select varieties from 2019 which were measured more than one year
a$variety<-as.factor(a$variety)
#a<-filter(a,variety %in% c("27","67","31","83"))
#droplevels(a$variety)

#add in brand names provided to Jason Wight
a_brand_key <- c("27"="AG38X8","67"="AG38X8", "31"="S39-G2X", "83"="SH3814 LL",
                 "32"="S39XT68","55"="7390ET")
a$brandline<-recode(a$variety, !!!a_brand_key)
#droplevels(a$brandline)
#add in groupings for treated or untreated seed coats
a_treat_key <- c("27"="treated", "67"="untreated", "31"= "treated", "83"="treated"
                 ,"32"="treated","55"="treated")
a$seed_treat<-recode(a$variety, !!!a_treat_key)
a$brandline_seedcoat<-interaction(a$brandline,a$seed_treat)
#clean up names and metrics
a<-mutate(a,
          nodule_count=functional_nodules+scenesced_nodules,
          lwp_bar=pressure_chamber*10)

a<-select(a,site,plot,indiv,
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
b$variety<-as.factor(b$variety)
b_brand_key <- c("2"="AG38X8","5"="AG38X8", "70"="S39-G2X", "50"="SH3814 LL")
b$brandline<-recode(b$variety, !!!b_brand_key)
b_treat_key <- c("5"="treated", "2"="untreated", "50"= "treated", "70"="treated")
b$seed_treat<-recode(b$variety, !!!b_treat_key)
b$brandline_seedcoat<-interaction(b$brandline,b$seed_treat)
#also add in plot sizes for 2020
b_length_key<- c("PH"=18,"W"=18,"C"=18,"K"=13)
b$plot_length<-recode(b$site,!!!b_length_key)
#clean up names and metrics
b<-mutate(b,planted_length=plot_length*2*3,
          est_num_plants=planted_length*plants_per_foot,
          est_trial_yield=est_num_plants*beans_g)
b<-select(b,site,plot,indiv,
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
c$variety<-as.factor(c$variety)
c_brand_key <- c("19"="AG38X8","59"="AG38X8", "82"="S39-G2X", "57"="SH3814 LL", "58"="SH3814 LL")
c$brandline<-recode(c$variety, !!!c_brand_key)
c_treat_key <- c("59"="treated", "19"="untreated", "82"= "treated", "58"="treated", "57"="untreated")
c$seed_treat<-recode(c$variety, !!!c_treat_key)
c$brandline_seedcoat<-interaction(c$brandline,c$seed_treat)
#also add in plot sizes
c_length_key<- c("PH"=18,"W"=13,"C"=13,"K"=13)
c$plot_length<-recode(c$site,!!!c_length_key)
#clean up names and metrics
c<-mutate(c,planted_length=plot_length*2*3,
          est_num_plants=planted_length*plants_per_foot,
          est_trial_yield=est_num_plants*beans_g)
c<-select(c,site,plot,indiv,
          seed_treat,brandline,brandline_seedcoat,
          row,
          pheno_stage,
          plant_biomass_g,beans_g,row_dens=initial,
          #root_biomass_g,
          nodule_count,stem_diameter_cm=root_collar_diameter_cm,
          wafer_count=wafer,discolored_count=purple,
          discolored_wrinkled_count=purple_wrinkled,
          wrinkled_count=wrinkled,healthy_count=healthy,
          leaf_area_cm=area_cm,wet_leaf_g=wet_mass,lwp_bar=pressure_chamber,
          dry_leaf_g=dry_mass,
          leaf_toughness=avg_toughness,leaf_trichomes_count=avg_trichomes,
          SLA,leaf_dry_matter,
          #nodules_per_g,
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
a<-add_column(a,.before="site",year=as.factor(2019))
b<-add_column(b,.before="site",year=as.factor(2020))
c<-add_column(c,.before="site",year=as.factor(2021))
long<-bind_rows(a,b,c)
write.csv(long,file="clean_data/clean_all_years_long.csv",row.names=F)

#wide join data
a<-select(a,-year)
b<-select(b,-year)
c<-select(c,-year)

ab<-left_join(b,a,suffix=c(".20",".19"),by=
                  c("site","plot","indiv","seed_treat","brandline",
                    "brandline_seedcoat","row"))

#paste ".21" onto end of all metrics before joining with other years
colnames(c)[8:38] <- paste(colnames(c[8:38]),"21",sep = ".")

wide<-left_join(c,ab,by=
                   c("site","plot","indiv","seed_treat","brandline",
                     "brandline_seedcoat","row"))

write.csv(wide,file="clean_data/clean_all_years_wide.csv",row.names=F)





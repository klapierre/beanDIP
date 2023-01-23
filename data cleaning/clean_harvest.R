################################################################################
##  clean_harvest.R: Cleaning and combining variety trial harvest metrics
##
##  Author: Kelsey McGurrin
################################################################################

####setup####
library(readxl)
library(tidyverse)

# working directory path (add yours if different)
setwd("~/Dropbox/bean_dip_2018-2024/field trials/data")

####2019####
#aboveground biomass, bean yield, bean counts
data <- read_excel("raw_data/2019_AGbiomass_counts.xlsx")
data$variety<-as.factor(data$variety)
data<-rename(data,plant_biomass_g=dry_g)

#add in yield data from overall variety trial
yield<-read_excel("raw_data/2019_statewide_yield.xlsx")
yield<-yield %>%
  select("Entry No.","KV full season bu/ac","CV full season bu/ac","W full season bu/ac","PH full season bu/ac") %>%
  rename(variety="Entry No.",K="KV full season bu/ac",C="CV full season bu/ac",W="W full season bu/ac",PH="PH full season bu/ac")
yield$K<-as.numeric(yield$K)
yield$PH<-as.numeric(yield$PH)
yield$variety<-as.factor(yield$variety)
yield<-pivot_longer(yield,cols = K:PH,names_to="site",values_to="trial_yield")
harvest<-left_join(data,yield,by=c("site", "variety"))

#add in root biomass data
#had same bag labelling issues as nodule data see Excel
roots<-read.csv("raw_data/beanDIP_2019_rootbiomass_corrected.csv")
roots<-rename(roots,variety=varietal, plot=block, indiv=plant)
level_key <- c(KD = "K", CV = "C")
roots$site<-recode(roots$site, !!!level_key)
roots$variety<-as.factor(roots$variety)
harvest<-left_join(harvest,roots,by=c("site", "variety","plot","indiv"))

#add in nodule data
#several corrections were needed to original .csv see Excel
nodules<-read.csv("raw_data/beanDIP_2019_nodules_corrected.csv")
nodules<-rename(nodules,variety=varietal, plot=block, indiv=plant)
nodules$site<-recode(nodules$site, !!!level_key)
nodules$variety<-as.factor(nodules$variety)
harvest$site<-as.factor(harvest$site)
harvest<-left_join(harvest,nodules,by=c("site", "variety","plot","indiv"))

#add in pheno stage at harvest so we can group by it
pheno<-read.csv("clean_data/clean_pheno_2019.csv")
pheno$variety<-as.factor(pheno$variety)
pheno<-filter(pheno, date %in% c("2019-09-16","2019-09-17"))
harvest<-left_join(harvest,pheno,by=c("site", "variety", "plot", "indiv"))

#remove extra columns
harvest<-select(harvest,-c(date.x,date.y,counter))

#output cleaned file, should have 204 rows total
write.csv(harvest,file="clean_data/clean_harvest_2019.csv",row.names = F)

####2020####
#aboveground biomass, bean yield, bean counts
data <- read_excel("raw_data/2020_AGbiomass_counts.xlsx")
data <- select(data,-date)
data<-rename(data,plant_biomass_g=dry_g)
#add in yield data from overall variety trial


#add in root biomass data
roots<-read_excel("raw_data/beanDIP 2020 root biomass weight data.xlsx")
roots<-rename(roots,site=Site,variety="Variety #", plot="Plot #", indiv="Plant #",root_biomass_g="Root biomass (g)")

#values weren't matching up, so audited data and reweighed some samples
roots<-arrange(roots,site,variety,plot,indiv)

#there were two C 2-1-1 recorded and no C 2-1-2
# C 2-1-2 .4341 g very close to original .437 so not replaced
roots$indiv[roots$root_biomass_g==.4379]<-2
# C 2-1-1   1.8855 g
roots$root_biomass_g[roots$site=="C" & roots$variety==2 & roots$plot==1 & roots$indiv==1]<-1.8855
#duplicate W 2-2-3 labels
# W 2-3-3 4.0427 g
# W 2-2-3 3.0547 g
#the W 2-2-3 that weighs ~ 4 needs to be changed to plot 3
roots$plot[roots$root_biomass_g==4.1040]<-3
#missing W 50-2-6
# W 50-2-6 2.2820 g
roots<-add_row(roots,site="W",variety=50,plot=2,indiv=6,root_biomass_g=2.2820)

#some labels entered incorrectly (indiv-plot instead of plot-indiv)
#C 50-4-1 should be 50-1-4
roots$plot[roots$root_biomass_g==1.7245]<-1
roots$indiv[roots$root_biomass_g==1.7245]<-4
#C 50-6-1 should be 50-1-6
roots$plot[roots$root_biomass_g==2.6944]<-1
roots$indiv[roots$root_biomass_g==2.6944]<-6
#C 70-4-1 should be 70-1-4
roots$plot[roots$root_biomass_g==1.1368]<-1
roots$indiv[roots$root_biomass_g==1.1368]<-4
#C 70-6-1 should be 70-1-6
roots$plot[roots$root_biomass_g==1.7141]<-1
roots$indiv[roots$root_biomass_g==1.7141]<-6
#W 50-6-2 should be 50-2-6 but already value 2.2820 for that indiv?
#similar weight so just eliminate row?
# roots$plot[roots$root_biomass_g==2.3676]<-2
# roots$indiv[roots$root_biomass_g==2.3676]<-6
roots[roots$root_biomass_g==2.3676,]<-NA
roots<-drop_na(roots)

#some issues from bag labeling mistakes
#There are two bags labeled PH 50-1-6 (and no 5-1-6) and two bags labeled C 50-3-1 (and no 50-1-3)
#missing C 70-1-3 but have duplicate C 70-3-1
#70-1-3 had larger aboveground biomass, so assign it larger root value
roots$plot[roots$root_biomass_g==2.7829]<-1
roots$indiv[roots$root_biomass_g==2.7829]<-3
#missing PH 5-1-6 but have duplicate PH 50-1-6
#5-1-6 had larger aboveground biomass, so assign it larger root value
roots$variety[roots$root_biomass_g==1.6433]<-5
#missing C 50-1-3 but have duplicate C 50-3-1
#50-1-3 had smaller aboveground biomass, so assign it smaller root value
roots$plot[roots$root_biomass_g==1.2912]<-1
roots$indiv[roots$root_biomass_g==1.2912]<-3

harvest<-left_join(data,roots,by=c("site", "variety", "plot", "indiv"))

#add in nodule data
nodules<-read_excel("raw_data/beanDIP 2020 nodule and root collar data.xlsx")
nodules<-rename(nodules,site=Site,variety="Variety #", plot="Plot #", indiv="Plant #",stem_diameter_cm="Root collar diameter (cm)",nodule_count="Number of nodules")

#more issues with rows not matching up
#harvest2<-anti_join(harvest,nodules,by=c("site", "variety", "plot", "indiv"))

#W 50-6-2 should be 50-2-6 
nodules$plot[nodules$stem_diameter_cm==0.6 & nodules$nodule_count==78]<-2
nodules$indiv[nodules$stem_diameter_cm==0.6 & nodules$nodule_count==78]<-6
#C 50-4-1 should be 50-1-4
nodules$plot[nodules$stem_diameter_cm==1 & nodules$nodule_count==27]<-1
nodules$indiv[nodules$stem_diameter_cm==1 & nodules$nodule_count==27]<-4
#C 70-4-1 should be 70-1-4
nodules$plot[nodules$stem_diameter_cm==.7 & nodules$nodule_count==44]<-1
nodules$indiv[nodules$stem_diameter_cm==.7 & nodules$nodule_count==44]<-4
#C 50-6-1 should be 50-1-6
nodules$plot[nodules$stem_diameter_cm==1.1 & nodules$nodule_count==48]<-1
nodules$indiv[nodules$stem_diameter_cm==1.1 & nodules$nodule_count==48]<-6
#C 70-6-1 should be 70-1-6
nodules$plot[nodules$stem_diameter_cm==1.1 & nodules$nodule_count==28]<-1
nodules$indiv[nodules$stem_diameter_cm==1.1 & nodules$nodule_count==28]<-6
#missing C 70-1-3 and duplicate 70-3-1
#70-1-3 larger biomass so assign it the larger nodule count
nodules$plot[nodules$stem_diameter_cm==1.0 & nodules$nodule_count==44]<-1
nodules$indiv[nodules$stem_diameter_cm==1.0 & nodules$nodule_count==44]<-3
#missing C 50-1-3 and duplicate 50-3-1
#50-1-3 smaller biomass so assign it the smaller stem diameter
nodules$plot[nodules$stem_diameter_cm==.8 & nodules$nodule_count==78]<-1
nodules$indiv[nodules$stem_diameter_cm==.8 & nodules$nodule_count==78]<-3
#missing C 5-3-4 and duplicate 5-3-6
#5-3-4 larger biomass so assign it the larger stem diameter
nodules$indiv[nodules$stem_diameter_cm==1.2 & nodules$nodule_count==78 & nodules$indiv==6]<-4
#missing PH 5-1-6 and duplicate 50-1-6
#5-1-6 larger biomass so assign it the larger nodule count
nodules$variety[nodules$stem_diameter_cm==.5 & nodules$nodule_count==37 & nodules$variety==50]<-5

harvest<-left_join(harvest,nodules,by=c("site", "variety", "plot", "indiv"))

#add in pheno stage at harvest so we can group by it
pheno<-read.csv("clean_data/clean_pheno_2020.csv")
harvest<-left_join(harvest,pheno,by=c("site", "variety", "plot", "indiv"))

#remove extra columns
harvest <- select(harvest,-c(Notes,date))

#output cleaned file, should have 216 rows total
write.csv(harvest,file="clean_data/clean_harvest_2020.csv",row.names = F)


####2021####
#aboveground biomass, bean yield, bean counts
data <- read_excel("raw_data/2021_AGbiomass_counts.xlsx")
data <- select(data,-date,-notes)
data<-rename(data,plant_biomass_g=dry_g)

#add in nodule data
nodules<-read_excel("raw_data/beanDIP 2021 nodule and root collar diameter.xlsx")
#duplicate PH 58-1-i3 and no 58-3-i1
#in aboveground biomass PH 58-1-i3 is largest of all plot 1 so the nodules row with nodule_count 78 and 0.9 RCD must be true value
#change row with nodule_count 66 and 0.4 RCD to plot 3, indiv 1
nodules$plot[nodules$site=="PH" & nodules$variety==58 & nodules$nodule_count==66 & nodules$root_collar_diameter_cm==0.4]<-3
nodules$indiv[nodules$site=="PH" & nodules$variety==58 & nodules$nodule_count==66 & nodules$root_collar_diameter_cm==0.4]<-1
harvest<-left_join(data,nodules,by=c("site", "variety", "plot", "indiv"))

#add in pheno stage at harvest so we can group by it
pheno<-read.csv("clean_data/clean_pheno_2021.csv")
#select only harvest date
pheno<-filter(pheno,sampling_round==4)
harvest<-left_join(harvest,pheno,by=c("site", "variety", "plot", "indiv"))
harvest<-select(harvest,-c(date,notes.x,notes.y,sampling_round))

#add in root biomass data
roots<-read_excel("raw_data/beanDIP 2021 root biomass weight data.xlsx")
roots<-roots %>% 
  rename(root_biomass_g=weight_g,notes="...6") %>%
  select(-notes)
#change label P to PH in site
roots$site[roots$site=="P"]<-"PH"
# duplicate PH 58-1-3 roots, and PH 58-3-1 root missing
# based on other metrics, 58-1-3 is the bigger number
roots$plot[roots$site=="PH" & roots$variety==58 & roots$root_biomass_g==0.9129]<-3
roots$indiv[roots$site=="PH" & roots$variety==58 & roots$root_biomass_g==0.9129]<-1
harvest<-left_join(harvest,roots,by=c("site", "variety", "plot", "indiv"))

#output cleaned file, should have 264 rows total
write.csv(harvest,file="clean_data/clean_harvest_2021.csv",row.names = F)

#### 2022 ####
#aboveground biomass
data <- read_excel("raw_data/beanDIP_2022_aboveground_soybean_biomass.xlsx")
data<-rename(data,plant_biomass_g=dry_mass_g,notes=`...7`)
data<-select(data,-notes)

#add in nodule data
nodules<-read_excel("raw_data/beanDIP 2022 nodule and root collar data.xlsx")
# duplicate values for PH plot 3 indiv 4 - missing value for W plot 3 indiv 4. based on other plant size values, PH is collar 9.6 - W is 11.7
nodules$site[nodules$plot==3 & nodules$indiv==4 & nodules$root_collar_diameter_mm == 11.7]<-"W"
nodules<-select(nodules,-c(notes,harvest_date))
harvest<-left_join(data,nodules,by=c("site", "variety", "plot", "indiv"))

#add in pheno stage at harvest so we can group by it
pheno<-read.csv("clean_data/clean_pheno_2022.csv")
#select only harvest date
pheno<-filter(pheno,sampling_round==4)
harvest<-left_join(harvest,pheno,by=c("site", "variety", "plot", "indiv"))
harvest<-select(harvest,-c(date,notes.x,notes.y,sampling_round))

#add in root biomass data
roots<-read_excel("raw_data/beanDIP 2022 root biomass weight data.xlsx")
roots<-rename(roots,root_biomass_g=weight_g)
#duplicate values for PH plot 1 indiv 1 - missing value for PH plot 3 indiv 1. based on other plant size values, 2.28 is plot 1 - 0.93 is plot 3
roots$plot[roots$site=="PH" & roots$indiv==1 & roots$root_biomass_g == 0.9363]<-3
harvest<-left_join(harvest,roots,by=c("site", "variety", "plot", "indiv"))

#output cleaned file, should have 72 rows total
write.csv(harvest,file="clean_data/clean_harvest_2022.csv",row.names = F)

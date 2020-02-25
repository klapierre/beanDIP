####setup####
library(readxl)
library(tidyverse)
setwd("~/Dropbox/bean_dip_2018-2024/data")

####data wrangling####
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




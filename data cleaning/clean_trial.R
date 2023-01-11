################################################################################
##  clean_trial.R: Cleaning trial data.
##
##  Author: Kelsey McGurrin
################################################################################


####setup####
library(readxl)
library(tidyverse)
setwd("~/Dropbox/bean_dip_2018-2024/field trials/data")

####2019####
#import each site data from raw
w <- read_excel("raw_data/StatewideSoyTrialResultsByPlot_2019.xlsx",sheet="Wye DC_2019")
c <- read_excel("raw_data/StatewideSoyTrialResultsByPlot_2019.xlsx",sheet="CV")
ph <- read_excel("raw_data/KV_PHFS _G3.xlsx",sheet="PH_FS_G3")
k <- read_excel("raw_data/KV_PHFS _G3.xlsx",sheet="KV_G3")

#select subset of columns and rows
our_vars<-c(31,32,55,27,67,83)

w <-w %>%
  select(plot=Rep,variety="Entry#",trial_yield=Yield) %>%
  add_column(site="W") %>%
  filter(variety %in% our_vars)

ph <-ph %>%
  select(plot=Rep,variety="Entry#",trial_yield="Yield (bu/a)") %>%
  add_column(site="PH")%>%
  filter(variety %in% our_vars)

k$variety<-as.numeric(k$`Entry#`)
k <-k %>%
  select(plot=Rep,variety,trial_yield=Yield) %>%
  add_column(site="K") %>%
  filter(variety %in% our_vars)

c <-c %>%
  select(plot=Rep,variety="Entry#",trial_yield=Yield) %>%
  add_column(site="C")%>%
  filter(variety %in% our_vars)

#join all sites together- should have 66 rows?
all<-bind_rows(k,c,w,ph)
all<-arrange(all,site,variety,plot)

#export
write.csv(all,file="clean_data/clean_trial_harvest_2019.csv",row.names = F)

####2020####
data<-read_excel("raw_data/soybean VT 2020 for Kelsey.xlsx")

data<-data %>%
  filter(TEST=="FS") %>%
  select(variety=Entry,plot=Rep,trial_yield=bu_per_ac,site=LOCATION) 
data$site<-recode(data$site,"KV"="K","CV"="C","WYE"="W")

#export
write.csv(data,file="clean_data/clean_trial_harvest_2020.csv",row.names = F)

####2021####
data<-read_excel("raw_data/2021 SVT yields for Kelsey.xlsx")

data<-data %>%
  filter(TEST=="FS") %>%
  select(variety=Entry,plot=Rep,trial_yield=bu_per_ac,site=LOCATION) 
data$site<-recode(data$site,"KV"="K","CV"="C","WYE"="W")

#export
write.csv(data,file="clean_data/clean_trial_harvest_2021.csv",row.names = F)

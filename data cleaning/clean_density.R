################################################################################
##  clean_density.R: Cleaning density data.
##
##  Author: Kelsey McGurrin
################################################################################

####setup####
library(readxl)
library(tidyverse)
setwd("~/Dropbox/bean_dip_2018-2024/field trials/data")

####2020####
#separate and clean density data
damage<- read_excel("raw_data/2020_damage_soil_data.xlsx",sheet = "damage_row")
damage<-rename(damage,samp_round="sampling round")

#include sampling round because have both initial and harvest data
dens<-select(damage,samp_round,site,variety,plot,row,density)

#clarksville variety 2, plot 1, row b was marked incorrectly at first sampling- use second round for initial density count
#look at what the initial count should be (15)
mixup<-filter(dens,site=="C" & variety==2 & plot==1 & row=="b")
#change values
dens$density[dens$site=="C" & dens$variety==2 & dens$plot==1 & dens$row=="b" & dens$samp_round==1]<-15
dens$density[dens$site=="C" & dens$variety==2 & dens$plot==1 & dens$row=="b" & dens$samp_round==2]<-NA

#make new df for corrected data with 1 column initial, 1 final
densb<-filter(dens,samp_round %in% c(1,2,5))
densb$timing<-recode(densb$samp_round,
                    "1"="initial","2"="initial","5"="final")

densb<-densb %>%
  group_by(site,variety,plot,row,timing) %>%
  summarise(density = mean(density,na.rm = T))

densb<-pivot_wider(densb,names_from = timing,values_from = density)

#check for changes in density
densb<-mutate(densb,diff=final-initial)

#output clean csv
write.csv(densb,file="clean_data/clean_density_2020.csv",row.names = F)

####2021####
#separate and clean density data
damage<- read_excel("raw_data/svt_2021_damage_soil_data.xlsx",sheet = "damage_row")
damage<-rename(damage,samp_round="sampling round")

#include sampling round because have both initial and harvest data
dens<-select(damage,samp_round,site,variety,plot,row,density)

#make new df for data with 1 column initial, 1 final
densb<-filter(dens,samp_round %in% c(1,4))
densb$timing<-recode(densb$samp_round,
                     "1"="initial","4"="final")

densb<-densb %>%
  group_by(site,variety,plot,row,timing) %>%
  summarise(density = mean(density,na.rm = T))

densb<-pivot_wider(densb,names_from = timing,values_from = density)

#check for changes in density
densb<-mutate(densb,diff=final-initial)

#output clean csv
write.csv(densb,file="clean_data/clean_density_2021.csv",row.names = F)

####2022####
#separate and clean density data
damage<- read_excel("raw_data/svt_2022_damage_soil_data.xlsx",sheet = "damage_row")
damage<-rename(damage,samp_round="sampling round")

#include sampling round because have both initial and harvest data
dens<-select(damage,samp_round,site,variety,plot,row,density)

#make new df for data with 1 column initial, 1 final
densb<-filter(dens,samp_round %in% c(1,4))
densb$timing<-recode(densb$samp_round,
                     "1"="initial","4"="final")

densb<-densb %>%
  group_by(site,variety,plot,row,timing) %>%
  summarise(density = mean(density,na.rm = T))

densb<-pivot_wider(densb,names_from = timing,values_from = density)

#check for changes in density
densb<-mutate(densb,diff=final-initial)

#output clean csv (should be 24 rows)
write.csv(densb,file="clean_data/clean_density_2022.csv",row.names = F)

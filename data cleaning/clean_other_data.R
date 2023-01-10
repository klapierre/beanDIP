################################################################################
##  clean_other_data.R: Cleaning phenology, damage, and soil moisture data.
##
##  Author: Kelsey McGurrin
################################################################################

####setup####
library(readxl)
library(tidyverse)
setwd("~/Dropbox/bean_dip_2018-2024/field trials/data") #kelsey wd
setwd('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\bean_dip_2018-2024\\field trials\\data') #kim wd


####2019####
#separate and clean damage data
damage<- read_excel("raw_data/2019_damage_soil_data.xlsx",sheet = "damage")
#could also split up insect community data someday?
write.csv(damage,file="clean_data/clean_damage_2019.csv",row.names = F)


#separate and clean phenology data
pheno<- read_excel("raw_data/2019_damage_soil_data.xlsx",sheet = "phenology")
pheno<-separate(pheno, pheno_stage, sep = "R",into = c("junk", "pheno_stage"))
pheno$junk<-NULL
write.csv(pheno,file="clean_data/clean_pheno_2019.csv",row.names = F)


#separate and clean soil moisture data
TDR<-read_excel("raw_data/2019_damage_soil_data.xlsx",sheet = "TDR")
TDR<-mutate(TDR,TDR_avg=(TDR1+TDR2+TDR3)/3)
write.csv(TDR,file="clean_data/clean_soilmoisture_2019.csv",row.names = F)

####2020####
#separate and clean damage data
damage<- read_excel("raw_data/2020_damage_soil_data.xlsx",sheet = "damage_row")
#could also split up insect community data someday?
write.csv(damage,file="clean_data/clean_damage_2020.csv",row.names = F)

#also have damage per plant this year

#separate and clean phenology data
pheno<- read_excel("raw_data/2020_damage_soil_data.xlsx",sheet = "phenology")
pheno<-separate(pheno, pheno_stage, sep = "R",into = c("junk", "pheno_stage"))
pheno$junk<-NULL
write.csv(pheno,file="clean_data/clean_pheno_2020.csv",row.names = F)

####2021####
#separate and clean damage data
damage<- read_excel("raw_data/svt_2021_damage_soil_data.xlsx",sheet = "damage_row")
#could also split up insect community data someday?
write.csv(damage,file="clean_data/clean_damage_2021.csv",row.names = F)

#also have damage per plant this year

#separate and clean phenology data
pheno<- read_excel("raw_data/svt_2021_damage_soil_data.xlsx",sheet = "phenology")
pheno<-separate(pheno, pheno_stage, sep = "R",into = c("junk", "pheno_stage"))
pheno$junk<-NULL
write.csv(pheno,file="clean_data/clean_pheno_2021.csv",row.names = F)


#separate and clean soil moisture data
TDR<-read_excel("raw_data/svt_2021_damage_soil_data.xlsx",sheet = "TDR")
write.csv(TDR,file="clean_data/clean_soilmoisture_2021.csv",row.names = F)


####2022####
#separate and clean damage data
damage<- read_excel("raw_data/svt_2022_damage_soil_data.xlsx",sheet = "damage_row")
#could also split up insect community data someday?
write.csv(damage,file="clean_data/clean_damage_2022.csv",row.names = F)

#also have damage per plant this year

#separate and clean phenology data
pheno<- read_excel("raw_data/svt_2022_damage_soil_data.xlsx",sheet = "phenology")
pheno<-separate(pheno, pheno_stage, sep = "R",into = c("junk", "pheno_stage"))
pheno$junk<-NULL
write.csv(pheno,file="clean_data/clean_pheno_2022.csv",row.names = F)


#separate and clean soil moisture data
TDR<-read_excel("raw_data/svt_2022_damage_soil_data.xlsx",sheet = "TDR")
write.csv(TDR,file="clean_data/clean_soilmoisture_2022.csv",row.names = F)

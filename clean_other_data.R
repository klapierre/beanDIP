####setup####
library(readxl)
library(tidyverse)
setwd("~/Dropbox/bean_dip_2018-2024/data")

####data wrangling####
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

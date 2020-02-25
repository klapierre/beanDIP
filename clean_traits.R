####setup####
library(readxl)
library(tidyverse)
setwd("~/Dropbox/bean_dip_2018-2024/data")

####data wrangling####
#metrics from jamie
traits<-read_excel("raw_data/2019_LeafTraits.xlsx")
site_key <- c("clarksville" = "C","keedysville" = "K","poplar hill"= "PH","wye"="W")
traits$site<-recode(traits$site, !!!site_key)
traits$site<-factor(traits$site, levels = c("K", "C", "W", "PH"))
traits$variety<-as.factor(traits$variety)
traits<-rename(traits,indiv=individual)

#add in subset of clean photosynq data from same days
photosynq<-read.csv("clean_data/clean_photosynq_2019.csv")
photosynq<-rename(photosynq,site=Site,variety=Variety,plot=Plot,indiv=Individual)
site_key <- c("Clarksville" = "C","Keedysville" = "K","Poplar Hill"= "PH","Wye"="W")
photosynq$site<-recode(photosynq$site, !!!site_key)
photosynq$variety<-as.factor(photosynq$variety)
photosynq$time<-as.character(photosynq$time)
photosynq<-separate(photosynq, "time", sep = " ", into = c("date", "newtime", "AMPM"))
photosynq<-unite(photosynq,newtime,AMPM,col="time",sep = " ")
photosynq<-filter(photosynq, date %in% c("07/29/2019","07/30/2019","07/31/2019","08/01/2019","08/02/2019"))

#combine and output
all<-left_join(photosynq,traits,by=c("site", "variety","plot","indiv"))
all<-select(all,site,variety,plot,indiv,area_cm:dry_mass,date=date.x,time,Ambient.Humidity:vH.)
write.csv(all,file="clean_data/clean_traits_2019.csv",row.names = F)

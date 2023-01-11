################################################################################
##  clean_traits.R: Cleaning trait data.
##
##  Author: Kelsey McGurrin
################################################################################

####setup####
library(readxl)
library(tidyverse)
setwd("~/Dropbox/bean_dip_2018-2024/field trials/data")

####2019####
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

####2020####
#pressure chamber
lwp<-read_excel("raw_data/2020_LeafTraits_pressureChamber.xlsx")
lwp$variety<-as.factor(lwp$variety)
site_key <- c("clarksville" = "C","keedysville" = "K","poplar hill"= "PH","wye"="W")
lwp$site<-recode(lwp$site, !!!site_key)
lwp$site<-factor(lwp$site, levels = c("K", "C", "W", "PH"))
lwp<-rename(lwp,indiv=individual,pressure_chamber=lwp_bar)
lwp<-select(lwp,-c(date,notes))

#leaf area
area<-read_excel("raw_data/LeafArea_BeanDip2020_FieldTrials.xlsx")
area$variety<-as.factor(area$variety)
site_key <- c("Clarksville" = "C","Keedysville" = "K","Poplar hill"= "PH","Wye"="W")
area$site<-recode(area$site, !!!site_key)
area<-rename(area,indiv=leaf)
area <- select(area,c(site,variety,plot,indiv,Area))
#there are duplicate readings for most leaves, so take average and remove controls
area<-area %>%
  group_by(site,variety,plot,indiv) %>%
  summarise(area_cm = mean(Area)) %>%
  filter(variety %in% c(2, 5, 50, 70))

#other manual traits
traits<-read_excel("raw_data/LeafTraits_BeanDip2020_fieldTrials.xlsx")
traits$variety<-as.factor(traits$variety)
site_key <- c("Clarksville" = "C","Keedysville" = "K","Poplar hill"= "PH","wye"="W")
traits$site<-recode(traits$site, !!!site_key)
traits$site<-factor(traits$site, levels = c("K", "C", "W", "PH"))
traits<-rename(traits,indiv=leaf,wet_mass=wet_mass_g,toughness_1=toughness1,toughness_2=toughness2,dry_mass=dry_mass_g)
traits<-select(traits,-picture)

#add in clean photosynq data
photosynq<-read.csv("clean_data/clean_photosynq_2020.csv")
photosynq$variety<-as.factor(photosynq$variety)
photosynq<-rename(photosynq,site=Site,variety=Variety,plot=Plot,indiv=Individual)
site_key <- c("Clarksville" = "C","Keedysville" = "K","Poplar Hill"= "PH","Wye"="W")
photosynq$site<-recode(photosynq$site, !!!site_key)
photosynq<-select(photosynq,Ambient.Humidity:Time.of.Day)

#combine and output
all<-left_join(lwp,area,by=c("site", "variety","plot","indiv"))
all<-left_join(all,traits,by=c("site", "variety","plot","indiv"))
all<-left_join(all,photosynq,by=c("site", "variety","plot","indiv"))

write.csv(all,file="clean_data/clean_traits_2020.csv",row.names = F)

####2021####
#pressure chamber
lwp<-read_excel("raw_data/2021_LeafTraits_pressureChamber.xlsx")
# lwp$variety<-as.factor(lwp$variety)
site_key <- c("clarksville" = "C","keedysville" = "K","poplar hill"= "PH","wye"="W")
lwp$site<-recode(lwp$site, !!!site_key)
lwp$site<-factor(lwp$site, levels = c("K", "C", "W", "PH"))
lwp<-rename(lwp,indiv=individual,pressure_chamber=lwp_bar)
lwp<-select(lwp,-c(date,notes))

#other manual traits
traits<-read_excel("raw_data/LeafTraits_BeanDip2021_fieldTrials.xlsx")
site_key <- c("keedysville" = "K","clarksville" = "C","wye"="W","poplar hill"= "PH")
traits$site<-recode(traits$site, !!!site_key)
traits<-rename(traits,indiv=leaf,wet_mass=wet_mass_g,toughness_1=toughness1,toughness_2=toughness2,dry_mass=dry_mass_g)


#add in clean photosynq data
photosynq<-read_csv("clean_data/clean_photosynq_2021.csv")
photosynq<-rename(photosynq,site=Site,variety=Variety,plot=Plot,indiv=Individual)
site_key <- c("Keedysville" = "K","Clarksville" = "C","Wye"="W","Poplar Hill"= "PH")
photosynq$site<-recode(photosynq$site, !!!site_key)
photosynq<-select(photosynq,Ambient.Humidity:Time.of.Day)

#combine and output
all<-left_join(lwp,photosynq,by=c("site", "variety","plot","indiv"))
all<-left_join(traits,all,by=c("site", "variety","plot","indiv"))

#still has 2 extra rows?
write.csv(all,file="clean_data/clean_traits_2021.csv",row.names = F)

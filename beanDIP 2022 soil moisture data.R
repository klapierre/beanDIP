#December 8, 2022
#Sarah Alley 
#beanDIP 2022 soil moisture data

#Notes on how to approach this
#First, review all HOBO files from the summer
#Look at the data to make sure it looks legit- abrupt drops or funky patterns can indicate failure
#If the logger was doing ok for a while and has a point where you can tell an error occured, you can save the data before that point
#Click on the graph where you'd like to crop, then right click, crop series and enter info for the point where you want it to crop (or something like that)
#Export the new cropped data as a csv 

#Notes
#K 73_3, PH 73_1. W 73_1, W 73_3 were total failures, no data coming from those loggers
#C 73_1, PH 73_3 were partial failures, logged successfully as of August 1, then failed 

#Still need to reformat the data to reuse code from the past 
#For soil moisture data, name columns the following: obs, date, time, am_pm, water_content 
#Add two new columns to the right of the date column
#Data, text to columns, delimited, check space (will separate into date, time, am/pm)
#Control + 1, use this to put the date column in mm/dd/yyyy format 

library(tidyverse)

install.packages("chron")
library(chron)

#Laptop


setwd("C:/Users/Sarah Alley/Dropbox (Smithsonian)/bean_dip_2018-2024/field trials/data/raw_data/HOBO soil moisture and temperature data/Soil Moisture 2022/Soil Moisture 2022 ready for R")

#Reading in data

C73_1<-read.csv("Cville_73_1_for_R.csv")%>%
  mutate(site="Clarksville", plot=1) %>% 
  rename(obs=ï..obs)

C73_2<-read.csv("Cville_73_2_for_R.csv")%>%
  mutate(site="Clarksville", plot=2) %>% 
  rename(obs=ï..obs)

C73_3<-read.csv("Cville_73_3_for_R.csv")%>%
  mutate(site="Clarksville", plot=3) %>% 
  rename(obs=ï..obs)

K73_1<-read.csv("Kville_73_1_for_R.csv")%>%
  mutate(site="Keedysville", plot=1) %>% 
  rename(obs=ï..obs)

K73_2<-read.csv("Kville_73_2_for_R.csv")%>%
  mutate(site="Keedysville", plot=2) %>% 
  rename(obs=ï..obs)

PH73_2<-read.csv("PopHill_73_2_for_R.csv")%>%
  mutate(site="Poplar Hill", plot=2) %>% 
  rename(obs=ï..obs)

PH73_3<-read.csv("PopHill_73_3_for_R.csv")%>%
  mutate(site="Poplar Hill", plot=3) %>% 
  rename(obs=ï..obs)

W73_2<-read.csv("Wye_73_2_for_R.csv")%>%
  mutate(site="Wye", plot=2) %>% 
  rename(obs=ï..obs)

#Combining into one big dataset

#Having issues with dates, when converting to 24 hour time, R is giving 0022 as the year instead of 2022 

#soilmoisture2022 <- C73_1%>%rbind(C73_2)%>%rbind(C73_3)%>%rbind(K73_1)%>%rbind(K73_2)%>%
  #rbind(PH73_2)%>%rbind(PH73_3)%>%rbind(W73_2)%>%
  #unite(time_am_pm, c(date, time, am_pm), sep=" ") %>% 
  #mutate(time_24hr=as.POSIXct(strptime(time_am_pm,"%m/%d/%Y %I:%M:%S %p")))


#soilmoisture2022 <- C73_1%>%rbind(C73_2)%>%rbind(C73_3)%>%rbind(K73_1)%>%rbind(K73_2)%>%
  #rbind(PH73_2)%>%rbind(PH73_3)%>%rbind(W73_2)%>%
  #unite(time_am_pm, c(date, time, am_pm), sep=" ") %>% 
  #mutate(time_24hr=as.Date(time_am_pm,"%m/%d/%Y %I:%M:%S %p"), origin = "2022-06-24")


#This one works!!! Change the format of the time before pasting on to date! 

soilmoisture2022_better <- C73_1%>%rbind(C73_2)%>%rbind(C73_3)%>%rbind(K73_1)%>%rbind(K73_2)%>%
  rbind(PH73_2)%>%rbind(PH73_3)%>%rbind(W73_2)%>%
  unite(time_am_pm, c(time, am_pm), sep=" ") %>% 
  mutate(time_24=format(strptime(time_am_pm,"%I:%M:%S %p"),"%H:%M:%S"))%>%
  unite(time_24_date, c(date, time_24), sep = " ")
  
  
  
  mutate(time_24hr=as.POSIXct(time_am_pm,"%m/%d/%Y %I:%M:%S %p"))


#Graph
#Will not work if line 100 is run, still need to fix x axis so dates are more visible 
ggplot(data=subset(soilmoisture2022_better,water_content > 0.1),aes(x=time_24_date, y=water_content,color=as.factor(plot))) + 
  geom_point()+
  facet_wrap(~site)+
  #scale_x_datetime(date_labels= "%m/%d/%Y %HH:%MM:%SS")+
  ylab('Water Content') + xlab('Date')+
  ggtitle("beanDIP 2022 Soil Moisture Data")+
  guides(col=guide_legend("Plot"))


ggplot(data=subset(soilmoisture2022_better,water_content > 0.1),aes(x=time_24_date, y=water_content,color=as.factor(plot))) + 
  geom_line()+
  facet_wrap(~site)+
  scale_x_datetime(date_labels= "%Y-%m-%d %HH:%MM:%SS")+
  ylab('Water Content') + xlab('Date')+
  ggtitle("beanDIP 2022 Soil Moisture Data")+
  guides(col=guide_legend("Plot"))

ggplot(data=subset(soilmoisture2022_better,water_content > 0.1),aes(x=time_24_date, y=water_content,color=as.factor(plot))) + 
  geom_line()+
  facet_wrap(~site)+
  ylab('Water Content') + xlab('Date')+
  ggtitle("beanDIP 2022 Soil Moisture Data")+
  guides(col=guide_legend("Plot"))

#This graph works (kinda)
#Date doesn't appear on x-axis, and looks like too many measurements at a given time point?
#Best I've been able to do, added the group = 1 statement 
ggplot(data=soilmoisture2022_better,aes(x=time_24_date, y=water_content, group = 1, color=as.factor(plot))) + 
  geom_line()+
  facet_wrap(~site)+
  ylab('Water Content') + xlab('Date')+
  ggtitle("beanDIP 2022 Soil Moisture Data")+
  guides(col=guide_legend("Plot"))
  

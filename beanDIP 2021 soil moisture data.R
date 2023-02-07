#beanDIP 2021 temperature data
#Sarah Alley
#January 5, 2021

library(tidyverse)

#For future reference, to format data in Excel for use in R:
#For soil moisture data, name columns the following: obs, date, time, am_pm, water_content 
#Add two new columns to the right of the date column
#Data, text to columns, delimited, check space (will separate into date, time, am/pm)
#Control + 1, use this to put the date column in mm/dd/yyyy format 


#A couple of the sensors had issues 
#PH 19-3: don't think this is worth including, logged for less than a day before it stopped 
#PH 19-1: break in the data, stopped logging (I know Kelsey checked them on 7/27, day it stopped)
#C 19-1: only logged up until 7/28 (a day Kelsey checked them), I wonder if this is the one that got run over and wouldn't connect to my computer? 
#K 19-3: unsure why this happened? 

#For C, K, and W, the loggers seem to be following the same pattern of precipitation and drying out events
#PH: same pattern, but values off 

#Laptop

setwd("C:/Users/Sarah Alley/Dropbox (Smithsonian)/bean_dip_2018-2024/field trials/data/raw_data/HOBO soil moisture and temperature data/Soil Moisture 2021/Ready for R")


C19_1<-read.csv("Clarksville_19-1_forR.csv")%>%
  mutate(site="Clarksville", plot=1) %>% 
  rename(obs=ï..obs)

C19_2<-read.csv("Clarksville_19-2_forR.csv")%>%
  mutate(site="Clarksville", plot=2) %>% 
  rename(obs=ï..obs)

C19_3<-read.csv("Clarksville_19-3_forR.csv")%>%
  mutate(site="Clarksville", plot=3) %>% 
  rename(obs=ï..obs)

K19_1<-read.csv("Keedysville_19-1_forR.csv")%>%
  mutate(site="Keedysville", plot=1) %>% 
  rename(obs=ï..obs)

K19_2<-read.csv("Keedsyville_19-2_forR.csv")%>%
  mutate(site="Keedysville", plot=2) %>% 
  rename(obs=ï..obs)

K19_3<-read.csv("Keedysville_19-3_forR.csv")%>%
  mutate(site="Keedysville", plot=3) %>% 
  rename(obs=ï..obs)

P19_1<-read.csv("Poplar_Hill_19-1_forR.csv")%>%
  mutate(site="Poplar Hill", plot=1) %>% 
  rename(obs=ï..obs)

P19_2<-read.csv("Poplar_Hill_19-2_forR.csv")%>%
  mutate(site="Poplar Hill", plot=2) %>% 
  rename(obs=ï..obs)

W19_1<-read.csv("Wye_19-1_forR.csv")%>%
  mutate(site="Wye", plot=1) %>% 
  rename(obs=ï..obs)

W19_2<-read.csv("Wye_19-2_forR.csv")%>%
  mutate(site="Wye", plot=2) %>% 
  rename(obs=ï..obs)

W19_3<-read.csv("Wye_19-3_forR.csv")%>%
  mutate(site="Wye", plot=3) %>% 
  rename(obs=ï..obs)

#Hm... so we formatted our data to have date, time, am_pm in separate columns (guess I didn't need to do that)
#Here it's using a unite statement to put them back together

soilmoisture <- C19_1%>%rbind(C19_2)%>%rbind(C19_3)%>%rbind(K19_1)%>%rbind(K19_2)%>%
  rbind(K19_3)%>%rbind(P19_1)%>%rbind(P19_2)%>%rbind(W19_1)%>%rbind(W19_2)%>%rbind(W19_3)%>%
  unite(time_am_pm, c(date, time, am_pm), sep=" ") %>% 
  mutate(time_24hr=as.POSIXct(strptime(time_am_pm,"%m/%d/%Y %I:%M:%S %p")))

#Here's a graph! 

ggplot(data=subset(soilmoisture,water_content > 0.1),aes(x=time_24hr, y=water_content,color=as.factor(plot))) + 
  geom_line()+
  scale_x_datetime(date_labels= "%Y-%m-%d %HH:%MM:%SS")+
  facet_wrap(~site)+
  ylab('Water Content') + xlab('Date')+
  ggtitle("beanDIP 2021 Soil Moisture Data")+
  guides(col=guide_legend("Plot"))

#This one the date is easier to look at 

ggplot(data=subset(soilmoisture,water_content > 0.1),aes(x=time_24hr, y=water_content,color=as.factor(plot))) + 
  geom_line()+
  facet_wrap(~site)+
  ylab('Water Content') + xlab('Date')+
  ggtitle("beanDIP 2021 Soil Moisture Data")+
  guides(col=guide_legend("Plot"))

#K 19-3 is pretty messed up, graph below dropping that data

soilmoisturewithoutK19_3 <- C19_1%>%rbind(C19_2)%>%rbind(C19_3)%>%rbind(K19_1)%>%rbind(K19_2)%>%
  rbind(P19_1)%>%rbind(P19_2)%>%rbind(W19_1)%>%rbind(W19_2)%>%rbind(W19_3)%>%
  unite(time_am_pm, c(date, time, am_pm), sep=" ") %>% 
  mutate(time_24hr=as.POSIXct(strptime(time_am_pm,"%m/%d/%Y %I:%M:%S %p")))

ggplot(data=subset(soilmoisturewithoutK19_3,water_content > 0.1),aes(x=time_24hr, y=water_content,color=as.factor(plot))) + 
  geom_line()+
  facet_wrap(~site)+
  ylab('Water Content') + xlab('Date')+
  ggtitle("beanDIP 2021 Soil Moisture Data")+
  guides(col=guide_legend("Plot"))

  
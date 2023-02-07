#beanDIP 2021 temperature data
#Sarah Alley
#January 4, 2021 


#For future reference, to format data in Excel for use in R:
#Make sure to delete any notes like "logged" in far left column
#For temperature data, name columns the following: obs, date, time, am_pm, temp, rh
#Create two blank columns beside the Date Time column
#Select that entire column 
#In the data tab, go to text to columns, delimited, next
#Check space as delimiter
#This should separate out the one column into separate columns for date, time, and am/pm
#For some reason 12:00 will appear next to all the dates, to fix, hold control +1, and select 

library(tidyverse)
library(ggplot2)


#Set working directory

#Laptop
setwd("C:/Users/Sarah Alley/Dropbox (Smithsonian)/bean_dip_2018-2024/field trials/data/raw_data/HOBO soil moisture and temperature data/Temperature 2021/Ready for R")


#Reading in data

Cs<-read.csv("Clarksville_Shielded_2021_forR.csv")%>%
  mutate(site="Clarksville",shielded="yes") %>% 
  rename(obs=ï..obs)

Cu<-read.csv("Clarksville_Unshielded_2021_forR.csv")%>%
  mutate(site="Clarksville",shielded="no") %>% 
  rename(obs=ï..obs)

Ks<-read.csv("Keedysville_Shielded_2021_forR.csv")%>%
  mutate(site="Keedysville",shielded="yes") %>% 
  rename(obs=ï..obs)

Ku<-read.csv("Keedysville_Unshielded_2021_forR.csv")%>%
  mutate(site="Keedysville",shielded="no") %>% 
  rename(obs=ï..obs)

Ws<-read.csv("Wye_Shielded_2021_forR.csv")%>%
  mutate(site="Wye",shielded="yes") %>% 
  rename(obs=ï..obs)

Wu<-read.csv("Wye_Unshielded_2021_forR.csv")%>%
  mutate(site="Wye",shielded="no") %>% 
  rename(obs=ï..obs)

PHs<-read.csv("Poplar_Hill_Shielded_2021_forR.csv")%>%
  mutate(site="Poplar Hill",shielded="yes") %>% 
  rename(obs=ï..obs)

PHu<-read.csv("Poplar_Hill_Unshielded_2021_forR.csv")%>%
  mutate(site="Poplar Hill",shielded="no") %>% 
  rename(obs=ï..obs)

#Linking data together, switching to 24 hr time instead of am/pm

temperature <- Cs%>%rbind(Cu)%>%rbind(Ks)%>%rbind(Ku)%>%rbind(Ws)%>%
  rbind(Wu)%>%rbind(PHs)%>%rbind(PHu)%>%
  unite(time_am_pm, c(date, time, am_pm), sep=" ") %>% 
  mutate(time_24hr=as.POSIXct(strptime(time_am_pm,"%m/%d/%Y %I:%M:%S %p")))

#Converting temperatures from C to F
temperatureF <- mutate(temperature, temp_F = temp*(9/5)+32)

temperatureFfinal <- select

#Exact same code worked for 2020 data, not sure why it's giving an error message now 
ggplot(data=temperatureFfinal,aes(x=time_24hr, y=temp_F,color=as.factor(shielded))) + 
  geom_line()+
  ggtitle("beanDIP 2021 Temperature Data")+
  guides(col=guide_legend("Shielded"))+
  ylab('Temperature (F)') + xlab('Date')+
  facet_wrap(~site)






#Last time we did this to pick the lower of the two temperature measurements at a given time at a given site for shielded vs. unshielded
#Since shielded is better for day, unshielded better for night 
#temperatureMin <-temperature%>%
  group_by(time_24hr, site)%>%
  summarize(temp_min=min(temp))%>%
  ungroup()
#I think last time I pre-set the loggers to all start at the same time, so the shielded and unshielded loggers take measurements at the exact same time
#Last time R was able to pick the lower value for a given time at a given site
#But this summer, I started the loggers in the field, which means the shielded and unshielded loggers were measuring at different times
#Since R has no way of associating them based on time, not sure how to select the lower value since measurements were not in pairs... 
  
#I need to add axis titles, make the legend look prettier
#I need to make this into degrees C instead of F to compare it to 2020 data? 
  

#Not sure why we did the scale_x_datetime thing last time... 

  ggplot(data=temperature,aes(x=time_24hr, y=temp,color=as.factor(shielded))) + 
    geom_line()+
    scale_x_datetime(date_labels= "%Y-%m-%d %HH:%MM:%SS")+
    facet_wrap(~site)+
    ylab('Temperature (F)') + xlab('Date')
  
#This graph leaves out the scale_x_datetime part and it's much easier to see the months... 
  
  
####TEMPERATURE####
  
  ggplot(data=temperature,aes(x=time_24hr, y=temp,color=as.factor(shielded))) + 
    geom_line()+
    facet_wrap(~site)+
    ylab('Temperature (F)') + xlab('Date')+
    ggtitle("beanDIP 2021 Temperatures")+
    guides(col=guide_legend("Shielded?"))

#Still need to work on converting this to degrees C... 
  
temperatureC %>%
  mutate(temp_C = ((temp-32)*(5/9)))
  


#subtract 32 and multiply by 5/9

####RELATIVE HUMIDITY####

#Going to use unshielded for humidity across the sites because outside humidity won't be affected by any moisture potentially trapped within the shield

ggplot(data=subset(temperature,shielded=="no"),aes(x=time_24hr, y=rh)) + 
  geom_line(color="black")+
  facet_wrap(~site)+
  ylab('Relative Humidity (%)') + xlab('Date')+
  ggtitle("beanDIP 2021 Relative Humidity (Unshielded Sensor)")

  

  




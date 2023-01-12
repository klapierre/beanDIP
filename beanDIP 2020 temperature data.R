#beanDIP 2020 temperature data

library(tidyverse)

setwd("C:\\Users\\Sarah Alley\\Documents\\Work Stuff\\HOBOware\\Temperature\\Temperature ready for R")

#Reading in data
#Ooh, ok yay, I got it to to work! Since whether it's shielded or unshielded is yes or no, a word, not a number, it needs to be in quotes. 
#My naming of the data files will be site initial with a s or a u following it for shielded or unshielded. 

Cs<-read.csv("Cshieldedcomplete for R.csv")%>%
  mutate(site="Clarksville",shielded="yes") %>% 
  rename(obs=ï..obs)

Cu<-read.csv("Cunshieldedcomplete for R.csv")%>%
  mutate(site="Clarksville",shielded="no") %>% 
  rename(obs=ï..obs)

Ks<-read.csv("Kshieldedcomplete for R.csv")%>%
  mutate(site="Keedysville",shielded="yes") %>% 
  rename(obs=ï..obs)

Ku<-read.csv("Kunshieldedcomplete for R.csv")%>%
  mutate(site="Keedysville",shielded="no") %>% 
  rename(obs=ï..obs)

Ws<-read.csv("Wshieldedcomplete for R.csv")%>%
  mutate(site="Wye",shielded="yes") %>% 
  rename(obs=ï..obs)

Wu<-read.csv("Wunshieldedcomplete for R.csv")%>%
  mutate(site="Wye",shielded="no") %>% 
  rename(obs=ï..obs)

PHs<-read.csv("PHshieldedcomplete forR.csv")%>%
  mutate(site="Poplar Hill",shielded="yes") %>% 
  rename(obs=ï..obs)

PHu<-read.csv("PHunshieldedcomplete for R.csv")%>%
  mutate(site="Poplar Hill",shielded="no") %>% 
  rename(obs=ï..obs)

#Now to link all the data together and to make the times 24-hr. 

temperature <- Cs%>%rbind(Cu)%>%rbind(Ks)%>%rbind(Ku)%>%rbind(Ws)%>%
  rbind(Wu)%>%rbind(PHs)%>%rbind(PHu)%>%
  unite(time_am_pm, c(date, time, am_pm), sep=" ") %>% 
  mutate(time_24hr=as.POSIXct(strptime(time_am_pm,"%m/%d/%Y %I:%M:%S %p")))

#Now to make graphs! We want shielded/unshielded on the same graph with one panel for each site. 
#Why does it say the 'subset' must be logical? I want the temperature data (named temp) from the data set temperature... I feel like this should work! 
#Would the factor here be "shielded"/"unshielded"? I think yes... 

#Subsetting removes rows, not columns, so this doesn't really work here. You tell R what you want to graph when you tell it what you want on the axes. 

ggplot(data=subset(temperature, temp),aes(x=time_24hr, y=temp,color=as.factor(shielded))) + 
  geom_line()+
  scale_x_datetime(date_labels= "%Y-%m-%d %HH:%MM:%SS")+
  facet_wrap(~site)

#Ok, well it works when I don't subset the data and tell it I want it to plot temperature as y=
#I'm also unsure of why it's not showing the unshielded data for Keedysville...
#Also, Keedysville and Clarksville have super high temperatures at the beignning because I set all the sensors to start running at 06/29 at 8 AM
#So they sat in a hot truck for a day before they were installed in the field, need to get rid of the first day's data 
#It's nice you don't really notice the two week gap of missing data too much! 

ggplot(data=temperature,aes(x=time_24hr, y=temp,color=as.factor(shielded))) + 
  geom_line()+
  scale_x_datetime(date_labels= "%Y-%m-%d %HH:%MM:%SS")+
  facet_wrap(~site)

#To show just one site

ggplot(data=subset(temperature,site=="Keedysville"),aes(x=time_24hr, y=temp,color=as.factor(shielded))) + 
  geom_line()+
  scale_x_datetime(date_labels= "%Y-%m-%d %HH:%MM:%SS")

#This is a way of getting rid of the dates before the sensors were actually installed in the field. 
 
ggplot(data=subset(temperature,time_24hr>as.Date("2020-06-30")),aes(x=time_24hr, y=temp,color=as.factor(shielded))) + 
  geom_line()+
  scale_x_datetime(date_labels= "%Y-%m-%d %HH:%MM:%SS")+
  facet_wrap(~site)

#Want to know minimum temperature at each site, will put it in a new column called temp_min
#always want to put ungroup after you group so things won't be messed up in the future 
temperatureMin <-temperature%>%
  group_by(time_24hr, site)%>%
  summarize(temp_min=min(temp))%>%
  ungroup()

#This graph shows the smaller of the two values at a given time for all the sites 
#The minimum value is the more accurate of the two temperatures, shielded better during the day, unshielded better at night
  
ggplot(data=temperatureMin, aes(x=time_24hr, y=temp_min))+
  geom_line()+
  scale_x_datetime(date_labels= "%Y-%m-%d %HH:%MM:%SS")+
  facet_wrap(~site)

#Down here I tried it for humidity... really difficult to look at  
ggplot(data=temperature,aes(x=time_24hr, y=rh,color=as.factor(shielded))) + 
  geom_line()+
  scale_x_datetime(date_labels= "%Y-%m-%d %HH:%MM:%SS")+
  facet_wrap(~site)

#Going to use unshielded for humidity across the sites because outside humidity won't be affected by any moisture potentially trapped within the shield
ggplot(data=subset(temperature,shielded=="no"),aes(x=time_24hr, y=rh)) + 
  geom_line(color="pink")+
  scale_x_datetime(date_labels= "%Y-%m-%d %HH:%MM:%SS")+
  facet_wrap(~site)


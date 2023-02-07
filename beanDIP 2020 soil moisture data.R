library(tidyverse)

#Laptop

setwd("C:/Users/Sarah Alley/Dropbox (Smithsonian)/bean_dip_2018-2024/field trials/data/raw_data/HOBO soil moisture and temperature data/Soil Moisture 2020/Soil moisture ready for R")

#Reading in data 
#Had to rename some of the observation columns because for some reason R thinks they were called obs=i, even though that's not what they're named in Excel file.  
K2_1<-read.csv("K2-1complete for R.csv") %>%
  mutate(site="Keedysville", plot=1)

#I think the sensor K 2-2 needs to be dropped, has a bunch of weird low values 

K2_2<-read.csv("K2-2complete for R.csv")%>%
  mutate(site="Keedysville", plot=2) %>% 
  rename(obs=ï..obs)

K2_3<-read.csv("K2-3complete for R.csv")%>%
  mutate(site="Keedysville", plot=3) %>%
  rename(obs=ï..obs)

C2_1<-read.csv("C2-1complete for R.csv") %>%
mutate(site="Clarksville", plot=1) %>%
  rename(obs=ï..obs)

C2_2<-read.csv("C2-2complete for R.csv") %>%
  mutate(site="Clarksville", plot=2) %>%
  rename(obs=ï..obs)

C2_3<-read.csv("C2-3complete for R.csv") %>%
  mutate(site="Clarksville", plot=3) %>%
  rename(obs=ï..obs)

W2_1<-read.csv("W2-1complete for R.csv") %>%
  mutate(site="Wye", plot=1) %>%
rename(obs=ï..obs)

W2_2<-read.csv("W2-2complete for R.csv") %>%
  mutate(site="Wye", plot=2) %>%
rename(obs=ï..obs)

W2_3<-read.csv("W2-3complete for R.csv") %>%
  mutate(site="Wye", plot=3) %>%
rename(obs=ï..obs)

#For Poplar Hill data (where all 3 sensors were either chewed or thrown out of the ground):
#2-1, the big chunk of missing data in the middle of the summer is due to the dates being thrown way off for some reason... 
#2-2, some weird really low data, but could still be cleaned up and useable in R?
#2-3, data plummets and ends early 8/22, but some still may be useable 


PH2_1<-read.csv("PH2-1useable for R.csv") %>%
  mutate(site="Poplar Hill", plot=1)%>% 
  rename(obs=ï..obs)

PH2_2<-read.csv("PH2-2useable for R.csv") %>%
  mutate(site="Poplar Hill", plot=2)%>% 
  rename(obs=ï..obs)

PH2_3<-read.csv("PH2-3useable for R.csv") %>%
  mutate(site="Poplar Hill", plot=3) %>%
  rename(obs=ï..obs)
  

#Combining soil moisture data from the 3 plots and making the am/pm times into 24 hour times 
#Originally I separated date/time/ampm into separate columns in excel, but that was actually unnecessary, in the code there's a step to put those columns back together (unite)
#R was being a wierdo about this POSIXct or POSIXlt format of the times, so that's what that step means 
#We got rid of the Wye plot 3 soil moisture data because it was so different than the other two sensors at the site 


soilmoisture <- K2_1%>%rbind(K2_2)%>%rbind(K2_3)%>%rbind(C2_1)%>%rbind(C2_2)%>%
  rbind(C2_3)%>%rbind(W2_1)%>%rbind(W2_2)%>%rbind(W2_3)%>%
  unite(time_am_pm, c(date, time, am_pm), sep=" ") %>% 
  mutate(time_24hr=as.POSIXct(strptime(time_am_pm,"%m/%d/%Y %I:%M:%S %p")))%>%
  mutate(drop=ifelse(site=="Wye"&plot==3,"drop", "keep"))

#Wanted to get rid of Keedysville 2-2 and didn't want to mess up original code, added it below
#Not sure why it's not happy since it's *exactly* the same as what Kim did for Wye, says object 'site' not found... 
#Should I remove Wye and Keedysville in the same mutate step? Not sure how I'd indicate which plots I want removed for what sites, though

soilmoisture <- K2_1%>%rbind(K2_2)%>%rbind(K2_3)%>%rbind(C2_1)%>%rbind(C2_2)%>%
  rbind(C2_3)%>%rbind(W2_1)%>%rbind(W2_2)%>%rbind(W2_3)%>%
  unite(time_am_pm, c(date, time, am_pm), sep=" ") %>% 
  mutate(time_24hr=as.POSIXct(strptime(time_am_pm,"%m/%d/%Y %I:%M:%S %p")))%>%
  mutate(drop=ifelse(site=="Wye"&plot==3,"drop", "keep"))
  mutate(drop=ifelse(site=="Keedysville"&plot==2, "drop","keep"))
  
#I'll come back to trying to get rid of the Keedysville data. I don't want to mess up original code, so I'm going to try to add the PH data to the rest below.
  #All I did was add more data... why does it have a problem with the unite statement now when it doesn't in the original code?
  #OK, it works, I just forgot to add a pipe on the end (eyeroll)

  soilmoisture1 <- K2_1%>%rbind(K2_2)%>%rbind(K2_3)%>%rbind(C2_1)%>%rbind(C2_2)%>%
    rbind(C2_3)%>%rbind(W2_1)%>%rbind(W2_2)%>%rbind(W2_3)%>%rbind(PH2_1)%>%rbind(PH2_2)%>%rbind(PH2_3)%>%
    unite(time_am_pm, c(date, time, am_pm), sep=" ") %>% 
    mutate(time_24hr=as.POSIXct(strptime(time_am_pm,"%m/%d/%Y %I:%M:%S %p")))%>%
    mutate(drop=ifelse(site=="Wye"&plot==3,"drop", "keep"))
  

#Making graphs! Facet wrap means is makes a separate panel for each site (this code with the dataset soilmoisture is just for W, C, and K)

ggplot(data=subset(soilmoisture1,drop!="drop"),aes(x=time_24hr, y=water_content,color=as.factor(plot))) + 
  geom_line()+
  scale_x_datetime(date_labels= "%Y-%m-%d %HH:%MM:%SS")+
  facet_wrap(~site)

#Easier to look at date in this format 
ggplot(data=subset(soilmoisture1,water_content > 0.1),aes(x=time_24hr, y=water_content,color=as.factor(plot))) + 
  geom_line()+
  facet_wrap(~site)+
  ylab('Water Content') + xlab('Date')+
  ggtitle("beanDIP 2020 Soil Moisture Data")+
  guides(col=guide_legend("Plot"))

#Below is how to make a graph that shows just one plot at one site 

ggplot(data=subset(soilmoisture, site=="Keedysville"& plot==2),aes(x=time_24hr, y=water_content)) + 
  geom_line()+
  scale_x_datetime(date_labels= "%Y-%m-%d %HH:%MM:%SS")

#Below is a graph that should have a pane for all four sites! (using new dataset soilmoisture1)
#Hm... ok, I guess it worked but the PH data is so crazy that it throws off the scale for all the rest... I need to figure out how to delete values <.1 in the dataset... 

ggplot(data=subset(soilmoisture1,drop!="drop"),aes(x=time_24hr, y=water_content,color=as.factor(plot))) + 
  geom_line()+
  scale_x_datetime(date_labels= "%Y-%m-%d %HH:%MM:%SS")+
  facet_wrap(~site)

#Attempt at not including values <.1 
#Yay, I did it! :)
#Hm... some of the Wye and PH data still looks a little weird, with the lines that go straight up to join up the the rest of the data
#I guess at PH that's because of the hole in the middle of some of the data where the chewing happened... but what's going on with Wye?
#Ohhhh... I need to somehow reintroduce the "drop" of the third plot at Wye that was off... not really sure how to include two subsets in my data? 
ggplot(data=subset(soilmoisture1,water_content > 0.1),aes(x=time_24hr, y=water_content,color=as.factor(plot))) + 
  geom_line()+
  scale_x_datetime(date_labels= "%Y-%m-%d %HH:%MM:%SS")+
  facet_wrap(~site)

#Maybe I can include the drop of plot 3 at Wye and the > .1 in the same statement... 
#Ok, yes I can, but the Wye plot still has that weird straight line?! 
ggplot(data=subset(soilmoisture1,water_content > 0.1&drop!="drop"),aes(x=time_24hr, y=water_content,color=as.factor(plot))) + 
  geom_line()+
  scale_x_datetime(date_labels= "%Y-%m-%d %HH:%MM:%SS")+
  facet_wrap(~site)




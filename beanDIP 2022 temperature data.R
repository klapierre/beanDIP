#December 12, 2022
#Sarah Alley
#beanDIP 2022 temperature data

library(tidyverse)

#Work computer
setwd("C:/Users/AlleyS/Dropbox (Smithsonian)/Work Stuff/HOBOware")

Cs <- read.csv("Temperature 2022/Ready for R/2022_Cville_sh_for_R.csv")%>%
  mutate(site="Clarksville",shielded="yes")%>%
  rename(obs=ï..obs)

Cu <- read.csv("Temperature 2022/Ready for R/2022_Cville_unsh_for_R.csv")%>%
  mutate(site="Clarksville",shielded="no") %>% 
  rename(obs=ï..obs)

Ks<- read.csv("Temperature 2022/Ready for R/2022_Kville_sh_for_R.csv")%>%
  mutate(site="Keedysville",shielded="yes") %>% 
  rename(obs=ï..obs)

Ku <- read.csv("Temperature 2022/Ready for R/2022_Kville_unsh_for_R.csv")%>%
  mutate(site="Keedysville",shielded="no") %>% 
  rename(obs=ï..obs)

PHs <- read.csv("Temperature 2022/Ready for R/2022_PopHill_sh_for_R.csv")%>%
  mutate(site="Poplar Hill",shielded="yes") %>% 
  rename(obs=ï..obs)

PHu <- read.csv("Temperature 2022/Ready for R/2022_PopHill_unsh_for_R.csv")%>%
  mutate(site="Poplar Hill",shielded="no") %>% 
  rename(obs=ï..obs)

Wu <- read.csv("Temperature 2022/Ready for R/2022_Wye_unsh_for_R.csv")%>%
  mutate(site="Wye",shielded="no") %>% 
  rename(obs=ï..obs)

Ws <- read.csv("Temperature 2022/Ready for R/2022_Wye_sh_for_R.csv")%>%
  mutate(site="Wye",shielded="yes") %>% 
  rename(obs=ï..obs)

temperature2022 <- Cs%>%rbind(Cu)%>%rbind(Ks)%>%rbind(Ku)%>%rbind(Ws)%>%
  rbind(Wu)%>%rbind(PHs)%>%rbind(PHu)%>%
  unite(time_am_pm, c(date, time, am_pm), sep=" ") %>% 
  mutate(time_24hr=as.POSIXct(strptime(time_am_pm,"%m/%d/%Y %I:%M:%S %p")))

temperatureMin <-temperature2022%>%
group_by(time_24hr, site)%>%
  summarize(temp_min=min(temp))%>%
  ungroup()

ggplot(data=temperature2022,aes(x=time_24hr, y=temp,color=as.factor(shielded))) + 
  geom_line()+
  facet_wrap(~site)+
  ylab('Temperature (F)') + xlab('Date')+
  ggtitle("beanDIP 2021 Temperatures")+
  guides(col=guide_legend("Shielded?"))



ggplot(data=temperature2022,aes(x=time_24hr, y=temp,color=as.factor(shielded))) + 
  geom_line()+
  scale_x_datetime(date_labels= "%Y-%m-%d %HH:%MM:%SS")+
  facet_wrap(~site)+
  ylab('Temperature (F)') + xlab('Date')

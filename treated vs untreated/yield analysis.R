#####################################################################
##  yield analysis.R
##  graphing and analyzing variety trial yield data from all years
##  treated/untreated pairs only
##
##  Author: Kelsey McGurrin 
#####################################################################

# working directory path (add yours if different)
setwd("~/Dropbox/bean_dip_2018-2024/field trials/data")

#### setup ####
library(tidyverse)

# cleaned data file
data<-read_csv("clean_data/clean_all_years_long.csv")
data$site<-factor(data$site,levels=c("K","C","W","PH"))


#### all varieties ####

# all points
ggplot(data,aes(x=site,y=beans_g,color=brandline_seedcoat))+
  geom_jitter()+
  facet_wrap(~year,nrow=1)

# calculate means and error bars
siteAvgDf<-data %>%
  group_by(site,year,brandline_seedcoat) %>%
  summarise(site_avg=mean(beans_g,na.rm=T),
            se=sd(beans_g,na.rm = T)/sqrt(n()))

ggplot(siteAvgDf,aes(x=site,y=site_avg,group=brandline_seedcoat,color=brandline_seedcoat))+
  geom_line(size=1.5)+
  geom_point()+
  geom_errorbar(aes(ymin=site_avg-se,
                    ymax=site_avg+se),width=0.2,size=1.5)+
  facet_wrap(~year,nrow=1)+
  labs(title = "Average Yield (beans g)")

#### AG38X8 2019-2021 ####

# all points (n= 6 plants per plot for untreated, 4 plants for treated)
AG<-data %>%
  filter(brandline=="AG38X8" & year== c(2019, 2020, 2021))

ggplot(AG,aes(x=site,y=beans_g,color=seed_treat,shape=seed_treat))+
  geom_jitter()+
  facet_wrap(~year,nrow=1)+
  labs(title = "AG38X8 Yield (beans g)")


# calculate means and error bars
siteAvgDf<-data %>%
  filter(brandline=="AG38X8" & year== c(2019,2020,2021)) %>%
  group_by(site,year,brandline,seed_treat) %>%
  summarise(site_avg=mean(beans_g,na.rm=T),
            se=sd(beans_g,na.rm = T)/sqrt(n()))

ggplot(siteAvgDf,aes(x=site,y=site_avg,group=seed_treat,linetype=seed_treat))+
  geom_line(size=1.5)+
  geom_point()+
  geom_errorbar(aes(ymin=site_avg-se,
                    ymax=site_avg+se),width=0.2,size=1.5)+
  facet_wrap(~year,nrow=1)+
  labs(title = "AG38X8 Average Yield (beans g)")



#### SH3814LL 2021 ####

# all points (n= 4 plants per plot)
SH<-data %>%
  filter(brandline=="SH3814 LL" & year== 2021)

ggplot(SH,aes(x=site,y=beans_g,color=seed_treat,shape=seed_treat))+
  geom_jitter()+
  facet_wrap(~year,nrow=1)+
  labs(title = "SH3814LL Yield (beans g)")


# calculate means and error bars
siteAvgDf<-data %>%
  filter(brandline=="SH3814 LL" & year== 2021) %>%
  group_by(site,year,brandline,seed_treat) %>%
  summarise(site_avg=mean(beans_g,na.rm=T),
            se=sd(beans_g,na.rm = T)/sqrt(n()))

ggplot(siteAvgDf,aes(x=site,y=site_avg,group=seed_treat,linetype=seed_treat))+
  geom_line(size=1.5)+
  geom_point()+
  geom_errorbar(aes(ymin=site_avg-se,
                    ymax=site_avg+se),width=0.2,size=1.5)+
  facet_wrap(~year,nrow=1)+
  labs(title = "SH3814LL Average Yield (beans g)")



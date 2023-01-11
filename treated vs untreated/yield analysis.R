#####################################################################
## yield analysis.R
## combining, graphing, and analyzing variety trial yield data from all years
##
## Author: Kelsey McGurrin 
#####################################################################

# working directory path (add yours if different)
setwd("~/Dropbox/bean_dip_2018-2024/field trials/data")

####setup####
library(tidyverse)

dat_19<-read_csv("clean_data/clean_harvest_2019.csv")
dat_19<-dat_19 %>%
  add_column(year=2019) %>%
  select(year, site, variety, plot, indiv, plant_biomass_g)

dat_20<-read_csv("clean_data/clean_harvest_2020.csv")
dat_20<-dat_20 %>%
  add_column(year=2020) %>%
  select(year, site, variety, plot, indiv, plant_biomass_g)

dat_21<-read_csv("clean_data/clean_harvest_2021.csv")
dat_21<-dat_21 %>%
  add_column(year=2021) %>%
  select(year, site, variety, plot, indiv, plant_biomass_g)

dat_22<-read_csv("clean_data/clean_harvest_2022.csv")
dat_22<-dat_22 %>%
  select(year, site, variety, plot, indiv, plant_biomass_g)

data<-rbind(dat_19,dat_20,dat_21,dat_22)
data$site<-factor(data$site,levels=c("K","C","W","PH"))

#### put in the real variety names and whether seed treated or untreated ####
data<-data %>% 
  mutate(variety_name = ifelse(
  year == 2019 & variety %in% c(67, 27),
  'AG38X8',
  ifelse(
    year == 2020 & variety %in% c(2, 5),
    'AG38X8',
    ifelse(
      year == 2021 & variety %in% c(19, 59),
      'AG38X8',
      ifelse(
        year == 2022 & variety %in% c(73),
        'AG38X8',
        ifelse(
          year == 2019 & variety == 83,
          'SH3814LL',
          ifelse(
            year == 2020 & variety == 50,
            'SH3814LL',
            ifelse(
              year == 2021 & variety %in% c(57, 58),
              'SH3814LL',
              ifelse(
                year == 2019 & variety == 31,
                'S39-G2X',
                ifelse(
                  year == 2020 & variety == 70,
                  'S39-G2X',
                  ifelse(
                    year == 2021 & variety == 82,
                    'S39-G2X',
                    ifelse(
                      year == 2019 & variety == 55,
                      '7390ET',
                      ifelse(year ==
                               2019 & variety == 32, '539XT', 'NA')
                    )
                  )
                )
              )
            )
          )
        )
      )
    )
  )
)) %>%
  mutate(seed_coat = ifelse(
    year == 2019 & variety == 67,
    'untrt',
    ifelse(
      year == 2020 & variety == 2,
      'untrt',
      ifelse(
        year == 2021 & variety %in% c(19, 57),
        'untrt',
        ifelse(year == 2022 &
                 variety == 73, 'untrt', 'trt')
      )
    )
  ))


#### peek at everything ####
ggplot(data, aes(x=seed_coat,y=plant_biomass_g,color=variety_name))+geom_point()+facet_wrap(~year)


#### subset our target variety ####
# our target variety doesn't have treated/untreated in 2022
# SH3814LL has treated/untreated only in 2021
AG<-data %>%
  filter(variety_name=="AG38X8" & year != 2022)

#overlayed line graphs across environment
#calculate means and error bars
siteAvgDf<-AG %>%
  group_by(site,year,variety_name,seed_coat) %>%
  summarise(site_avg=mean(plant_biomass_g,na.rm=T),
            se=sd(plant_biomass_g,na.rm = T)/sqrt(n()))

ggplot(siteAvgDf, aes(x=site,y=site_avg,group=seed_coat,linetype=seed_coat,color=variety_name))+geom_line(size=1.5)+geom_point(size=2)+
  geom_errorbar(aes(ymin=site_avg-se,
                    ymax=site_avg+se),width=0.2,size=1.5)+
  facet_wrap(~year,nrow=1)+
  ylab("Average Aboveground Biomass (g)")

#### subset SH3814LL ####
# SH3814LL has treated/untreated only in 2021
SH<-data %>%
  filter(variety_name == "SH3814LL" & year == 2021)

#overlayed line graphs across environment
#calculate means and error bars
siteAvgDf<-SH %>%
  group_by(site,year,variety_name,seed_coat) %>%
  summarise(site_avg=mean(plant_biomass_g,na.rm=T),
            se=sd(plant_biomass_g,na.rm = T)/sqrt(n()))

ggplot(siteAvgDf, aes(x=site,y=site_avg,group=seed_coat,linetype=seed_coat,color=variety_name))+geom_line(size=1.5)+geom_point(size=2)+
  geom_errorbar(aes(ymin=site_avg-se,
                    ymax=site_avg+se),width=0.2,size=1.5)+
  facet_wrap(~year,nrow=1)+
  ylab("Average Aboveground Biomass (g)")

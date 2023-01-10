################################################################################
##  treated vs untreated_damage data.R: Examining field soybean trial insect damage data.
##
##  Authors: Kimberly Komatsu
##  Date created: January 10, 2023
################################################################################

#packages
library(tidyverse)

#set wd
setwd('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\bean_dip_2018-2024\\field trials\\data\\clean_data')

#import data
damage2019 <- read.csv('clean_damage_2019.csv') %>% 
  mutate(pct_pucker='NA', year=2019)
damage2020 <- read.csv('clean_damage_2020.csv') %>% 
  mutate(pct_pucker='NA', year=2020)
damage2021 <- read.csv('clean_damage_2021.csv') %>% 
  mutate(year=2021)
damage2022 <- read.csv('clean_damage_2022.csv') %>% 
  select(-...18, -plot_id) %>% 
  mutate(pct_pucker='NA', year=2022) %>% 
  rename(pct_damage_stippled=stip_pct)

damage <- rbind(damage2019, damage2020, damage2021, damage2022) %>% #bind all data together
  # put in the real variety names and whether treated or untreated
  mutate(variety_name=ifelse(year==2019 & variety %in% c(67,27), 'AG38X8',
                             ifelse(year==2020 & variety %in% c(2,5), 'AG38X8',
                                    ifelse(year==2021 & variety %in% c(19,59), 'AG38X8',
                                           ifelse(year==2022 & variety %in% c(73), 'AG38X8',
          ifelse(year==2019 & variety==83, 'SH3814LL',
                 ifelse(year==2020 & variety==50, 'SH3814LL',
                        ifelse(year==2021 & variety %in% c(57,58), 'SH3814LL',
          ifelse(year==2019 & variety==31, 'S39-G2X',
                 ifelse(year==2020 & variety==70, 'S39-G2X',
                        ifelse(year==2021 & variety==82, 'S39-G2X',
          ifelse(year==2019 & variety==55, '7390ET',
          ifelse(year==2019 & variety==32, '539XT', 'NA'))))))))))))) %>% 
  mutate(seed_coat=ifelse(year==2019 & variety==67, 'untrt',
                          ifelse(year==2020 & variety==2, 'untrt',
                                 ifelse(year==2021 & variety %in% c(19,57), 'untrt',
                                        ifelse(year==2022 & variety==73, 'untrt', 'trt'))))) %>% 
  group_by(year, sampling.round, site, variety_name, seed_coat, plot) %>% 
  summarise_at(vars('density', 'chew_pct', 'pct_leaves_stippled', 'pct_damage_stippled', 'pct_yellow_leaf', 'pct_pucker'), funs(mean), na.rm=T) %>% 
  ungroup()

ggplot(data=)

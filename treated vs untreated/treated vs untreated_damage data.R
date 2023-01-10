################################################################################
##  treated vs untreated_damage data.R: Examining field soybean trial insect damage data.
##
##  Authors: Kimberly Komatsu
##  Date created: January 10, 2023
################################################################################

library(tidyverse)

setwd('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\bean_dip_2018-2024\\field trials\\data\\clean_data')

damage2019 <- read.csv('clean_damage_2019.csv') %>% 
  mutate(pct_pucker='NA')
damage2020 <- read.csv('clean_damage_2020.csv') %>% 
  mutate(pct_pucker='NA')
damage2021 <- read.csv('clean_damage_2021.csv')
damage2022 <- read.csv('clean_damage_2022.csv') %>% 
  select(-...18, -sampling.round) %>% 
  mutate(pct_pucker='NA')


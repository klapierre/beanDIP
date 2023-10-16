################################################################################
##  insects.R: Examining field soybean trial insect community data.
##
##  Authors: Kim Komatsu, Karin Burghardt
##  Date created: January 10, 2023
################################################################################

library(tidyverse)


setwd("C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\bean_dip_2018-2024\\field trials\\data\\raw_data") #Kim's path
setwd("~/Dropbox (Smithsonian)/bean_dip_2018-2024/field trials/data/raw_data") #Karin's path

#### Get insect taxa to classify to functional group ####
insects <- read.csv('beanDIP_insectdata_allyears_kjk.csv') %>% 
  select(Unique, sampling.round, Year, date, site, variety, plot, row, treatment, treatment_agg, notes, insect.names) %>% 
  separate(insect.names, into=c('names.a', 'names.b', 'names.c', 'names.d', 'names.e', 'names.f', 'names.g', 'names.h', 'names.i', 'names.j', 'names.k', 'names.l', 'names.m', 'names.n'), sep=',') %>% 
  pivot_longer(names.a:names.n, names_to='drop', values_to='names') %>% 
  filter(!is.na(names), names!='') %>% 
  separate(names, into = c("number", "taxa"), sep = "(?=[a-z +]+)(?<=[0-9])") %>% 
  left_join(read.csv('beanDIP_insectdata_allyears - InsectID_kjk3.csv')) %>%  # join functional groups
  group_by()

# Generate a list of unique taxa names to export and assign to functional groups (joined above)
insectNames <- insects %>% 
  select(taxa) %>% 
  unique()

#### Functional group responses to seed coats ####


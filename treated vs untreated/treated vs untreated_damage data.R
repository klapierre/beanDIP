################################################################################
##  treated vs untreated_damage data.R: Examining field soybean trial insect damage data.
##
##  Authors: Kimberly Komatsu
##  Date created: January 10, 2023
################################################################################

#### setup ####
#packages
library(tidyverse)

#set wd
#setwd('C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\bean_dip_2018-2024\\field trials\\data\\clean_data') # Kim path
setwd('~/Dropbox/bean_dip_2018-2024/field trials/data/clean_data') # Kelsey path

#import data
damage2019 <- read.csv('clean_damage_2019.csv') %>% 
  mutate(pct_pucker='NA', year=2019)
damage2020 <- read.csv('clean_damage_2020.csv') %>% 
  mutate(pct_pucker='NA', year=2020)
damage2021 <- read.csv('clean_damage_2021.csv') %>% 
  mutate(year=2021)
# damage2022 <- read.csv('clean_damage_2022.csv') %>% # don't include 2022 because we only have data for untreated
#   select(-...18, -plot_id) %>% 
#   mutate(pct_pucker='NA', year=2022) %>% 
#   rename(pct_damage_stippled=stip_pct)

damage <- rbind(damage2019, damage2020, damage2021) %>% #bind all data together
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
  filter(variety_name=="AG38X8" & year== c(2019,2020,2021) # only varieties and years where we have both treated/untreated
         | variety_name=="SH3814LL" & year== 2021) %>%
  mutate(stip_pct=(pct_leaves_stippled*as.numeric(pct_damage_stippled))/100,.after=chew_pct,.keep="unused") %>% #calculate stippling
  select(-c(density,mammals,insects,insect.names,pct_pucker)) %>% # drop unnecessary columns
  group_by(year, sampling.round, site, variety_name, seed_coat, plot) %>% # average rows a and b to get one value per plot
  summarise(across(chew_pct:pct_yellow_leaf, ~ mean(.x, na.rm = TRUE))) %>%
  ungroup()

# so sites display properly from west to east
damage$site<-factor(damage$site, levels = c("K", "C", "W", "PH"))

#### AG38X8 only - panels ####
AG_siteAvgDf<-damage %>%
  filter(variety_name=="AG38X8") %>%
  group_by(sampling.round,site,year,seed_coat) %>%
  summarise(site_avg_chew=mean(chew_pct,na.rm=T),
            se_avg_chew=sd(chew_pct,na.rm = T)/sqrt(n()),
            site_avg_stip=mean(stip_pct,na.rm=T),
            se_avg_stip=sd(stip_pct,na.rm=T)/sqrt(n()))

# percent chew panels
ggplot(AG_siteAvgDf, aes(x=sampling.round,y=site_avg_chew,colour=seed_coat,group=seed_coat))+
  facet_grid(rows=vars(year),cols=vars(site))+geom_point()+geom_line()+
  geom_errorbar(aes(ymin=site_avg_chew-se_avg_chew,
                    ymax=site_avg_chew+se_avg_chew),width=0.2)+
  labs(title="AG38X8",y="Avg % Chew ± SE")

# percent stippling panels
ggplot(AG_siteAvgDf, aes(x=sampling.round,y=site_avg_stip,colour=seed_coat,group=seed_coat))+
  facet_grid(rows=vars(year),cols=vars(site))+geom_point()+geom_line()+
  geom_errorbar(aes(ymin=site_avg_stip-se_avg_stip,
                    ymax=site_avg_stip+se_avg_stip),width=0.2)+
  labs(title="AG38X8",y="Avg % Stippling ± SE")

#### SH only - panels ####

# Keedysville roundup damage - skipped survey sampling.round 3, round 4 only measured surviving stalks
# mean density went from ~ 16 stems to ~ 6 stems

SH_siteAvgDf<-damage %>%
  filter(variety_name=="SH3814LL") %>%
  group_by(sampling.round,site,seed_coat) %>%
  summarise(site_avg_chew=mean(chew_pct,na.rm=T),
            se_avg_chew=sd(chew_pct,na.rm = T)/sqrt(n()),
            site_avg_stip=mean(stip_pct,na.rm=T),
            se_avg_stip=sd(stip_pct,na.rm=T)/sqrt(n()))

# percent chew panels
ggplot(SH_siteAvgDf, aes(x=sampling.round,y=site_avg_chew,colour=seed_coat,group=seed_coat))+
  facet_wrap(~site,nrow=1)+geom_point()+geom_line()+
  geom_errorbar(aes(ymin=site_avg_chew-se_avg_chew,
                    ymax=site_avg_chew+se_avg_chew),width=0.2)+
  labs(title="SH3814LL",y="Avg % Chew ± SE")

# percent stippling panels
ggplot(SH_siteAvgDf, aes(x=sampling.round,y=site_avg_stip,colour=seed_coat,group=seed_coat))+
  facet_wrap(~site,nrow=1)+geom_point()+geom_line()+
  geom_errorbar(aes(ymin=site_avg_stip-se_avg_stip,
                    ymax=site_avg_stip+se_avg_stip),width=0.2)+
  labs(title="SH3814LL",y="Avg % Stippling ± SE")

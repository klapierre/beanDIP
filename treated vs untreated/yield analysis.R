#####################################################################
## yield analysis.R
## graphing and analyzing variety trial yield data from all years
##
## Author: Kelsey McGurrin 
## Date created: Jan 10, 2023
#####################################################################

# paths for working directories (add yours)
setwd("~/Documents/beanDIP/non_git")

####setup####
library(tidyverse)
library(gridExtra)
library(GGally)

data<-read_csv("clean_all_years_long.csv")
data$site<-factor(data$site,levels=c("K","C","W","PH"))


####seed coat comparison####
coats<-data %>%
  filter(brandline == "AG38X8")

#overlayed line graphs across environment
#calculate means and error bars
siteAvgDf<-coats %>%
  group_by(site,year,brandline,seed_treat) %>%
  summarise(site_avg=mean(beans_g,na.rm=T),
            se=sd(beans_g,na.rm = T)/sqrt(n()))

ggplot(siteAvgDf, aes(x=site,y=site_avg,group=seed_treat,linetype=seed_treat,color=brandline))+geom_line(size=1.5)+geom_point(size=2)+
  geom_errorbar(aes(ymin=site_avg-se,
                    ymax=site_avg+se),width=0.2,size=1.5)+
  facet_wrap(~year,nrow=1)+
  ylab("Site Average Beans (g)")

####variety comparison####
varieties<-data %>%
  filter(seed_treat=="treated")

#overlayed line graphs across environment
#calculate means and error bars
siteAvgDf<-varieties %>%
  group_by(site,year,brandline) %>%
  summarise(site_avg=mean(leaf_trichomes_count,na.rm=T),
            se=sd(leaf_trichomes_count,na.rm = T)/sqrt(n()))

ggplot(siteAvgDf, aes(x=site,y=site_avg,group=brandline,color=brandline))+geom_line(size=1.5)+geom_point(size=2)+
  geom_errorbar(aes(ymin=site_avg-se,
                    ymax=site_avg+se),width=0.2,size=1.5)+
  facet_wrap(~year,nrow=1)+
  ylab("site_avg_trichomes")

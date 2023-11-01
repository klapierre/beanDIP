################################################################################
##  beanDIP_monocultureGreenshouse.R: Examining greenhouse trial with monocultures of rhizobial strains.
##
##  Authors: Kim Komatsu
##  Date created: October 31, 2023
################################################################################

library(codyn)
library(nlme)
library(emmeans)
library(GLMMadaptive)
library(performance)
library(tidyverse)


setwd("C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\bean_dip_2018-2024\\greenhouse trials\\UMD Summer 2022 Greenhouse_monoculture") #Kim's path


theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=20),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=20),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_text(size=16), legend.text=element_text(size=16))

###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))
barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}


options(contrasts=c('contr.sum','contr.poly'))

#not in function
`%!in%` = Negate(`%in%`)



###########################################################################
###########################################################################

#### Read in data ####

# treatment data
trt <- read.csv('Datasheets\\Tray and Pot Assignments\\beanDIP_GH2022_Reference_Analysis.csv')

# leaf traits
traits <- trt %>% 
  left_join(read.csv('Datasheets\\Leaf Traits\\Leaf_traits_all.csv')) %>% 
  # left_join(read.csv('Datasheets\\Leaf Traits\\LeafTraits_trichome data_clean.csv')) %>% #not enough trichome data for insect added plants
  # select(-notes) %>% 
  left_join(read.csv('Datasheets\\Feeding Trial\\Feeding Trial_GH2022_caterpillar_weights_clean.csv')) %>% 
  select(-notes) %>% 
  left_join(read.csv('Datasheets\\Feeding Trial\\Feeding Trial_Herbivory Data_LeafByte_clean.csv')) %>% 
  filter(!is.na(initital_weight..g.))


#### Exploratory graphs ####

ggplot(data=subset(traits, Thickness<3), aes(x=Thickness, y=relative_growth_rate)) +
  geom_point() +
  geom_smooth(method='lm')

ggplot(data=subset(traits, Thickness<3), aes(x=Thickness, y=percent_consumed)) +
  geom_point() +
  geom_smooth(method='lm')

ggplot(data=traits, aes(x=Relative.Chlorophyll, y=relative_growth_rate)) +
  geom_point() +
  geom_smooth(method='lm')

ggplot(data=traits, aes(x=Relative.Chlorophyll, y=percent_consumed)) +
  geom_point() +
  geom_smooth(method='lm')

ggplot(data=subset(traits, specific_leaf_area<2.5), aes(x=specific_leaf_area, y=relative_growth_rate)) +
  geom_point() +
  geom_smooth(method='lm')

ggplot(data=subset(traits, specific_leaf_area<2.5), aes(x=specific_leaf_area, y=percent_consumed)) +
  geom_point() +
  geom_smooth(method='lm')
#SLA look into more

ggplot(data=traits, aes(x=leaf_dry_matter_content, y=relative_growth_rate)) +
  geom_point() +
  geom_smooth(method='lm')
#LDMC look into more

ggplot(data=traits, aes(x=leaf_dry_matter_content, y=percent_consumed)) +
  geom_point() +
  geom_smooth(method='lm')
#LDMC look into more

ggplot(data=traits, aes(x=nodule_count, y=relative_growth_rate)) +
  geom_point() +
  geom_smooth(method='lm')

ggplot(data=traits, aes(x=nodule_count, y=percent_consumed)) +
  geom_point() +
  geom_smooth(method='lm')

ggplot(data=traits, aes(x=shoot_weight_g, y=relative_growth_rate)) +
  geom_point() +
  geom_smooth(method='lm')
#shoot weight look into more

ggplot(data=traits, aes(x=shoot_weight_g, y=percent_consumed)) +
  geom_point() +
  geom_smooth(method='lm')
#shoot weight look into more

ggplot(data=traits, aes(x=root_weight_g, y=relative_growth_rate)) +
  geom_point() +
  geom_smooth(method='lm')
#root weight look into more

ggplot(data=traits, aes(x=root_weight_g, y=percent_consumed)) +
  geom_point() +
  geom_smooth(method='lm')
#root weight look into more

ggplot(data=traits, aes(x=as.numeric(plant_height), y=relative_growth_rate)) +
  geom_point() +
  geom_smooth(method='lm')
#height look into more

ggplot(data=traits, aes(x=as.numeric(plant_height), y=percent_consumed)) +
  geom_point() +
  geom_smooth(method='lm')

ggplot(data=traits, aes(x=PhiNPQ, y=relative_growth_rate)) +
  geom_point() +
  geom_smooth(method='lm')

ggplot(data=traits, aes(x=PhiNPQ, y=percent_consumed)) +
  geom_point() +
  geom_smooth(method='lm')

ggplot(data=traits, aes(x=percent_consumed, y=relative_growth_rate)) +
  geom_point() +
  geom_smooth(method='lm')
#percent consumed tightly linked to relative growth rate

#### Statistics ####

hist(traits$relative_growth_rate)


# #Does strain identity matter for leaf damage?
# summary(strainModel <- lme(percent_consumed ~ as.factor(strain_num)*water_regimen,
#                               data=traits,
#                               random=~1|bench,
#                               na.action=na.omit))
# anova.lme(strainModel, type='sequential') #interaction among selected from, water regimen, and strain ID
# emmeans(strainModel, ~strain_num*water_regimen, adjust="tukey") 
# 
# #graph
# ggplot(data=barGraphStats(data=subset(traits, !is.na(percent_consumed)), variable="percent_consumed", byFactorNames=c("strain_num", 'water_regimen')), aes(x=strain_num, y=mean, color=water_regimen)) +
#   geom_point(stat='identity', position=position_dodge(0.9)) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9))
# ### if strain came from watered, then plant was consumed less in drought and more in watered; if strain came from drought, then plant was consumed intermediate amt in both watering treatments


#Does strain origin matter for insect growth?
summary(insectGrowthModel <- lme(relative_growth_rate ~ selected_from*water_regimen*strain_num,
                            data=traits,
                            random=~1|bench,
                            na.action=na.omit))
anova.lme(insectGrowthModel, type='sequential') 
emmeans(insectGrowthModel, ~selected_from*water_regimen, adjust="tukey") 

#graph
ggplot(data=barGraphStats(data=subset(traits, !is.na(relative_growth_rate)), variable="relative_growth_rate", byFactorNames=c("selected_from", 'water_regimen')), aes(x=selected_from, y=mean, fill=water_regimen)) +
  geom_bar(stat='identity', position=position_dodge(0.9)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9))


#Does strain origin matter for leaf damage?
summary(herbivoryModel <- lme(percent_consumed ~ selected_from*water_regimen*strain_num,
                            data=traits,
                            random=~1|bench,
                            na.action=na.omit))
anova.lme(herbivoryModel, type='sequential') #interaction among selected from, water regimen, and strain ID
emmeans(herbivoryModel, ~selected_from*water_regimen, adjust="tukey") 

#graph
ggplot(data=barGraphStats(data=subset(traits, !is.na(percent_consumed)), variable="percent_consumed", byFactorNames=c("selected_from", 'water_regimen')), aes(x=selected_from, y=mean, fill=water_regimen)) +
  geom_bar(stat='identity', position=position_dodge(0.9)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  ylab('Leaf Area Consumed (%)') + xlab('Original Treatment') +
  scale_x_discrete(breaks=c('ambient','Drought'), labels=c("Control", "Drought")) +
  scale_fill_manual(values=c('#f37735', '#0091FF')) +
  theme(legend.position='none')
### if strain came from watered, then plant was consumed less in drought and more in watered; if strain came from drought, then plant was consumed intermediate amt in both watering treatments

#do an RDA of traits with points colored by caterpillar growth rate
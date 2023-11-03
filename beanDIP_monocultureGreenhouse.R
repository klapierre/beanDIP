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
library(vegan)
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
  filter(!is.na(initital_weight..g.)) %>% 
  mutate(strain_ID=ifelse(selected_from=='ambient', paste('A', strain_num, sep=''),
                                                    paste('D', strain_num, sep='')))

#average traits
traitsAvg <- traits %>% 
  mutate(plant_height=as.numeric(plant_height)) %>% 
  select(-ECSt.mAU, -gH., -vH., -date_initial, -date_final, -Date..year.month.day., -Time, -Sample.Number) %>% 
  group_by(strain_num, strain_ID, selected_from, water_regimen) %>% 
  summarise_at(vars(Ambient.Humidity:Scale.Length..cm.), mean, na.rm=T) %>% 
  ungroup()

# #Trichomes
# hist(traits$trichomes_total)
# 
# summary(strainModel <- lme(trichomes_total ~ as.factor(strain_num)*water_regimen,
#                            data=traits,
#                            random=~1|bench,
#                            na.action=na.omit))
# anova.lme(strainModel, type='sequential') #interaction among selected from, water regimen, and strain ID
# emmeans(strainModel, ~strain_num*water_regimen, adjust="tukey")
# 
# #graph
# ggplot(data=barGraphStats(data=subset(traits, !is.na(trichomes_total)), variable="trichomes_total", byFactorNames=c("strain_num", 'water_regimen')), aes(x=strain_num, y=mean, color=water_regimen)) +
#   geom_point(size=3, position=position_dodge(0.9)) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
#   ylab('Leaf Area Consumed (%)') + xlab('Rhizobial Strain') +
#   scale_x_continuous(breaks=seq(1,24,1)) +
#   scale_color_manual(values=c('#f37735', '#0091FF')) +
#   theme(legend.position='none')



ggplot(data=barGraphStats(data=traitsAvg, variable="specific_leaf_area", byFactorNames=c('strain_ID', 'water_regimen')), 
       aes(x=strain_ID, y=mean, color=water_regimen)) +
  geom_point(size=3, position=position_dodge(0.9)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  ylab('') + xlab('Strain ID') +
  scale_x_discrete(breaks=c('ambient','Drought'), labels=c("Control", "Drought")) +
  scale_color_manual(values=c('#f37735', '#0091FF')) +
  theme(legend.position='none')

ggplot(data=subset(traits, leaf_dry_matter_content<3), aes(x=leaf_dry_matter_content, y=relative_growth_rate)) +
  geom_point(size=2.5) +
  geom_smooth(method='lm', se=F, color='black')


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
  geom_point(aes(color=as.factor(strain_num))) +
  geom_smooth(method='lm', se=F) #+
  # scale_color_manual(values=c('#f37735', '#0091FF'))
#LDMC look into more

ggplot(data=subset(traitsAvg, strain_ID %!in% c('D8', 'D19', 'A13', 'A1')), aes(x=leaf_dry_matter_content, y=relative_growth_rate, label=as.factor(strain_ID))) +
  geom_point(aes(color=water_regimen, shape=selected_from)) +
  geom_smooth(method='lm', se=F) +
  geom_text(hjust=0, vjust=0)

ggplot(data=subset(traitsAvg, strain_ID %!in% c('D8', 'D19', 'A13', 'A1')), aes(x=water_regimen, y=percent_consumed, label=as.factor(strain_ID))) +
  geom_point(aes(color=selected_from, shape=selected_from)) +
  geom_smooth(method='lm', se=F) +
  geom_text(hjust=0, vjust=0) +
  facet_wrap(~strain_ID)

ggplot(data=barGraphStats(data=subset(traits, strain_ID %!in% c('D8', 'D19', 'A13', 'A1')), variable="relative_growth_rate", byFactorNames=c("strain_ID", 'water_regimen')), aes(x=strain_ID, y=mean, color=water_regimen)) +
  geom_point(size=5, position=position_dodge(0.9)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  ylab('Caterpillar Relative Growth Rate') + xlab('Strain Identity') +
  scale_color_manual(values=c('#f37735', '#0091FF')) +
  scale_x_discrete(limits=c('A11','A5','A6','A12','A15','A16','A4','A18','A20',
                            'D21','D2','D3','D7','D9','D14','D23','D22','D17','D24')) +
  theme(legend.position='none')
  

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

ggplot(data=traits, aes(x=as.numeric(plant_height), y=relative_growth_rate, color=as.factor(strain_num))) +
  geom_point() +
  geom_smooth(method='lm', se=F) #+
  # scale_color_manual(values=c('#f37735', '#0091FF'))
#height look into more

ggplot(data=traits, aes(x=strain_num, y=as.numeric(plant_height), color=as.factor(strain_num))) +
  geom_point() #+
  # geom_smooth(method='lm', se=F) #+
# scale_color_manual(values=c('#f37735', '#0091FF'))
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



#### Resppnse of consumers to strains ####

hist(traits$relative_growth_rate)
hist(sqrt(traits$percent_consumed))

#Does strain identity matter for leaf damage?
summary(strainModel <- lme(sqrt(percent_consumed) ~ as.factor(strain_num)*water_regimen,
                              data=traits,
                              random=~1|bench,
                              na.action=na.omit))
anova.lme(strainModel, type='sequential') #interaction among selected from, water regimen, and strain ID
emmeans(strainModel, ~strain_num*water_regimen, adjust="tukey")

#graph
ggplot(data=barGraphStats(data=subset(traits, !is.na(percent_consumed)), variable="percent_consumed", byFactorNames=c("strain_num", 'water_regimen')), aes(x=strain_num, y=mean, color=water_regimen)) +
  geom_point(size=3, position=position_dodge(0.9)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  ylab('Leaf Area Consumed (%)') + xlab('Rhizobial Strain') +
  scale_x_continuous(breaks=seq(1,24,1)) +
  scale_color_manual(values=c('#f37735', '#0091FF')) +
  theme(legend.position='none')


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
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  ylab('Caterpillar Relative Growth Rate') + xlab('Original Treatment') +
  scale_x_discrete(breaks=c('ambient','Drought'), labels=c("Control", "Drought")) +
  scale_fill_manual(values=c('#f37735', '#0091FF')) +
  theme(legend.position='none')


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
  ylab('') + xlab('Original Treatment') +
  scale_x_discrete(breaks=c('ambient','Drought'), labels=c("Control", "Drought")) +
  scale_fill_manual(values=c('#f37735', '#0091FF')) +
  theme(legend.position='none')
### if strain came from watered, then plant was consumed less in drought and more in watered; if strain came from drought, then plant was consumed intermediate amt in both watering treatments



#### RDA of traits ####

#histograms
hist((traits$Phi2)^3)
hist(log(traits$PhiNPQ))
hist((traits$Relative.Chlorophyll))
hist(log(traits$Thickness))
hist((traits$leaves_num))
hist((as.numeric(traits$plant_height)))
hist(log(traits$specific_leaf_area))
hist(log(traits$leaf_dry_matter_content))
hist((traits$nodule_count))
hist((traits$shoot_weight_g))
hist((traits$root_weight_g))


#by original and current water treatments
rdaData <- traitsAvg %>% 
  mutate(plant_height=as.numeric(plant_height)) %>% 
  select(#bench,pot_id,
         water_regimen, selected_from, strain_ID, 
         Phi2, PhiNPQ, Relative.Chlorophyll, Thickness, leaves_num, plant_height, specific_leaf_area, leaf_dry_matter_content, 
         nodule_count, shoot_weight_g, root_weight_g, percent_consumed, relative_growth_rate) %>% 
  na.omit() %>% 
  mutate(Phi2=Phi2^3,
         PhiNPQ=log(PhiNPQ),
         Thickness=log(Thickness),
         specific_leaf_area=log(specific_leaf_area),
         leaf_dry_matter_content=log(leaf_dry_matter_content)) %>%
  mutate(Phi2=scale(Phi2),
         PhiNPQ=scale(PhiNPQ),
         Relative.Chlorophyll=scale(Relative.Chlorophyll),
         Thickness=scale(Thickness),
         leaves_num=scale(leaves_num),
         plant_height=scale(plant_height),
         specific_leaf_area=scale(specific_leaf_area),
         leaf_dry_matter_content=scale(leaf_dry_matter_content),
         nodule_count=scale(nodule_count),
         shoot_weight_g=scale(shoot_weight_g),
         root_weight_g=scale(root_weight_g))

#treatments
explanatoryData1 <- rdaData[,c('water_regimen','selected_from','strain_ID')]
#herbivory
explanatoryData2 <- rdaData[,c('percent_consumed','relative_growth_rate','strain_ID')]
#both
explanatoryData3 <- rdaData[,c('water_regimen','selected_from','strain_ID','percent_consumed','relative_growth_rate')]

#all traits
responseData1 <- rdaData[,c('Phi2','PhiNPQ','Relative.Chlorophyll','Thickness','leaves_num','plant_height',
                           'specific_leaf_area','leaf_dry_matter_content','shoot_weight_g','root_weight_g')]
#subset of traits
responseData2 <- rdaData[,c('plant_height','specific_leaf_area','leaf_dry_matter_content','shoot_weight_g','root_weight_g')]

summary(rdaModel <- rda(responseData1, explanatoryData3, na.action='na.exclude'))
screeplot(rdaModel)
plot(rdaModel)
variableScores <- scores(rdaModel, display='species')
sampleScores <- scores(rdaModel, display='sites')

rdaResults <- cbind(rdaData, sampleScores) %>% 
  mutate(combo_trt=paste(selected_from, water_regimen, sep='::'))

ggplot(data=subset(rdaResults, RDA1>-2), aes(x=RDA1, y=RDA2, color=water_regimen, shape=selected_from, label=as.factor(strain_ID))) +
  geom_point(size=4) +
  geom_text(hjust=0, vjust=0) +
  scale_color_manual(values=c('#f37735', '#0091FF')) +
  geom_segment(aes(x = 0, y = 0, xend = -0.53282055, yend = -0.4170669), arrow = arrow(length = unit(0.5, "cm")), color='black') + #Phi2
  geom_segment(aes(x = 0, y = 0, xend = -0.55942330, yend = 0.3589058), arrow = arrow(length = unit(0.5, "cm")), color='black') + #PhiNPQ
  geom_segment(aes(x = 0, y = 0, xend = 0.0068049, yend = 0.3017377), arrow = arrow(length = unit(0.5, "cm")), color='black') + #relative chlorophyll
  geom_segment(aes(x = 0, y = 0, xend = 0.21778685, yend = 0.3262338), arrow = arrow(length = unit(0.5, "cm")), color='black') + #thickness
  geom_segment(aes(x = 0, y = 0, xend = -0.92432598, yend = -0.3778270), arrow = arrow(length = unit(0.5, "cm")), color='black') + #num leaves
  geom_segment(aes(x = 0, y = 0, xend = -1.14559406, yend = -0.0577383), arrow = arrow(length = unit(0.5, "cm")), color='black') + #plant height
  geom_segment(aes(x = 0, y = 0, xend = 0.06758576, yend = -0.8175265), arrow = arrow(length = unit(0.5, "cm")), color='black') + #SLA
  geom_segment(aes(x = 0, y = 0, xend = 0.28780961, yend = -0.7295169), arrow = arrow(length = unit(0.5, "cm")), color='black') + #LDMC
  geom_segment(aes(x = 0, y = 0, xend = 0.50997683, yend = 0.2703377), arrow = arrow(length = unit(0.5, "cm")), color='black') + #shoot weight
  geom_segment(aes(x = 0, y = 0, xend = 0.41451056, yend = 0.1492567), arrow = arrow(length = unit(0.5, "cm")), color='black') + #root weight
  xlab('RDA1 (27.9%)') + ylab('RDA2 (16.7%)')

ggplot(data=subset(rdaResults, RDA1>-2), aes(x=RDA1, y=relative_growth_rate)) +
  geom_point(size=3, aes(color=water_regimen, shape=selected_from)) +
  geom_smooth(method='lm', se=F, color='black') +
  scale_color_manual(values=c('#f37735', '#0091FF')) +
  xlab('RDA1 (27.9%)') + ylab('Caterpillar Relative Growth Rate')

ggplot(data=rdaResults, aes(x=RDA2, y=relative_growth_rate)) +
  geom_point(size=3, aes(color=water_regimen, shape=selected_from)) +
  geom_smooth(method='lm', se=F, color='black') +
  scale_color_manual(values=c('#f37735', '#0091FF')) +
  xlab('RDA2 (16.7%)') + ylab('Caterpillar Relative Growth Rate')

summary(lm(relative_growth_rate ~ RDA1 + RDA2, data=subset(rdaResults, RDA1>-2))) #R2=0.073



#### Response of individual traits ####

#SLA
summary(slaModel <- lme(specific_leaf_area ~ selected_from*water_regimen,
                              data=rdaData,
                              random=~1|bench,
                              na.action=na.omit))
anova.lme(slaModel, type='sequential') #interaction among selected from, water regimen, and strain ID
emmeans(slaModel, ~selected_from*water_regimen, adjust="tukey") 

ggplot(data=barGraphStats(data=rdaResults, variable="specific_leaf_area", byFactorNames=c("selected_from", 'water_regimen')), 
       aes(x=selected_from, y=mean, fill=water_regimen)) +
  geom_bar(stat='identity', position=position_dodge(0.9)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  ylab('') + xlab('Original Treatment') +
  scale_x_discrete(breaks=c('ambient','Drought'), labels=c("Control", "Drought")) +
  scale_fill_manual(values=c('#f37735', '#0091FF')) +
  theme(legend.position='none')

#LDMC
summary(ldmcModel <- lme(leaf_dry_matter_content ~ selected_from*water_regimen,
                        data=rdaData,
                        random=~1|bench,
                        na.action=na.omit))
anova.lme(ldmcModel, type='sequential') #interaction among selected from, water regimen, and strain ID
emmeans(ldmcModel, ~selected_from*water_regimen, adjust="tukey") 

ggplot(data=barGraphStats(data=rdaResults, variable="leaf_dry_matter_content", byFactorNames=c("selected_from", 'water_regimen')), 
       aes(x=selected_from, y=mean, fill=water_regimen)) +
  geom_bar(stat='identity', position=position_dodge(0.9)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  ylab('') + xlab('Original Treatment') +
  scale_x_discrete(breaks=c('ambient','Drought'), labels=c("Control", "Drought")) +
  scale_fill_manual(values=c('#f37735', '#0091FF')) +
  theme(legend.position='none')

#shoot weight
summary(shootModel <- lme(shoot_weight_g ~ selected_from*water_regimen,
                        data=rdaData,
                        random=~1|bench,
                        na.action=na.omit))
anova.lme(shootModel, type='sequential') #interaction among selected from, water regimen, and strain ID
emmeans(shootModel, ~selected_from*water_regimen, adjust="tukey") 

ggplot(data=barGraphStats(data=rdaResults, variable="shoot_weight_g", byFactorNames=c("selected_from", 'water_regimen')), 
       aes(x=selected_from, y=mean, fill=water_regimen)) +
  geom_bar(stat='identity', position=position_dodge(0.9)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  ylab('') + xlab('Original Treatment') +
  scale_x_discrete(breaks=c('ambient','Drought'), labels=c("Control", "Drought")) +
  scale_fill_manual(values=c('#f37735', '#0091FF')) +
  theme(legend.position='none')

#root weight
summary(rootModel <- lme(root_weight_g ~ selected_from*water_regimen,
                        data=rdaData,
                        random=~1|bench,
                        na.action=na.omit))
anova.lme(rootModel, type='sequential') #interaction among selected from, water regimen, and strain ID
emmeans(rootModel, ~selected_from*water_regimen, adjust="tukey") 

ggplot(data=barGraphStats(data=rdaResults, variable="root_weight_g", byFactorNames=c("selected_from", 'water_regimen')), 
       aes(x=selected_from, y=mean, fill=water_regimen)) +
  geom_bar(stat='identity', position=position_dodge(0.9)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  ylab('') + xlab('Original Treatment') +
  scale_x_discrete(breaks=c('ambient','Drought'), labels=c("Control", "Drought")) +
  scale_fill_manual(values=c('#f37735', '#0091FF')) +
  theme(legend.position='none')

#height
summary(heightModel <- lme(plant_height ~ selected_from*water_regimen,
                        data=rdaData,
                        random=~1|bench,
                        na.action=na.omit))
anova.lme(heightModel, type='sequential') #interaction among selected from, water regimen, and strain ID
emmeans(heightModel, ~selected_from*water_regimen, adjust="tukey") 

ggplot(data=barGraphStats(data=rdaResults, variable="plant_height", byFactorNames=c("selected_from", 'water_regimen')), 
       aes(x=selected_from, y=mean, fill=water_regimen)) +
  geom_bar(stat='identity', position=position_dodge(0.9)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  ylab('') + xlab('Original Treatment') +
  scale_x_discrete(breaks=c('ambient','Drought'), labels=c("Control", "Drought")) +
  scale_fill_manual(values=c('#f37735', '#0091FF')) +
  theme(legend.position='none')

################################################################################
##  insects.R: Examining field soybean trial insect community data.
##
##  Authors: Kim Komatsu, Karin Burghardt
##  Date created: January 10, 2023
################################################################################

library(codyn)
library(nlme)
library(emmeans)
library(GLMMadaptive)
library(performance)
library(tidyverse)


setwd("C:\\Users\\kjkomatsu\\Dropbox (Smithsonian)\\bean_dip_2018-2024\\field trials\\data\\raw_data") #Kim's path
# setwd("~/Dropbox (Smithsonian)/bean_dip_2018-2024/field trials/data/raw_data") #Karin's path


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

#### Get insect taxa to classify to functional group ####
insects <- read.csv('beanDIP_insectdata_allyears_kjk.csv') %>% 
  select(Unique, sampling.round, Year, date, site, variety, plot, row, treatment, treatment_agg, notes, insect.names) %>% 
  separate(insect.names, into=c('names.a', 'names.b', 'names.c', 'names.d', 'names.e', 'names.f', 'names.g', 'names.h', 'names.i', 'names.j', 'names.k', 'names.l', 'names.m', 'names.n'), sep=',') %>% 
  pivot_longer(names.a:names.n, names_to='drop', values_to='names') %>% 
  filter(!is.na(names), names!='') %>% 
  separate(names, into = c("number", "taxa"), sep = "(?=[a-z +]+)(?<=[0-9])") %>% 
  mutate(number=as.integer(number)) %>% 
  left_join(read.csv('beanDIP_insectdata_allyears - InsectID_kjk3.csv')) %>%  # join functional groups
  filter(functional_group!='drop') %>% 
  rename(unique=Unique, sampling_round=sampling.round, year=Year) %>% 
  group_by(unique, sampling_round, year, date, site, variety, plot, row, treatment, treatment_agg, functional_group) %>% 
  summarise(count=sum(number)) %>% 
  ungroup() %>% 
  pivot_wider(names_from=functional_group, values_from=count, values_fill=0) %>% #fill with 0 when the functional group was not present in a plot
  pivot_longer(sucking:none, names_to='functional_group', values_to='count')

# # Generate a list of unique taxa names to export and assign to functional groups (joined above)
# insectNames <- insects %>% 
#   select(taxa) %>% 
#   unique()

#### Check functional group count normality ####
with(subset(insects, functional_group=='chewing'), hist(((count))))
with(subset(insects, functional_group=='sucking'), hist((count)))
with(subset(insects, functional_group=='multiple'), hist((count))) #just ants #too few to do stats
with(subset(insects, functional_group=='pollinator'), hist((count)))
with(subset(insects, functional_group=='predator'), hist((count)))
with(subset(insects, functional_group=='detritivore'), hist((count))) #too few to do stats
with(subset(insects, functional_group=='borer'), hist((count))) #too few to do stats
with(subset(insects, functional_group=='seed'), hist((count))) #too few to do stats
with(subset(insects, functional_group=='other'), hist((count))) #too few to do stats
with(subset(insects, functional_group=='unknown'), hist((count))) #too few to do stats


#### ANOVA functional group responses to seed coats ####

#chewing herbivores
summary(chewingModel <- lme(count ~ treatment_agg,
                            data=subset(insects, functional_group=='chewing' & year!=2022 & treatment %in% c('untreated_focal', 'treated_focal') & sampling_round==4),
                            random=~1|site/row/plot,
                            weights=varIdent(),
                            correlation=corCompSymm(form=~year/sampling_round|site/row/plot),
                            control=lmeControl(returnObject=T)))
anova.lme(chewingModel, type='sequential') 
emmeans(chewingModel, ~treatment_agg, adjust="tukey") 

# #zero-inflated - chewing herbivores
# summary(chewingModel <- mixed_model(count ~ treatment_agg,
#                                     data=subset(insects, functional_group=='chewing' & year!=2022 & treatment %in% c('untreated_focal', 'treated_focal')),
#                                     random=~1|site/row/plot,
#                                     family=zi.poisson(),
#                                     zi_fixed=~treatment_agg #,
#                                     # correlation=corCompSymm(form=~year/sampling_round|site/row/plot),
#                                     # control=lmeControl(returnObject=T)
#         ))
# anova.lme(chewingModel, type='sequential') 
# emmeans(chewingModel, ~treatment_agg, adjust="tukey") 


#sucking herbivores
summary(suckingModel <- lme(count ~ treatment_agg,
                            data=subset(insects, functional_group=='sucking' & year!=2022 & treatment %in% c('untreated_focal', 'treated_focal')),
                            random=~1|site/row/plot,
                            weights=varIdent(),
                            correlation=corCompSymm(form=~year/sampling_round|site/row/plot),
                            control=lmeControl(returnObject=T)))
anova.lme(suckingModel, type='sequential') #marginal more in untreated
emmeans(suckingModel, ~treatment_agg, adjust="tukey") 


#all herbivores
summary(herbivoreModel <- lme(count ~ treatment_agg,
                              data=subset(insects, functional_group %in% c('sucking', 'chewing') & year!=2022 & treatment %in% c('untreated_focal', 'treated_focal')),
                              random=~1|site/row/plot,
                              weights=varIdent(),
                              correlation=corCompSymm(form=~year/sampling_round|site/row/plot),
                              control=lmeControl(returnObject=T)))
anova.lme(herbivoreModel, type='sequential')
emmeans(herbivoreModel, ~treatment_agg, adjust="tukey") 


#pollinators
summary(pollinatorModel <- lme(count ~ treatment_agg,
                               data=subset(insects, functional_group=='pollinator' & year!=2022 & treatment %in% c('untreated_focal', 'treated_focal')),
                               random=~1|site/row/plot,
                               weights=varIdent(),
                               correlation=corCompSymm(form=~year/sampling_round|site/row/plot),
                               control=lmeControl(returnObject=T)))
anova.lme(pollinatorModel, type='sequential')
emmeans(pollinatorModel, ~treatment_agg, adjust="tukey") 


#predators
summary(predatorModel <- lme(count ~ treatment_agg,
                               data=subset(insects, functional_group=='predator' & year!=2022 & treatment %in% c('treated_focal','untreated_focal')),
                               random=~1|site/row/plot,
                               weights=varIdent(),
                               correlation=corCompSymm(form=~year/sampling_round|site/row/plot),
                               control=lmeControl(returnObject=T)))
anova.lme(predatorModel, type='sequential') 
emmeans(predatorModel, ~treatment_agg, adjust="tukey") 



#### Figures ####
#interactions
ggplot(data=subset(insects, functional_group %!in% c('borer', 'seed', 'detritivore', 'multiple', 'other', 'none', 'unknown') & treatment %in% c('treated_focal','untreated_focal')
                   & site=='C' & year!=2022),
       aes(x=sampling_round, y=count, color=treatment)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  facet_grid(rows=vars(functional_group), cols=vars(year), scales='free_y')

ggplot(data=subset(insects, functional_group %!in% c('borer', 'seed', 'detritivore', 'multiple', 'other', 'none', 'unknown') & treatment %in% c('treated_focal','untreated_focal')
                   & site=='K' & year!=2022),
       aes(x=sampling_round, y=count, color=treatment)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  facet_grid(rows=vars(functional_group), cols=vars(year), scales='free_y')

ggplot(data=subset(insects, functional_group %!in% c('borer', 'seed', 'detritivore', 'multiple', 'other', 'none', 'unknown') & treatment %in% c('treated_focal','untreated_focal')
                   & site=='W' & year!=2022),
       aes(x=sampling_round, y=count, color=treatment)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  facet_grid(rows=vars(functional_group), cols=vars(year), scales='free_y')

ggplot(data=subset(insects, functional_group %!in% c('borer', 'seed', 'detritivore', 'multiple', 'other', 'none', 'unknown') & treatment %in% c('treated_focal','untreated_focal')
                   & site=='PH' & year!=2022),
       aes(x=sampling_round, y=count, color=treatment)) +
  geom_point() +
  geom_smooth(method='lm', se=F) +
  facet_grid(rows=vars(functional_group), cols=vars(year), scales='free_y')


#no interactions
ggplot(data=barGraphStats(data=subset(insects, functional_group %!in% c('borer', 'seed', 'other', 'none', 'unknown') & treatment %in% c('treated_focal','untreated_focal') & year!=2022), 
                          variable="count", byFactorNames=c("functional_group", "treatment")),
       aes(x=treatment, y=mean, fill=treatment)) +
  geom_bar(stat='identity', position=position_dodge(0.9)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, position=position_dodge(0.9)) +
  facet_wrap(~functional_group, scales='free')
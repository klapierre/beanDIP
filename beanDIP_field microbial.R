################################################################################
##  bacterial community analysis.R: Sequence data from beanDIP field trials
##
##  Author: Kimberly Komatsu
##  Date created: August 5, 2025
################################################################################

# NOTE: ignoring soils for now (no rhizobia amplified in the soil samples, only in nodules)


#install and load package
# if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
# devtools::install_github("jbisanz/qiime2R")

library(qiime2R)
library(openxlsx)
library(performance)
library(PerformanceAnalytics)
library(nlme)
library(emmeans)
library(car)
library(vegan)
library(grid)
library(tidyverse)

'%notin%' <- negate('%in%')

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

###setting the graph look
theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_text(size=20), legend.text=element_text(size=20),
             strip.text.x=element_text(size=20), strip.text.y=element_text(size=20))


#set working directory
setwd('C:\\Users\\kjkomatsu\\Box\\RESEARCH - kjkomatsu - Komatsu Lab Research\\beanDIP\\20250708_beanDIP_fieldAmplicon_kim version')


##### Import Data #####
#sequence variants
noduleSVs <- read_qza("20250710_beanDIP_fieldNodules\\deblur_output\\deblur_table_final.qza")

#metadata
metadata <- read.table("beanDIP_metadata_final.txt", header = TRUE, sep = "\t")

# soilsPlate <- read.xlsx("soil_beanDIP_Alley_DNA_extraction_sample_layout.xlsx") %>%
#   mutate(indiv=0) %>% 
#   select(plate, well, site, block, varietal, indiv, month, year)

nodulesPlate <- read.xlsx("nodule_beanDIP_Alley_DNA_extraction_sample_layout.xlsx") %>% 
  mutate(month=8) %>% 
  select(plate, well, site, block, varietal, indiv, month, year)

plateInfo <- nodulesPlate %>% 
  filter(site!='BLANK', site!='?') %>% 
  mutate(block=as.integer(block),
         indiv=as.integer(indiv),
         year=as.integer(year),
         entry=as.integer(varietal),
         site=ifelse(site=='PH', 'P', site)) %>% 
  select(-varietal)

varietalInfo <- read.csv('C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\bean_dip_2018-2024\\field trials\\data\\clean_data\\varietal_info.csv')

metadataFull <- metadata %>% 
  left_join(plateInfo) %>% 
  # rename(SampleID=sample.id) %>%
  mutate(site=ifelse(site=='PH', 'P', site),
         # entry=as.integer(varietal),
         year=as.integer(year)) %>%
  filter(!is.na(entry)) %>% 
  select(-varietal) %>% 
  left_join(varietalInfo)

# write.table(metadataFull, 'beanDIP_metadata_final.txt', sep = "\t", row.names = FALSE)


#taxonomy
noduletaxonomy <- read_qza("20250710_beanDIP_fieldNodules\\taxa\\classification.qza")
noduletaxonomy<-parse_taxonomy(noduletaxonomy$data)
noduletaxonomy2 <- noduletaxonomy %>% 
  rownames_to_column("Feature.ID")


##### Shannon's Entropy #####
noduleshannon <- read_qza("20250710_beanDIP_fieldNodules\\diversity\\shannon_vector.qza")$data %>%
  mutate(type='nodule') %>%
  rownames_to_column("SampleID")  #creates sample name column to merge with metadata

shannonTrt <- noduleshannon %>%
  left_join(metadataFull)

# write.csv(shannonTrt, 'beanDIP_fieldNodule_shannons.csv', row.names=F)

# shannonTrt$grazing_category=factor(shannonTrt$grazing_category,levels=c('heavy', 'stable', 'destock'))

shannonNoduleModel <- lme(shannon_entropy ~ indiv,
          data=subset(shannonTrt, type=='nodule' & varietal=='AG38X8' & treated=='untreated'),
          random=~1|site/year) #no diff across individuals
anova(shannonNoduleModel)

ggplot(data=subset(shannonTrt, type=='nodule' & varietal=='AG38X8' & treated=='untreated'), aes(x=shannon_entropy)) +
  geom_density()

shannonSiteModel <- lme(shannon_entropy ~ site,
          data=subset(shannonTrt, type=='nodule'& varietal=='AG38X8' & treated=='untreated'),
          random=~1|block/year) 

anova(shannonSiteModel) #no diff across sites


##### Faith's PD #####
nodulepd <- read_qza("20250710_beanDIP_fieldNodules\\diversity\\faith_pd_vector.qza")$data %>% 
  mutate(type='nodule')

pdTrt <- nodulepd %>% 
  rownames_to_column() %>% 
  rename(SampleID=rowname) %>%  #creates sample name column to merge with metadata
  left_join(metadataFull)

# write.csv(pdTrt, 'beanDIP_fieldNodules_faithPD.csv', row.names=F)


faithpdNoduleModel <- lme(faith_pd ~ indiv,
                          data=subset(pdTrt, type=='nodule' & varietal=='AG38X8' & treated=='untreated'),
                          random=~1|site/year) #no diff across individuals
anova(faithpdNoduleModel)



##### ATV Richness #####
noduleRichness <- read_qza("20250710_beanDIP_fieldNodules\\diversity\\observed_features_vector.qza")$data %>% 
  mutate(type='nodule')

richnessTrt <- noduleRichness %>% 
  rownames_to_column() %>% 
  rename(SampleID=rowname) %>%  #creates sample name column to merge with metadata
  left_join(metadataFull)

# write.csv(richnessTrt, 'beanDIP_fieldNodules_richness.csv', row.names=F)


richnessNoduleModel <- lme(observed_features ~ indiv,
                          data=subset(richnessTrt, type=='nodule' & varietal=='AG38X8' & treated=='untreated'),
                          random=~1|site/year) #no diff across individuals
anova(richnessNoduleModel)

ggplot(data=subset(richnessTrt, type=='nodule' & varietal=='AG38X8' & treated=='untreated'), aes(x=observed_features)) +
  geom_density()

ggplot(data=subset(richnessTrt, type=='nodule' & varietal=='AG38X8' & treated=='untreated'), aes(x=site, y=observed_features)) +
  geom_boxplot()

richnessSiteModel <- lme(observed_features ~ site,
                        data=subset(richnessTrt, type=='nodule'& varietal=='AG38X8' & treated=='untreated'),
                        random=~1|block/year) 

anova(richnessSiteModel) #sig diff across sites
emmeans(richnessSiteModel, pairwise~site)


ggplot(data=subset(richnessTrt, type=='nodule'), aes(x=observed_features)) +
  geom_density()

richnessSimple <- richnessTrt %>% 
  filter(type=='nodule' & varietal=='AG38X8' & treated=='untreated') %>% 
  group_by(year, site, block, treated, varietal, entry) %>% 
  summarise(richness=mean(observed_features)) %>% 
  ungroup()


##### combine with field data #####

yieldAll <- read.csv('C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\bean_dip_2018-2024\\field trials\\data\\clean_data\\clean_all_years_long.csv') %>% 
  rename(block=plot,
         treated=seed_treat,
         varietal=brandline) %>% 
  mutate(site=ifelse(site=='PH','P',site)) %>% 
  group_by(year, site, block, varietal, treated) %>% 
  summarise_at(vars(healthy_count, beans_g, plant_biomass_g), mean, na.rm = TRUE) %>% 
  ungroup() %>% 
  left_join(richnessSimple) %>% 
  filter(varietal=='AG38X8', treated=='untreated')

ggplot(data=subset(yieldAll), aes(x=richness, y=healthy_count)) +
  geom_point() +
  geom_smooth(method='lm', se=F)

ggplot(data=subset(yieldAll), aes(x=richness, y=beans_g)) +
  geom_point() +
  geom_smooth(method='lm', se=T, color='black') +
  xlab('Mean Rhizobial ASV Richness') + ylab('Mean Bean Mass (g)') +
  scale_x_continuous(breaks = seq(0, 10, by = 1))

ggplot(data=subset(yieldAll), aes(x=richness, y=plant_biomass_g)) +
  geom_point() +
  geom_smooth(method='lm', se=T, color='black') +
  xlab('Mean Rhizobial ASV Richness') + ylab('Mean Plant Mass (g)') +
  scale_x_continuous(breaks = seq(0, 10, by = 1))



damage2019 <- read.csv('C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\bean_dip_2018-2024\\field trials\\data\\clean_data\\clean_damage_2019.csv') %>% 
  mutate(year=2019)

damage2020 <- read.csv('C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\bean_dip_2018-2024\\field trials\\data\\clean_data\\clean_damage_2020.csv')%>% 
  mutate(year=2020)

damage2021 <- read.csv('C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\bean_dip_2018-2024\\field trials\\data\\clean_data\\clean_damage_2021.csv') %>% 
  select(-pct_pucker) %>% 
  mutate(year=2021)

damage2022 <- read.csv('C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\bean_dip_2018-2024\\field trials\\data\\clean_data\\clean_damage_2022.csv') %>% 
  select(-plot_id, -`...18`) %>% 
  rename(pct_damage_stippled=stip_pct) %>% 
  mutate(year=2022)

damageAll <- rbind(damage2019,damage2020,damage2021,damage2022) %>% 
  rename(block=plot,
         entry=variety) %>% 
  mutate(site=ifelse(site=='PH','P',site)) %>% 
  group_by(sampling.round, date, site, entry, block) %>% 
  summarise_at(vars(density, chew_pct, pct_leaves_stippled, 
                    pct_damage_stippled), mean, na.rm=T) %>% 
  ungroup() %>% 
  left_join(richnessSimple) %>% 
  filter(treated=='untreated', varietal=='AG38X8')

# write.csv(damageAll, 'C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\bean_dip_2018-2024\\field trials\\data\\clean_data\\insect_count_kim.csv')

ggplot(data=damageAll, aes(x=as.factor(sampling.round), y=chew_pct)) +
  geom_boxplot() +
  facet_wrap(~year)

ggplot(data=subset(damageAll, sampling.round==4), aes(x=richness, y=chew_pct)) +
  geom_point(size=3) +
  geom_smooth(method='lm',  se=T, color='black') +
  xlab('Mean Rhizobial ASV Richness') + ylab('Percent Leaf Damage') +
  scale_x_continuous(breaks = seq(0, 10, by = 2))
#export 700x700

ggplot(data=subset(damageAll, sampling.round==4), aes(x=richness, y=pct_leaves_stippled)) +
  geom_point() +
  geom_smooth(method='lm', se=T, color='black') +
  xlab('Mean Rhizobial ASV Richness') + ylab('Percent Leaves Stippled') +
  scale_x_continuous(breaks = seq(0, 10, by = 2))

ggplot(data=subset(damageAll, sampling.round==4), aes(x=richness, y=pct_damage_stippled)) +
  geom_point() +
  geom_smooth(method='lm', se=T, color='black') +
  xlab('Mean Rhizobial ASV Richness') + ylab('Percent Stipple Damage') +
  scale_x_continuous(breaks = seq(0, 10, by = 2))




##### FROM BRENDAN'S CODE - looking at monoculture greenhouse results, only control watering 

## beanDIP 2022 Greenhouse
#10/25/2022
#Brendan Randall

library(lme4)
library(lmerTest)
library(tidyverse)

# Calling data frames
setwd("C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\bean_dip_2018-2024\\greenhouse trials\\UMD Summer 2022 Greenhouse_monoculture\\Datasheets")

reference <- read.csv('Tray and Pot Assignments\\beanDIP_2022singlestrain_reference_analysis.csv')

# Predictor variable data frames 
caterpillar_weights <- read.csv('Feeding Trial\\Feeding Trial_GH2022_caterpillar_weights_clean.csv')
herbivory <- read.csv('Feeding Trial\\Feeding Trial_Herbivory Data_LeafByte_clean.csv')
trichomes <- read.csv('Leaf Traits\\LeafTraits_trichomes_clean.csv')
sla <- read.csv('Leaf Traits\\Leaf Traits_sla_clean.csv')
photosynq <- read.csv('Leaf Traits\\Leaf Traits_photosynq_clean.csv')
# heightleaves <- read.csv('LeafTraits_heightleaves.csv')
harvest <- read.csv('Harvest\\Harvest Summer 2022 Greenhouse.csv')


#Merging data frames

weightdata <-merge(reference, caterpillar_weights, by="pot_id", all.x = FALSE, all.y = TRUE)
weightdata <- na.omit(weightdata)
weightdata

herbivorydata <- merge(reference, herbivory, by="pot_id")

trichomedata <- merge(reference, trichomes, by="pot_id")


sladata <- merge(reference, sla, by="pot_id", all.x = FALSE, all.y = FALSE)
sladata <- na.omit(sladata)

photosynqdata <- merge(reference, photosynq, by="pot_id")
photosynqdata <- na.omit(photosynqdata)

# heightleavesdata <- merge(reference, heightleaves, by ="pot_id")
# heightleavesdata <- na.omit(heightleavesdata)

harvestdata <- merge(reference, harvest, by="pot_id")


#Sorting data frames according to tray number and tray position

order.weight <- with(weightdata, order(tray_number, tray_position))
weightdata <- weightdata[order.weight, ]
#Herbivory data already in tray and position order (sample number column)

#****Basic GLMMs caterpillar weights*****

weightdata$bench <- as.character(weightdata$bench)
weightdata$block <- as.factor(weightdata$block)
weightdata$sub_block <- as.factor(weightdata$sub_block)
weightdata$tray_number <- as.factor(weightdata$tray_number)
weightdata$tray_position <- as.factor(weightdata$tray_position)
weightdata$water_regimen <- as.factor(weightdata$water_regimen)
weightdata$strain_num <- as.character(weightdata$strain_num)
weightdata$species <- as.factor(weightdata$species)


summary(weightdata)

weight_model_water <- lmer(relative_growth_rate ~ water_regimen + (1|sub_block), data=weightdata)
summary(weight_model_water)
ranef(weight_model_water)
anova(weight_model_water)
plot(weight_model_water)


tab_model(weight_model_water)
plot(weight_model_water)
hist(residuals(weight_model_water))
summary(weight_model_water)
tab_model(weight_model_water)




weight_model_strain <- lmer(relative_growth_rate ~ water_regimen + strain_num + water_regimen:strain_num + (1|sub_block), data=weightdata)
summary(weight_model_strain)
ranef(weight_model_strain)
anova(weight_model_strain)
plot_model(weight_model_strain, type = "int", mdrt.values = "meansd")



tab_model(weight_model_strain)
plot(weight_model_strain)
hist(residuals(weight_model_strain))

#* Model comparison using log-likelihood with drought treatment nested in strain_ID, vary slopes


#****Basic GLMMs herbivory****

herbivorydata$bench <- as.character(herbivorydata$bench)
herbivorydata$block <- as.factor(herbivorydata$block)
herbivorydata$sub_block <- as.factor(herbivorydata$sub_block)
herbivorydata$tray_number <- as.factor(herbivorydata$tray_number)
herbivorydata$tray_position <- as.factor(herbivorydata$tray_position)
herbivorydata$water_regimen <- as.factor(herbivorydata$water_regimen)
herbivorydata$strain_num <- as.character(herbivorydata$strain_num, levels=unique(herbivorydata$strain_num))
herbivorydata$species <- as.factor(herbivorydata$species)

##consumed leaf area
herbivory_model_water_cons <- lmer(consumed_leaf_area ~ water_regimen + (1|sub_block), data=herbivorydata)
summary(herbivory_model_water_cons)
ranef(herbivory_model_water_cons)
anova(herbivory_model_water_cons)
plot(herbivory_model_water_cons)
hist(residuals(herbivory_model_water_cons))

herbivory_model_strain_cons <- lmer(consumed_leaf_area ~ water_regimen + strain_num + water_regimen:strain_num + (1|block:sub_block), data=herbivorydata)
summary(herbivory_model_strain_cons)
ranef(herbivory_model_strain_cons)
anova(herbivory_model_strain_cons)
plot(herbivory_model_strain_cons)
hist(residuals(herbivory_model_strain_cons))

## percent consumed

herbivory_model_water_per <- lmer(percent_consumed ~ water_regimen + (1|block:sub_block), data=herbivorydata)
summary(herbivory_model_water_per)
ranef(herbivory_model_water_per)
anova(herbivory_model_water_per)
plot(herbivory_model_water_per)
hist(residuals(herbivory_model_water_per))


herbivory_model_water_per <- lmer(percent_consumed ~ water_regimen + strain_num + water_regimen:strain_num + (1|block:sub_block), data=herbivorydata)
summary(herbivory_model_water_per)
ranef(herbivory_model_water_per)
anova(herbivory_model_water_per)
plot(herbivory_model_water_per)
hist(residuals(herbivory_model_water_per))

## Total Leaf Area

herbivory_model_leaf_area <- lmer(total_leaf_area ~ water_regimen + strain_num + water_regimen:strain_num + (1|sub_block), data=herbivorydata)
summary(herbivory_model_leaf_area)
ranef(herbivory_model_leaf_area)
anova(herbivory_model_leaf_area)
plot(herbivory_model_leaf_area)
hist(residuals(herbivory_model_leaf_area))

#**** Plasticity linear models- herbivory and caterpillars****


#****Basic GLMMS SLA****

sladata$bench <- as.character(sladata$bench)
sladata$block <- as.factor(sladata$block)
sladata$sub_block <- as.factor(sladata$sub_block)
sladata$tray_number <- as.factor(sladata$tray_number)
sladata$tray_position <- as.factor(sladata$tray_position)
sladata$water_regimen <- as.factor(sladata$water_regimen)
sladata$strain_num <- as.character(sladata$strain_num, levels=unique(sladata$strain_num))
sladata$species <- as.factor(sladata$species)

slamodel <- lmer(specific_leaf_area ~ water_regimen + strain_num + water_regimen:strain_num +(1|sub_block), data=sladata)
summary(slamodel)
ranef(slamodel)
anova(slamodel)
plot(slamodel)
hist(residuals(slamodel))

#****Basic GLMMS LDMC****
ldmcmodel <- lmer(leaf_dry_matter_content ~ water_regimen +strain_num +water_regimen:strain_num +(1|sub_block), data=sladata)
summary(ldmcmodel)
ranef(ldmcmodel)
anova(ldmcmodel)
plot(ldmcmodel)
hist(residuals(ldmcmodel))

#**** Basic GLMMS- PhotosynQ ****
photosynqdata$bench <- as.character(photosynqdata$bench)
photosynqdata$block <- as.factor(photosynqdata$block)
photosynqdata$sub_block <- as.factor(photosynqdata$sub_block)
photosynqdata$tray_number <- as.factor(photosynqdata$tray_number)
photosynqdata$tray_position <- as.factor(photosynqdata$tray_position)
photosynqdata$water_regimen <- as.factor(photosynqdata$water_regimen)
photosynqdata$strain_num <- as.character(photosynqdata$strain_num, levels=unique(sladata$strain_num))
photosynqdata$species <- as.factor(photosynqdata$species)
photosynqdata$plant_height <- as.numeric(photosynqdata$plant_height)

#*Leaf Temperature Differential****#
leaftempdiffmodel <- lmer(Leaf.Temperature.Differential ~ water_regimen + strain_num + water_regimen:strain_num + (1|sub_block), data=photosynqdata)
summary(leaftempdiffmodel)
ranef(leaftempdiffmodel)
anova(leaftempdiffmodel)
plot(leaftempdiffmodel)
hist(residuals(leaftempdiffmodel))

#* NPQt *#
NPQtmodel <- lmer(NPQt ~ water_regimen + strain_num + water_regimen:strain_num + (1|sub_block), data=photosynqdata)
summary(NPQtmodel)
ranef(NPQtmodel)
anova(NPQtmodel)
plot(NPQtmodel)
hist(residuals(NPQtmodel))

#*PhiNPQ non photochemical quenching*#
PhiNPQmodel <- lmer(PhiNPQ ~ water_regimen + strain_num + water_regimen:strain_num + (1|sub_block), data=photosynqdata)
summary(PhiNPQmodel)
ranef(PhiNPQmodel)
anova(PhiNPQmodel)
plot(PhiNPQmodel)
hist(residuals(PhiNPQmodel))

PhiNPQmodel <- lmer(PhiNPQ ~ water_regimen + strain_num + water_regimen:strain_num + (1|sub_block), data=photosynqdata)

PhiNPQmodel2 <- lmer(PhiNPQ ~ water_regimen + strain_num +(1|water_regimen:strain_num) + (1|sub_block), data=photosynqdata, REML= TRUE)

PhiNPQmodel3 <- lmer(PhiNPQ ~ water_regimen + strain_num + (1|sub_block), data=photosynqdata, REML= TRUE)

summary(PhiNPQmodel2)
ranef(PhiNPQmodel2)
anova(PhiNPQmodel2)
anova(PhiNPQmodel)
anova(PhiNPQmodel2, PhiNPQmodel3, test= "Chisq", REML=TRUE)
anova(PhiNPQmodel3, PhiNPQmodel2, test= "Chisq", REML=TRUE)

#*PhiNO*#
PhiNOmodel <- lmer(PhiNO ~ water_regimen + strain_num + water_regimen:strain_num + (1|sub_block), data=photosynqdata)
summary(PhiNOmodel)
ranef(PhiNOmodel)
anova(PhiNOmodel)
plot(PhiNOmodel)
hist(residuals(PhiNOmodel))

#*Phi2*#
Phi2model <- lmer(Phi2 ~ water_regimen + strain_num + water_regimen:strain_num + (1|sub_block), data=photosynqdata)
summary(Phi2model)
ranef(Phi2model)
anova(Phi2model)
plot(Phi2model)
hist(residuals(Phi2model))

# Relative Chlorophyll

chlorophyllmodel <- lmer(Relative.Chlorophyll ~ water_regimen + strain_num + water_regimen:strain_num + (1|sub_block), data=photosynqdata)
summary(chlorophyllmodel)
ranef(chlorophyllmodel)
anova(chlorophyllmodel)
plot(chlorophyllmodel)
hist(residuals(chlorophyllmodel))

#* Leaf thickness *# 

thicknessmodel <- lmer(Thickness ~ water_regimen + strain_num + water_regimen:strain_num + (1|sub_block), data=photosynqdata)
summary(thicknessmodel)
ranef(thicknessmodel)
anova(thicknessmodel)
plot(thicknessmodel)
hist(residuals(thicknessmodel))

#* Leaf Angle *#

anglemodel <- lmer(Leaf.Angle ~ water_regimen + strain_num + water_regimen:strain_num + (1|sub_block), data=photosynqdata)
summary(anglemodel)
ranef(anglemodel)
anova(anglemodel)
plot(anglemodel)
hist(residuals(anglemodel))

#**** Basic GLMMS- Height and Number of Leaves ****


#* plant height (cm)*#- boundary error, introducing NAs by coercion for some reason. When you try to run model as is, it claims it's a factor when it should be numeric. After converting to numeric, 

heightmodel <- lmer(plant_height ~ water_regimen + strain_num + water_regimen:strain_num +(1|sub_block), data=heightleavesdata)
summary(heightmodel)
ranef(heightmodel)
anova(heightmodel)
plot(heightmodel)
hist(residuals(heightmodel))

#* number of trifoliate leaflets*# - also not working well. Boundary errors and no intercepts for random effects. Could try extracting data from the photosynQ datasheet into another datasheet. 

leafletmodel <- lmer(trifoliate_leaflet_num ~ water_regimen + strain_num + water_regimen:strain_num +(1|sub_block), data=photosynqdata)
summary(leafletmodel)
ranef(leafletmodel)
anova(leafletmodel)
plot(leafletmodel)
hist(residuals(leafletmodel))

#**** Basic GLMMs Trichomes****

trichomedata$bench <- as.character(trichomedata$bench)
trichomedata$block <- as.factor(trichomedata$block)
trichomedata$sub_block <- as.factor(trichomedata$sub_block)
trichomedata$tray_number <- as.factor(trichomedata$tray_number)
trichomedata$tray_position <- as.factor(trichomedata$tray_position)
trichomedata$water_regimen <- as.factor(trichomedata$water_regimen)
trichomedata$strain_num <- as.character(trichomedata$strain_num, levels=unique(sladata$strain_num))
trichomedata$species <- as.factor(trichomedata$species)

trichomemodel <- lmer(trichomes_total ~ water_regimen + strain_num + water_regimen:strain_num +(1|sub_block), data=trichomedata)
summary(trichomemodel)
ranef(trichomemodel)
anova(trichomemodel)
plot(trichomemodel)
hist(residuals(trichomemodel))

trichomemodel <- lmer(trichomes_total ~ water_regimen + strain_pres + water_regimen:strain_pres +(1|sub_block), data=trichomedata)
summary(trichomemodel)
ranef(trichomemodel)
anova(trichomemodel)
plot(trichomemodel)
hist(residuals(trichomemodel))

#****Basic GLMMS Harvest****

harvestdata$bench <- as.character(harvestdata$bench)
harvestdata$block <- as.factor(harvestdata$block)
harvestdata$sub_block <- as.factor(harvestdata$sub_block)
harvestdata$tray_number <- as.factor(harvestdata$tray_number)
harvestdata$tray_position <- as.factor(harvestdata$tray_position)
harvestdata$water_regimen <- as.factor(harvestdata$water_regimen)
harvestdata$strain_num <- as.character(harvestdata$strain_num, levels=unique(harvestdata$strain_num))
harvestdata$species <- as.factor(harvestdata$species)

#Shoot biomass model#
shootbiomodel <- lmer(shoot_weight_g ~ water_regimen + strain_num + water_regimen:strain_num +(1|sub_block) + (1|herb_treat), data=harvestdata)
summary(shootbiomodel)
ranef(shootbiomodel)
anova(shootbiomodel)
plot(shootbiomodel)
hist(residuals(shootbiomodel))

#Shoot biomass with herbivory as fixed effect
shootbiomodel <- lmer(shoot_weight_g ~ water_regimen + strain_num + water_regimen:strain_num +(1|sub_block) + (1|herb_treat), data=harvestdata)
summary(shootbiomodel)
ranef(shootbiomodel)
anova(shootbiomodel)
plot(shootbiomodel)
hist(residuals(shootbiomodel))

#Root biomass model#
rootbiomodel <- lmer(root_weight_g ~ water_regimen + strain_num + water_regimen:strain_num +(1|sub_block) +(1|herb_treat), data=harvestdata)
summary(rootbiomodel)
ranef(rootbiomodel)
anova(rootbiomodel)
plot(rootbiomodel)
hist(residuals(rootbiomodel))

#Pod biomass model#
podbiomodel <- lmer(pod_weight_g ~ water_regimen + strain_num + water_regimen:strain_num +(1|sub_block), data=harvestdata)
summary(podbiomodel)
ranef(podbiomodel)
anova(podbiomodel)
plot(podbiomodel)
hist(residuals(podbiomodel))

#Root-shoot ratio model#
rootshootmodel <- lmer(root_shoot_ratio ~ water_regimen + strain_num + water_regimen:strain_num +(1|sub_block), data=harvestdata)
summary(rootshootmodel)
ranef(rootshootmodel)
anova(rootshootmodel)
plot(rootshootmodel)
hist(residuals(rootshootmodel))

#Total biomass model#
totbiomodel <- lmer(tot_biomass_g ~ water_regimen + strain_num + water_regimen:strain_num +(1|sub_block), data=harvestdata)
summary(totbiomodel)
ranef(totbiomodel)
anova(totbiomodel)
plot(totbiomodel)
hist(residuals(totbiomodel))

#Nodule number model#
nodnummodel <- lmer(nodule_count ~ water_regimen + strain_num + water_regimen:strain_num +(1|sub_block), data=harvestdata)
summary(nodnummodel)
ranef(nodnummodel)
anova(nodnummodel)
plot(nodnummodel)
hist(residuals(nodnummodel))


#ggplot2 colorpallette
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Basic ggplot- caterpillar relative growth rate:watering regimen

weightdata.summary <-weightdata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(RGR=mean(relative_growth_rate),se=sd(relative_growth_rate)/sqrt(n()))
weightdata.summary


caterpillargrowthrate <- ggplot(weightdata.summary, aes(x=water_regimen, y=RGR, ymin = RGR-se, ymax = RGR+se, color=water_regimen))


caterpillargrowthrate + geom_pointrange() + geom_errorbar(width = 0.5) + geom_point(size = 5)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+
  ylab("% Weight Gain")+ xlab("Water Regimen")+scale_color_manual('Water Regimen', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+ ylim(0.80, 0.95)+ theme(panel.grid = element_blank(), axis.text.y = element_text(size=rel(2)),
                                                                                                                                                                             axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
                                                                                                                                                                             axis.title = element_text(size=rel(2)),
                                                                                                                                                                             strip.text = element_text(size=rel(2),face="bold"))+ xlab("Strain ID")+scale_color_manual('% Weight Gain', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+
  ylab("% Weight Gain")+ xlab("Water Regimen")+scale_color_manual('Strain ID', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))+ggtitle("No Significant Differences in % Weight Gain by Water Regimen or Strain Identity")+theme(plot.title.position = 'plot', plot.title = element_text(size=18, hjust = 1, face="bold"))


#Basic ggplot- caterpillar relative growth rate: strain #

weightdata.summary <- weightdata %>%
  filter(water_regimen!='drought') %>% 
  group_by(strain_num) %>%
  dplyr::summarise(
    RGR = mean(relative_growth_rate), se =sd(relative_growth_rate)/sqrt(n()), na.rm = TRUE) %>% 
  ungroup()

weightdata.summary

ggplot(weightdata.summary, aes(x=strain_num, y=RGR, ymin = RGR-se, ymax = RGR+se)) +
  geom_pointrange() + geom_errorbar(width = 0.8) + geom_point(size = 5) +
  # theme(panel.grid = element_blank(),
        # axis.text.y = element_text(size=rel(3)),
        # axis.text.x = element_text(size=rel(3),angle=0,hjust=0.5,vjust=1),
        # axis.title = element_text(size=rel(3)),
        # strip.text = element_text(size=rel(3),face="bold"))+
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                              "11", "12", "13", "14", "15", "16", "17", "18", "19", 
                              "20", "21", "22", "23", "24"))+
  ylab("Relative Growth Rate") + xlab("Strain ID")
#export at 1000 x 400

herbivorydata.summary <- herbivorydata %>%
  filter(water_regimen!='drought') %>% 
  group_by(strain_num) %>%
  dplyr::summarise(
    Consumed_Leaf_Area =mean(consumed_leaf_area), se =sd(consumed_leaf_area)/sqrt(n()), na.rm=TRUE) %>% 
  ungroup()

ggplot(herbivorydata.summary, aes(x=strain_num, y=Consumed_Leaf_Area, 
                                  ymin = Consumed_Leaf_Area-se, ymax = Consumed_Leaf_Area+se)) +
  geom_pointrange() + geom_errorbar(width = 0.8) + geom_point(size = 5) +
  # theme(panel.grid = element_blank(),
  # axis.text.y = element_text(size=rel(3)),
  # axis.text.x = element_text(size=rel(3),angle=0,hjust=0.5,vjust=1),
  # axis.title = element_text(size=rel(3)),
  # strip.text = element_text(size=rel(3),face="bold"))+
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                              "11", "12", "13", "14", "15", "16", "17", "18", "19", 
                              "20", "21", "22", "23", "24"))+
  ylab("Percent Leaf Damage") + xlab("Strain ID")
#export at 1000x400

herbivorydata.summary2 <- herbivorydata %>%
  # filter(water_regimen!='drought') %>%
  group_by(strain_num, water_regimen) %>%
  dplyr::summarise(
    Consumed_Leaf_Area =mean(consumed_leaf_area), se =sd(consumed_leaf_area)/sqrt(n()), na.rm=TRUE) %>% 
  ungroup()

ggplot(herbivorydata.summary2, aes(x=strain_num, y=Consumed_Leaf_Area, 
                                  ymin = Consumed_Leaf_Area-se, ymax = Consumed_Leaf_Area+se,
                                  colour=water_regimen)) +
  geom_pointrange(position=position_dodge(width=0.5)) + 
  geom_errorbar(width = 0.8, position=position_dodge(width=0.5)) + 
  geom_point(size = 5, position=position_dodge(width=0.5)) +
  # theme(panel.grid = element_blank(),
  # axis.text.y = element_text(size=rel(3)),
  # axis.text.x = element_text(size=rel(3),angle=0,hjust=0.5,vjust=1),
  # axis.title = element_text(size=rel(3)),
  # strip.text = element_text(size=rel(3),face="bold"))+
  scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                              "11", "12", "13", "14", "15", "16", "17", "18", "19", 
                              "20", "21", "22", "23", "24")) +
  scale_color_manual('Water Regimen', values = c("#FF6E00", "#0091FF")) +
  ylab("Percent Leaf Damage") + xlab("Strain ID") +
  theme(legend.position='bottom')
#export at 1000x400

#Basic ggplot- caterpillar relative growth rate: strain#- grouped by watering regimen

weightdata.summary <- weightdata %>%
  group_by(strain_num, water_regimen) %>%
  dplyr::summarise(
    RGR = mean(relative_growth_rate), se =sd(relative_growth_rate)/sqrt(n()), na.rm = TRUE)

weightdata.summary

caterpillargrowthrate <- ggplot(weightdata.summary, aes(x=strain_num, y=RGR, colour=water_regimen, ymin = RGR-se, ymax = RGR+se))
caterpillargrowthrate + geom_pointrange(position=position_dodge(width=0.5), width = 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5) , size = 10) +scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24"))+
  ylab("% Weight Gain")+scale_color_manual('Water Regimen', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="right")+ theme(panel.grid = element_blank(), axis.text.y = element_text(size=rel(3)),
                                                                                                                                     axis.text.x = element_text(size=rel(3),angle=0,hjust=0.5,vjust=1),
                                                                                                                                     axis.title = element_text(size=rel(3)),
                                                                                                                                     strip.text = element_text(size=rel(3),face="bold"))+scale_color_manual('% Weight Gain', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(3)),
        axis.text.x = element_text(size=rel(3),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(3)),
        strip.text = element_text(size=rel(3),face="bold"))+
  +scale_color_manual('Strain ID', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))



#Basic ggplot- caterpillar relative growth rate: strain#- grouped by watering regimen only significant strains
weightdata.summary <- weightdata %>%
  group_by(strain_num, water_regimen) %>%
  dplyr::summarise(
    RGR = mean(relative_growth_rate), se =sd(relative_growth_rate)/sqrt(n()), na.rm = TRUE)

weightdata.summary

caterpillargrowthrate <- ggplot(weightdata.summary, aes(x=strain_num, y=RGR, colour=water_regimen, ymin = RGR-se, ymax = RGR+se))
caterpillargrowthrate + geom_pointrange() + geom_errorbar(width = 0.8) + geom_point(size = 5)+scale_x_discrete(limits = c("1", "5", "10", "11", "13", "16",  "19", "21"))+scale_color_manual('Water Regimen', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="right")+ ylab("% Weight Gain")+xlab("Strain ID")+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))

caterpillargrowthrate_sep <- ggplot(weightdata.summary, aes(x=water_regimen, y=RGR, colour=water_regimen, ymin = RGR-se, ymax = RGR+se))+scale_color_manual('Water Regimen', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")


#Basic ggplot- caterpillar herbivory (Consumed Leaf Area)

herbivorydata.summary <- herbivorydata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(
    Consumed_Leaf_Area =mean(consumed_leaf_area), se =sd(consumed_leaf_area)/sqrt(n()), na.rm=TRUE
  )

herbivorydata.summary

consumed_leaf_area <- ggplot(herbivorydata.summary, aes(x=water_regimen, y=Consumed_Leaf_Area, color = water_regimen, ymin = Consumed_Leaf_Area-se, ymax = Consumed_Leaf_Area+se))
consumed_leaf_area + geom_pointrange() + geom_errorbar(width = 0.2) + geom_point(size = 5)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(3)),
        axis.text.x = element_text(size=rel(3),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(3)),
        strip.text = element_text(size=rel(3),face="bold"))+
  ylab("Leaf Area Consumed")+ xlab("Water Regimen")+scale_y_continuous(limits = c(2, 3))+scale_color_manual('Water Regimen', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")


#Basic ggplot- caterpillar herbivory (% Leaf area removed)

herbivorydata.summary <- herbivorydata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(
    Percent_Consumed = mean(percent_consumed), se = sd(percent_consumed)/sqrt(n()), na.rm = TRUE
  )
herbivorydata.summary


percent_consumed <- ggplot(herbivorydata.summary, aes(x=water_regimen, y=Percent_Consumed, color=water_regimen, ymin = Percent_Consumed-se, ymax = Percent_Consumed+se))
percent_consumed + geom_pointrange(position=position_dodge(width=0.5), width = 0.2, alpha=0.5) + geom_point(position=position_dodge(width=0.5) , size = 8)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(3)),
        axis.text.x = element_text(size=rel(3),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(3)),
        strip.text = element_text(size=rel(3),face="bold"))+
  ylab("% Leaf Area Consumed")+ xlab("Water Treatment")+scale_color_manual('Water Treatment', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))+scale_y_continuous(limits = c(0, 10))

percent_consumed <- ggplot(herbivorydata.summary, aes(x=water_regimen, y=Percent_Consumed, color=water_regimen, ymin = Percent_Consumed-se, ymax = Percent_Consumed+se))
percent_consumed + geom_boxplot()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(3)),
        axis.text.x = element_text(size=rel(3),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(3)),
        strip.text = element_text(size=rel(3),face="bold"))+
  ylab("% Leaf Area Consumed")+ xlab("Water Treatment")+scale_color_manual('Water Treatment', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))

percent_consumed<-ggplot(herbivorydata.summary, aes(x = factor(water_regimen), y= Percent_Consumed, fill=water_regimen))
percent_consumed+geom_boxplot(width=0.3)+theme_classic()

percent_consumed <- ggplot(herbivorydata.summary, aes(x=water_regimen, y=Percent_Consumed, color=water_regimen, ymin = Percent_Consumed-se, ymax = Percent_Consumed+se))
percent_consumed + geom_pointrange(position=position_dodge(width=0.5), width = 1) + geom_point(position=position_dodge(width=0.5) , size = 10)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(3)),
        axis.text.x = element_text(size=rel(3),angle=0,hjust=0.5),
        axis.title = element_text(size=rel(3)),
        strip.text = element_text(size=rel(3),face="bold"))+
  ylab("% Leaf Area Consumed")+ xlab("Water Treatment")+scale_color_manual('Water Treatment', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))


herbivorydata.summary <- herbivorydata %>%
  group_by(water_regimen, strain_num) %>%
  dplyr::summarise(
    Percent_Consumed = mean(percent_consumed), se = sd(percent_consumed)/sqrt(n()), na.rm = TRUE
  )
herbivorydata.summary

percent_consumed <- ggplot(herbivorydata.summary, aes(x=strain_num, y=Percent_Consumed, color=water_regimen, ymin = Percent_Consumed-se, ymax = Percent_Consumed+se))
percent_consumed + geom_pointrange() + geom_errorbar(width = 0.8) + geom_point(size = 5)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24"))+
  ylab("% Leaf Area Consumed")+ xlab("Strain ID")+scale_color_manual('Water Regimen', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+scale_color_manual('% Weight Gain', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+
  ylab("% Leaf Area Consumed")+ xlab("Rhizobia Strain ID")+scale_color_manual('Strain ID', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))+
  theme(text = element_text(family = "Times New Roman"))


#Basic ggplot- Total leaf area




#Basic ggplot- trichome density

trichomedata.summary <- trichomedata %>%
  group_by(water_regimen, strain_pres) %>%
  dplyr::summarise(
    Trichomes_Total = mean(trichomes_total), se = sd(trichomes_total)/sqrt(n())
  )
trichomedata.summary

trichome <- ggplot(trichomedata.summary, aes(x=strain_pres, y=Trichomes_Total, color=water_regimen, ymin = Trichomes_Total-se, ymax = Trichomes_Total+se))



trichome + geom_errorbar(position=position_dodge(width=0.5), width = 0.2) + geom_point(position=position_dodge(width=0.5) , size = 5) +  theme(panel.grid = element_blank(),
                                                                                                                                               axis.text.y = element_text(size=rel(2)),axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1), axis.title = element_text(size=rel(2)), strip.text = element_text(size=rel(2),face="bold"))+
  ylab("Trichome Density")+ xlab("Strain Presence")+
  scale_color_manual('Water Regimen', values = c("#FF6E00", "#0091FF"))+ scale_y_continuous(breaks = seq(200, 500, by = 50))+theme(legend.position="none")+theme(axis.line = element_line(size = 2))+ theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+ theme(legend.text=element_text(size=12))+theme(legend.title = element_text(size=12))

trichome_plot <- trichome + geom_pointrange(position=position_dodge(width=0.5), width = 0.2) + geom_point(position=position_dodge(width=0.5) , size = 6) +  theme(panel.grid = element_blank(),
                                                                                                                                                                  axis.text.y = element_text(size=rel(3)),axis.text.x = element_text(size=rel(3),angle=0,hjust=0.5,vjust=1), axis.title = element_text(size=rel(3)), strip.text = element_text(size=rel(3),face="bold"))+
  ylab("Trichome Density")+ xlab("Rhizobia Presence")+
  scale_color_manual('Water Regimen', values = c("#FF6E00", "#0091FF"))+ scale_y_continuous(breaks = seq(200, 500, by = 50))+theme(legend.position="none")+theme(axis.line = element_line(size = 2))+ theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+ theme(legend.text=element_text(size=12))+theme(legend.title = element_text(size=12))

trichome_plot

#Basic ggplot- specific leaf area grouped by strain ID and watering regimen 

sla.summary <- sladata %>%
  group_by(water_regimen, strain_num) %>%
  dplyr::summarise(
    Specific_leaf_area = mean(specific_leaf_area), se = sd(specific_leaf_area)/sqrt(n())
  )
sla.summary

sla_plot <- ggplot(sla.summary, aes(x=strain_num, y=Specific_leaf_area, color=water_regimen, ymin = Specific_leaf_area-se, ymax = Specific_leaf_area+se))
sla_plot + geom_pointrange() + geom_point(size = 5)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24"))+
  ylab("Specific Leaf Area")+ xlab("Strain ID")+scale_color_manual('Water Regimen', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+scale_color_manual('Specific Leaf Area', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+
  ylab("Specific Leaf Area (mm^2/mg)")+ xlab("Rhizobia Strain ID")+scale_color_manual('Strain ID', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))+
  theme(text = element_text(family = "Times New Roman"))

#Basic ggplot- specific leaf area grouped by watering regimen only

sla.summary <- sladata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(
    Specific_leaf_area = mean(specific_leaf_area), se = sd(specific_leaf_area)/sqrt(n())
  )
sla.summary

sla_pl <- ggplot(sla.summary, aes(x=water_regimen, y=Specific_leaf_area, color=water_regimen, ymin = Specific_leaf_area-se, ymax = Specific_leaf_area+se))
sla_plot <-sla_pl + geom_pointrange()  + geom_point(size = 5)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+scale_color_manual('Water Treatment', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))+theme(legend.position="none")+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+
  ylab("Specific Leaf Area (mm^2/mg)")+ xlab("Water Treatment")+scale_color_manual('water_regimen', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))

sla_plot

#Basic ggplot- LDMC grouped by strain_ID and watering regimen
ldmc.summary <- sladata %>%
  group_by(water_regimen, strain_num) %>%
  dplyr::summarise(
    Leaf_dry_matter_content = mean(leaf_dry_matter_content), se = sd(leaf_dry_matter_content)/sqrt(n())
  )
ldmc.summary

ldmc_plot <- ggplot(ldmc.summary, aes(x=strain_num, y=Leaf_dry_matter_content, color=water_regimen, ymin = Leaf_dry_matter_content-se, ymax = Leaf_dry_matter_content+se))
ldmc_plot + geom_pointrange() + geom_point(size = 5)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24"))+
  ylab("Specific Leaf Area")+ xlab("Strain ID")+scale_color_manual('Water Regimen', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+scale_color_manual('Leaf Dry Matter Content', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+
  ylab("Leaf Dry Matter Content (mg)")+ xlab("Rhizobia Strain ID")+scale_color_manual('Strain ID', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))

# Basic ggplot- LDMC grouped by water regimen
ldmc.summary <- sladata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(
    Leaf_dry_matter_content = mean(leaf_dry_matter_content), se = sd(leaf_dry_matter_content)/sqrt(n())
  )
ldmc.summary

ldmc_plot <- ggplot(ldmc.summary, aes(x=water_regimen, y=Leaf_dry_matter_content, color=water_regimen, ymin = Leaf_dry_matter_content-se, ymax = Leaf_dry_matter_content+se))
ldmc_plot_water <- ldmc_plot + geom_pointrange()  + geom_point(size = 5)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))+theme(legend.position="none")+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+
  ylab("Leaf Dry Matter Content (mg)")+ xlab("Water Treatment")+scale_color_manual('Water Treatment', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))

ldmc_plot_water

#Basic ggplot- PhiNPQ grouped by strain_ID and water_regimen
PhiNPQ.summary <- photosynqdata %>%
  group_by(water_regimen, strain_num) %>%
  dplyr::summarise(
    Phinonphotoquench = mean(PhiNPQ), se = sd(PhiNPQ)/sqrt(n())
  )
PhiNPQ.summary

PhiNPQ_plot <- ggplot(PhiNPQ.summary, aes(x=strain_num, y=Phinonphotoquench, color=water_regimen, ymin = Phinonphotoquench-se, ymax =Phinonphotoquench+se))
PhiNPQ_plot + geom_pointrange() + geom_point(size = 5)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24"))+
  ylab("PhiNPQ")+ xlab("Strain ID")+scale_color_manual('Water Regimen', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+ theme(legend.position="none")+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+
  ylab("PhiNPQ")+ xlab("Rhizobia Strain ID")+scale_color_manual('Strain ID', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))

#Basic ggplot- PhiNPQ grouped by watering regimen
PhiNPQ.summary <- photosynqdata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(
    Phinonphotoquench = mean(PhiNPQ), se = sd(PhiNPQ)/sqrt(n())
  )
PhiNPQ.summary


PhiNPQ_pl <- ggplot(PhiNPQ.summary, aes(x=water_regimen, y=Phinonphotoquench, color=water_regimen, ymin = Phinonphotoquench-se, ymax = Phinonphotoquench+se))
PhiNPQ_plot <- PhiNPQ_pl + geom_pointrange()  + geom_point(size = 5)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+scale_color_manual('Water Treatment', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))+theme(legend.position="none")+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+
  ylab("PhiNPQ")+ xlab("Water Treatment")+scale_color_manual('water_regimen', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))

PhiNPQ_plot

# Basic ggplot- Leaf temperature differential grouped by strain_ID and water_regimen

leaftempdiff.summary <- photosynqdata %>%
  group_by(water_regimen, strain_num) %>%
  dplyr::summarise(
    Leaf_Temperature_Differential = mean(Leaf.Temperature.Differential), se = sd(Leaf.Temperature.Differential)/sqrt(n())
  )
leaftempdiff.summary

leaftempdiff_plot <- ggplot(leaftempdiff.summary, aes(x=strain_num, y=Leaf_Temperature_Differential, color=water_regimen, ymin = Leaf_Temperature_Differential-se, ymax =Leaf_Temperature_Differential+se))
leaftempdiff_plot + geom_pointrange() + geom_point(size = 5)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24"))+
  ylab("Leaf Temperature Differential")+ xlab("Strain ID")+scale_color_manual('Water Regimen', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+ theme(legend.position="none")+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+
  ylab("Leaf Temperature Differential")+ xlab("Rhizobia Strain ID")+scale_color_manual('Strain ID', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))

# Basic ggplot Leaf temperature differential grouped by water regimen
leaftempdiff.summary <- photosynqdata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(
    Leaf_Temperature_Differential = mean(Leaf.Temperature.Differential), se = sd(Leaf.Temperature.Differential)/sqrt(n())
  )
leaftempdiff.summary

leaftempdiff_plot <- ggplot(leaftempdiff.summary, aes(x=water_regimen, y=Leaf_Temperature_Differential, color=water_regimen, ymin = Leaf_Temperature_Differential-se, ymax =Leaf_Temperature_Differential+se))

leaftempdiff_pl <- leaftempdiff_plot + geom_pointrange()  + geom_point(size = 5)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+scale_color_manual('Water Treatment', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))+theme(legend.position="none")+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+
  ylab("Leaf Temperature Differential (°C)")+ xlab("Water Treatment")+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))

leaftempdiff_pl

# Basic ggplot NPQt grouped by water regimen

NPQt.summary <- photosynqdata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(
    NPQt_fin = mean(NPQt), se = sd(NPQt)/sqrt(n())
  )
NPQt.summary

NPQt_pl <- ggplot(NPQt.summary, aes(x=water_regimen, y=NPQt_fin, color=water_regimen, ymin = NPQt_fin-se, ymax =NPQt_fin+se))

NPQt_plot <- NPQt_pl + geom_pointrange()  + geom_point(size = 5)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+scale_color_manual('Water Treatment', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))+theme(legend.position="none")+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+
  ylab("NPQt")+ xlab("Water Treatment")+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))

NPQt_plot

# Basic ggplot PhiNO grouped by watering regimen

PhiNO.summary <- photosynqdata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(
    PhiNO_fin = mean(PhiNO), se = sd(PhiNO)/sqrt(n())
  )
PhiNO.summary

PhiNO_pl <- ggplot(PhiNO.summary, aes(x=water_regimen, y=PhiNO_fin, color=water_regimen, ymin = PhiNO_fin-se, ymax = PhiNO_fin+se))

PhiNO_plot <- PhiNO_pl + geom_pointrange()  + geom_point(size = 5)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+scale_color_manual('Water Treatment', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))+theme(legend.position="none")+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+
  ylab("PhiNO")+ xlab("Water Treatment")+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))

PhiNO_plot

# Basic ggplot Phi2 grouped by watering regimen

Phi2.summary <- photosynqdata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(
    Phi2_fin = mean(Phi2), se = sd(Phi2)/sqrt(n())
  )
Phi2.summary

Phi2_pl <- ggplot(Phi2.summary, aes(x=water_regimen, y=Phi2_fin, color=water_regimen, ymin = Phi2_fin-se, ymax = Phi2_fin+se))

Phi2_plot <- Phi2_pl + geom_pointrange()  + geom_point(size = 5)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+scale_color_manual('Water Treatment', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))+theme(legend.position="none")+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+
  ylab("Phi2")+ xlab("Water Treatment")+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))

Phi2_plot

# Basic ggplot relative chlorophyll grouped by watering treatment 

chlorophyll.summary <- photosynqdata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(
    chlorophyll_fin = mean(Relative.Chlorophyll), se = sd(Relative.Chlorophyll)/sqrt(n())
  )
chlorophyll.summary

chlorophyll_pl <- ggplot(chlorophyll.summary, aes(x=water_regimen, y=chlorophyll_fin, color=water_regimen, ymin = chlorophyll_fin-se, ymax = chlorophyll_fin+se))

chlorophyll_plot <- chlorophyll_pl + geom_pointrange()  + geom_point(size = 5)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+scale_color_manual('Water Treatment', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))+theme(legend.position="none")+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+
  ylab("Relative Chlorophyll")+ xlab("Water Treatment")+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))

chlorophyll_plot

# Basic ggplot leaf thickness grouped by water regimen 

thickness.summary <- photosynqdata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(
    thickness_fin = mean(Thickness), se = sd(Thickness)/sqrt(n())
  )
thickness.summary

thickness_pl <- ggplot(thickness.summary, aes(x=water_regimen, y=thickness_fin, color=water_regimen, ymin = thickness_fin-se, ymax = thickness_fin+se))

thickness_plot <- thickness_pl + geom_pointrange()  + geom_point(size = 5)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+scale_color_manual('Water Treatment', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))+theme(legend.position="none")+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+
  ylab("Leaf Thickness")+ xlab("Water Treatment")+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))

thickness_plot

# Basic ggplot shoot biomass grouped by strain_ID and water regimen

shoot.summary <- harvestdata %>%
  group_by(water_regimen, strain_num) %>%
  dplyr::summarise(
    shootbiomass = mean(shoot_weight_g), se = sd(shoot_weight_g)/sqrt(n())
  )
shoot.summary

shootbio_plot <- ggplot(shoot.summary, aes(x=strain_num, y=shootbiomass, color=water_regimen, ymin = shootbiomass-se, ymax =shootbiomass+se))
shootbio_plot + geom_pointrange() + geom_point(size = 5)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24"))+
  ylab("Shoot Biomass (g)")+ xlab("Strain ID")+scale_color_manual('Water Regimen', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+ theme(legend.position="none")+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+
  ylab("Shoot Biomass (g)")+ xlab("Rhizobia Strain ID")+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))

# Basic ggplot- Shoot biomass grouped by just watering regimen
shoot.summary <- harvestdata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(
    shootbiomass = mean(shoot_weight_g), se = sd(shoot_weight_g)/sqrt(n())
  )
shoot.summary

shootbio_plot <- ggplot(shoot.summary, aes(x=water_regimen, y=shootbiomass, color=water_regimen, ymin = shootbiomass-se, ymax =shootbiomass+se))
shootbio <- shootbio_plot + geom_pointrange()  + geom_point(size = 5)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+scale_color_manual('Water Treatment', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))+theme(legend.position="none")+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+
  ylab("Shoot Biomass (g)")+ xlab("Water Treatment")+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))

shootbio
# Basic ggplot- root biomass grouped by strain_ID and watering regimen

root.summary <- harvestdata %>%
  group_by(water_regimen, strain_num) %>%
  dplyr::summarise(
    rootbiomass = mean(root_weight_g), se = sd(root_weight_g)/sqrt(n())
  )
root.summary

rootbio_plot <- ggplot(root.summary, aes(x=strain_num, y=rootbiomass, color=water_regimen, ymin = rootbiomass-se, ymax =rootbiomass+se))
rootbio_plot + geom_pointrange() + geom_point(size = 5)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24"))+
  ylab("Root Biomass (g)")+ xlab("Strain ID")+scale_color_manual('Water Regimen', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+ theme(legend.position="none")+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+
  ylab("Root Biomass (g)")+ xlab("Rhizobia Strain ID")+scale_color_manual('Strain ID', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))

# Basic ggplot root biomass grouped by watering regimen
root.summary <- harvestdata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(
    rootbiomass = mean(root_weight_g), se = sd(root_weight_g)/sqrt(n())
  )
root.summary

rootbio_plot <- ggplot(root.summary, aes(x=water_regimen, y=rootbiomass, color=water_regimen, ymin = rootbiomass-se, ymax =rootbiomass+se))
rootbio_plot + geom_pointrange()  + geom_point(size = 5)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+scale_color_manual('Water Treatment', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))+theme(legend.position="none")+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+
  ylab("Root Biomass (g)")+ xlab("Water Treatment")+scale_color_manual('water_regimen', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))

#Basic ggplot total biomass grouped by watering regimen
totbio.summary <- harvestdata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(
    tot_biomass = mean(tot_biomass_g), se = sd(tot_biomass_g)/sqrt(n())
  )
totbio.summary

totbio_plot <- ggplot(totbio.summary, aes(x=water_regimen, y=tot_biomass, color=water_regimen, ymin = tot_biomass-se, ymax =tot_biomass+se))
totbio <- totbio_plot + geom_pointrange()  + geom_point(size = 5)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+scale_color_manual('Water Treatment', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))+theme(legend.position="none")+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+
  ylab("Total Biomass (g)")+ xlab("Water Treatment")+scale_color_manual('water_regimen', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))

# Basic ggplot root:shoot ratio grouped by watering regimen

rootshoot.summary <- harvestdata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(
    rootshootratio = mean(root_shoot_ratio), se = sd(root_shoot_ratio)/sqrt(n())
  )
rootshoot.summary

rootshoot_pl <- ggplot(rootshoot.summary, aes(x=water_regimen, y=rootshootratio, color=water_regimen, ymin = rootshootratio-se, ymax =rootshootratio+se))
rootshoot_plot <- rootshoot_pl + geom_pointrange()  + geom_point(size = 5)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+scale_color_manual('Water Treatment', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))+theme(legend.position="none")+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+
  ylab("Root:Shoot Ratio")+ xlab("Water Treatment")+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))+ scale_y_continuous(breaks = seq(1.05, 1.25, by = .05))

rootshoot_plot

# Basic ggplot nodule number grouped by watering regimen

nodule.summary <- harvestdata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(
    nodule_fin = mean(nodule_count), se = sd(nodule_count)/sqrt(n())
  )
nodule.summary

nodule_pl <- ggplot(nodule.summary, aes(x=water_regimen, y=nodule_fin, color=water_regimen, ymin = nodule_fin-se, ymax = nodule_fin+se))

nodule_plot <- nodule_pl + geom_pointrange()  + geom_point(size = 5)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+scale_color_manual('Water Treatment', values = c("#FF6E00", "#0091FF"))+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))+theme(legend.position="none")+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(2)),
        axis.text.x = element_text(size=rel(2),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(2)),
        strip.text = element_text(size=rel(2),face="bold"))+
  ylab("Nodule Number")+ xlab("Water Treatment")+ theme(legend.position="none")+theme(axis.text.x=element_text(colour="black")) + theme(axis.text.y=element_text(colour="black"))+theme(axis.line = element_line(size = 3))

nodule_plot

##############Random Slope Intercept Modeling##################





##### FROM BRENDAN'S CODE - looking at polyculture greenhousue results, only control watering #####
library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(directlabels)
library(vegan)
library(tidyverse)


setwd("C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\bean_dip_2018-2024\\greenhouse trials\\beanDIP 2023 Diversity Greenhouse Experiment\\Datasheets")
getwd()

#Reference dataframe
reference <- read.csv('Reference\\div2023_reference.csv')

#List of predictor variable dataframes
caterpillar_weights <- read.csv('Feeding Trial/div2023_feedingtrial_clean.csv')
herbivory <- read.csv('Feeding Trial/div2023_leafbyte_clean.csv')
punch <- read.csv('Leaf Traits/clean data/div2023_leafpunch_clean.csv')
photosynq <- read.csv('Leaf Traits/clean data/div2023_photosynqaugust_clean.csv')
trichomes <- read.csv('Leaf Traits/clean data/div2023_trichomes_clean.csv')
biomass <- read.csv('Harvest/div2023_harvest_clean.csv')



# Read in dataframes

```{r read in dataframes}

reference <- read.csv('Reference/div2023_reference.csv')

weightdata <- reference %>%
  left_join(read.csv('Feeding Trial/div2023_feedingtrial_clean.csv')) %>%
  filter(culture_id != "control") %>%
  filter(!is.na(relative_growth_rate))

herbivorydata <- reference %>%
  left_join(read.csv('Feeding Trial/div2023_leafbyte_clean.csv')) %>%
  filter(culture_id != "control") %>%
  filter(!is.na(percent_consumed))

punchdata <- reference %>%
  left_join(read.csv('Leaf Traits/clean data/div2023_leafpunch_clean.csv')) %>%
  filter(culture_id != "control") %>%
  filter(!is.na(SLA))

photosynqdata <- reference %>%
  left_join(read.csv('Leaf Traits/clean data/div2023_photosynqaugust_clean.csv')) %>%
  filter(culture_id != "control") %>%
  filter(!is.na(Phi2))

trichomedata <- reference %>%
  left_join(read.csv('Leaf Traits/clean data/div2023_trichomes_clean.csv')) %>%
  filter(culture_id !="control") %>%
  filter(!is.na(trichome_density))

harvestdata <- reference %>%
  left_join(read.csv('Harvest/div2023_harvest_clean.csv')) %>%
  filter(culture_id != "control") %>%
  filter(!is.na(totbio.g))

```

```{r subset dataframes into separate monoculture, biculture, and polyculture dataframes}

```

# Convert fixed effects in dataframes

```{r }
######Strain diversity is coded as a numeric. There are more two categories, and given our prediction that strain diversity would change traits along a numeric gradient, and that the number does matter, we can code as a numeric. 

###Caterpillar weight dataframe
weightdata$bench <- as.factor(weightdata$bench)
weightdata$block <- as.factor(weightdata$block)
weightdata$sub_block <- as.factor(weightdata$sub_block)
weightdata$tray_number <- as.factor(weightdata$tray_number)
weightdata$tray_position <- as.factor(weightdata$tray_position)
weightdata$water_regimen <- as.factor(weightdata$water_regimen)
weightdata$herb_treat <- as.factor(weightdata$herb_treat)
weightdata$strain_treatment <- as.factor(weightdata$strain_treatment)
weightdata$strain_diversity <- as.numeric(weightdata$strain_diversity)
weightdata$performance_treatment <- as.factor(weightdata$performance_treatment)
weightdata$culture_id <- as.factor(weightdata$culture_id)

#Herbivory dataframe
herbivorydata$bench <- as.character(herbivorydata$bench)
herbivorydata$sub_block <- as.factor(herbivorydata$sub_block)
herbivorydata$block <- as.factor(herbivorydata$block)
herbivorydata$tray_number <- as.factor(herbivorydata$tray_number)
herbivorydata$tray_position <- as.factor(herbivorydata$tray_position)
herbivorydata$water_regimen <- as.factor(herbivorydata$water_regimen)
herbivorydata$herb_treat <- as.factor(herbivorydata$herb_treat)
herbivorydata$strain_treatment <- as.factor(herbivorydata$strain_treatment)
herbivorydata$strain_diversity <- as.numeric(herbivorydata$strain_diversity)
herbivorydata$performance_treatment <- as.factor(herbivorydata$performance_treatment)

#Leaf punch dataframe
punchdata$bench <- as.character(punchdata$bench)
punchdata$block <- as.factor(punchdata$block)
punchdata$sub_block <- as.factor(punchdata$sub_block)
punchdata$tray_number <- as.factor(punchdata$tray_number)
punchdata$tray_position <- as.factor(punchdata$tray_position)
punchdata$water_regimen <- as.factor(punchdata$water_regimen)
punchdata$herb_treat <- as.factor(punchdata$herb_treat)
punchdata$strain_treatment <- as.factor(punchdata$strain_treatment)
punchdata$strain_diversity <- as.numeric(punchdata$strain_diversity)
punchdata$performance_treatment <- as.factor(punchdata$performance_treatment)

#Trichome dataframe
trichomedata$bench <- as.character(trichomedata$bench)
trichomedata$block <- as.factor(trichomedata$block)
trichomedata$sub_block <- as.factor(trichomedata$sub_block)
trichomedata$tray_number <- as.factor(trichomedata$tray_number)
trichomedata$tray_position <- as.factor(trichomedata$tray_position)
trichomedata$water_regimen <- as.factor(trichomedata$water_regimen)
trichomedata$herb_treat <- as.factor(trichomedata$herb_treat)
trichomedata$strain_treatment <- as.factor(trichomedata$strain_treatment)
trichomedata$strain_diversity <- as.numeric(trichomedata$strain_diversity)
trichomedata$performance_treatment <- as.factor(trichomedata$performance_treatment)

#Photosynq dataframe
photosynqdata$bench <- as.character(photosynqdata$bench)
photosynqdata$block <- as.factor(photosynqdata$block)
photosynqdata$sub_block <- as.factor(photosynqdata$sub_block)
photosynqdata$tray_number <- as.factor(photosynqdata$tray_number)
photosynqdata$tray_position <- as.factor(photosynqdata$tray_position)
photosynqdata$water_regimen <- as.factor(photosynqdata$water_regimen)
photosynqdata$herb_treat <- as.factor(photosynqdata$herb_treat)
photosynqdata$strain_treatment <- as.factor(photosynqdata$strain_treatment)
photosynqdata$strain_diversity <- as.numeric(photosynqdata$strain_diversity)
photosynqdata$performance_treatment <- as.factor(photosynqdata$performance_treatment)

#Harvest dataframe
harvestdata$bench <- as.character(harvestdata$bench)
harvestdata$block <- as.factor(harvestdata$block)
harvestdata$sub_block <- as.factor(harvestdata$sub_block)
harvestdata$tray_number <- as.factor(harvestdata$tray_number)
harvestdata$tray_position <- as.factor(harvestdata$tray_position)
harvestdata$water_regimen <- as.factor(harvestdata$water_regimen)
harvestdata$herb_treat <- as.factor(harvestdata$herb_treat)
harvestdata$strain_treatment <- as.factor(harvestdata$strain_treatment)
harvestdata$strain_diversity <- as.numeric(harvestdata$strain_diversity)
harvestdata$performance_treatment <- as.factor(harvestdata$performance_treatment)
harvestdata$nodule_count <- as.numeric(harvestdata$nodule_count)

```

# Caterpillar GLMMs

```{r }

###Caterpillar RGR. Significant main effect of watering and interaction. 

rgr_model <- lmer(relative_growth_rate ~ water_regimen + strain_diversity + water_regimen:strain_diversity +(1|block), data=weightdata)
summary(rgr_model)
ranef(rgr_model)
anova(rgr_model)
plot(rgr_model)
hist(residuals(rgr_model))

rgr_model_trt <- lmer(relative_growth_rate ~ water_regimen + strain_treatment + water_regimen:strain_treatment +(1|block) + (1|performance_treatment), data=weightdata)
summary(rgr_model_trt)
ranef(rgr_model_trt)
anova(rgr_model_trt)
plot(rgr_model_trt)
hist(residuals(rgr_model_trt))

#Compare estimated marginal means. When I try it this way, it doesn't work with strain diversity as a numeric. Difference between watering is significant in monocultures, but less in bicultures and polycultures. Estimate marginal means for well-watered and droughted groups for each diversity level. #

weightdata$strain_diversity_factor <- factor(weightdata$strain_diversity, 
                                             levels = c(1, 2, 4), 
                                             ordered = TRUE)

rgr_model_factor <- lmer(relative_growth_rate ~ water_regimen + strain_diversity_factor + water_regimen:strain_diversity_factor +(1|block), data=weightdata)
summary(rgr_model_factor)
ranef(rgr_model_factor)
anova(rgr_model_factor)
plot(rgr_model_factor)
hist(residuals(rgr_model_factor))

emm_rgr <- emmeans(rgr_model_factor, pairwise ~ water_regimen|strain_diversity_factor)
emm_rgr

#Change contrasts to look at differences between diversity levels
emmeans_results <- emmeans(rgr_model_factor, ~ strain_diversity_factor * water_regimen)

# View the emmeans table
emmeans_results

# Run pairwise comparison specifically for watering treatment levels only
pairs(emmeans_results, by = "water_regimen")

#Consumed leaf area (cm^2) Watering significant. NS interaction very similar pattern to the RGR interaction  
cons_herbivory_model <- lmer(consumedarea.cm2 ~ water_regimen + strain_diversity + water_regimen:strain_diversity + (1|block), data=herbivorydata)
summary(cons_herbivory_model)
ranef(cons_herbivory_model)
anova(cons_herbivory_model)
plot(cons_herbivory_model)
hist(residuals(cons_herbivory_model))

#Consume leaf area model with strain diversity as factor to do emmeans analysis
herbivorydata$strain_diversity_factor <- factor(herbivorydata$strain_diversity, 
                                                levels = c(1, 2, 4), 
                                                ordered = TRUE)

herb_model_factor <- lmer(consumedarea.cm2 ~ water_regimen + strain_diversity_factor + water_regimen:strain_diversity_factor +(1|block), data=herbivorydata)
summary(herb_model_factor)
ranef(herb_model_factor)
anova(herb_model_factor)
plot(herb_model_factor)
hist(residuals(herb_model_factor))

emm_herb <- emmeans(herb_model_factor, pairwise ~ water_regimen|strain_diversity_factor)
emm_herb

#Change contrasts to look at differences between diversity levels
emmeans_herbresults <- emmeans(herb_model_factor, ~ strain_diversity_factor * water_regimen)

# View the emmeans table
emmeans_herbresults

# Run pairwise comparison specifically for watering treatment levels only
pairs(emmeans_herbresults, by = "water_regimen")

### Percent leaf area removed. Nothing significant with strain diversity as numeric. 
percent_herbivory_model <- lmer(percent_consumed ~ water_regimen + strain_diversity + water_regimen:strain_diversity  + (1|sub_block), data=herbivorydata)
summary(percent_herbivory_model)
ranef(percent_herbivory_model)
anova(percent_herbivory_model)
plot(percent_herbivory_model)
hist(residuals(percent_herbivory_model))

#Total leaf area (cm^2). Nothing significant  
leafarea_model <- lmer(totalarea.cm2 ~ water_regimen + strain_diversity + water_regimen:strain_diversity + (1|block), data=herbivorydata)
summary(leafarea_model)
ranef(leafarea_model)
anova(leafarea_model)
plot(leafarea_model)
hist(residuals(leafarea_model))

#Nothing significant when comparing means
herbivorydata$strain_diversity <- factor(herbivorydata$strain_diversity, 
                                         levels = c(1, 2, 4), 
                                         labels = c("monoculture", "biculture", "polyculture"), 
                                         ordered = TRUE)

#Nothing significant when comparing emmeans, but a similar trend
emm_herb <- emmeans(percent_herbivory_model, pairwise ~ water_regimen|strain_diversity)
emm_herb

#Nothing related to herbivores is significant when comparing the means and controls are removed.

bartlett.test(relative_growth_rate ~ interaction(water_regimen, strain_diversity), data = weightdata)

```


```{r}
#Caterpillar RGR by strain diversity and watering regimen 

weightdata.summary <- weightdata %>%
  group_by(strain_diversity, water_regimen) %>%
  dplyr::summarise(RGR=mean(relative_growth_rate), se=sd(relative_growth_rate)/sqrt(n()))
weightdata.summary

weightdata.summary$water_regimen <- factor(weightdata.summary$water_regimen, levels = c("well-watered", "drought"))

caterpillargrowthrateint_plot <- ggplot(weightdata.summary, aes(x=as.numeric(strain_diversity), y=RGR, ymin = RGR-se, ymax = RGR+se, color=water_regimen))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) + ylab(bquote('Caterpillar Relative Growth Rate ' (mg/day)))+scale_color_manual('Water Regimen', values = c("#0091FF", "#FF6E00")) +xlab("Strain Diversity") + theme_bw() + theme(panel.grid = element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_y_continuous(limits = c(0.1, 0.5), breaks = seq(0.1, 0.5, 0.1)) + theme(legend.position="none") 

#+ annotate("text", x = 1, y = 0.4, label = "***", size = 6)

#Caterpillar relative growth rate plot by strain treatment and watering 
weightdata.summary <- weightdata %>%
  group_by(strain_treatment, water_regimen) %>%
  dplyr::summarise(RGR=mean(relative_growth_rate), se=sd(relative_growth_rate)/sqrt(n()))
weightdata.summary

caterpillargrowthrate_strain <- ggplot(weightdata.summary, aes(x=strain_treatment, y=RGR, ymin = RGR-se, ymax = RGR+se, color= water_regimen))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) + scale_x_discrete(limits = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "BI1", "BI2", "BI3", "BI4", "BI5", "BI6", "BI7", "BI8", "POLY1", "POLY2", "POLY3", "POLY4", "POLY5", "POLY6", "POLY7", "POLY8")) + ylab(bquote('Caterpillar Relative Growth Rate ' (mg/day))) + theme(legend.position="right") + scale_color_manual('Water Regimen', values = c("#FF6E00", "#0091FF"))

caterpillargrowthrate_strain
caterpillargrowthrate_strain + theme_bw() + theme(panel.grid = element_blank()) +
  scale_y_continuous(limits = c(-0.2, 0.8), breaks = seq(-0.2, 0.8, 0.1)) 


#Caterpillar relative growth rate plot by strain treatment
weightdata.summary <- weightdata %>%
  group_by(strain_treatment) %>%
  dplyr::summarise(RGR=mean(relative_growth_rate), se=sd(relative_growth_rate)/sqrt(n()))
weightdata.summary

caterpillargrowthrate_strain <- ggplot(weightdata.summary, aes(x=strain_treatment, y=RGR, ymin = RGR-se, ymax = RGR+se, ))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) + scale_x_discrete(limits = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "BI1", "BI2", "BI3", "BI4", "BI5", "BI6", "BI7", "BI8", "POLY1", "POLY2", "POLY3", "POLY4", "POLY5", "POLY6", "POLY7", "POLY8")) + ylab(bquote('Caterpillar Relative Growth Rate ' (mg/day))) + theme(legend.position="right") 

caterpillargrowthrate_strain
caterpillargrowthrate_strain + theme_bw() + theme(panel.grid = element_blank()) +
  scale_y_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.1)) 

#Consumed leaf area by watering regimen standard dot plot
herbivorydata.summary <- herbivorydata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(Consumed_Leaf_Area=mean(consumedarea.cm2), se=sd(consumedarea.cm2)/sqrt(n()))
herbivorydata.summary

herbivorydata.summary$water_regimen <- factor(herbivorydata.summary$water_regimen, levels = c("well-watered", "drought"))

herbivory_plot <- ggplot(herbivorydata.summary, aes(x=as.factor(water_regimen), y=Consumed_Leaf_Area, ymin = Consumed_Leaf_Area-se, ymax = Consumed_Leaf_Area+se, color=water_regimen))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) + ylab(bquote('Consumed Leaf Area ' (cm^2)))+scale_color_manual('Water Regimen', values = c("#0091FF", "#FF6E00")) + theme(legend.position="right")+xlab("Watering Treatment")

herbivory_plot + theme_bw() + theme(panel.grid = element_blank()) +
  scale_y_continuous(limits = c(2, 3.5), breaks = seq(2, 3.5, 0.25)) 

#Consumed leaf area by watering regimen and strain diversity standard dot plot
herbivorydata.summary <- herbivorydata %>%
  group_by(water_regimen, strain_diversity_factor) %>%
  dplyr::summarise(Consumed_Leaf_Area=mean(consumedarea.cm2), se=sd(consumedarea.cm2)/sqrt(n())) %>% 
  ungroup() %>% 
  mutate(strain_diversity_factor=as.numeric(strain_diversity_factor)) %>% 
  mutate(strain_diversity_factor=ifelse(strain_diversity_factor==3, 4, strain_diversity_factor)) 
herbivorydata.summary

herbivorydata.summary$water_regimen <- factor(herbivorydata.summary$water_regimen, levels = c("well-watered", "drought"))

ggplot(herbivorydata.summary, 
                            aes(x=strain_diversity_factor, y=Consumed_Leaf_Area, 
                                ymin = Consumed_Leaf_Area-se, ymax = Consumed_Leaf_Area+se,
                                color=water_regimen))+
  geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + 
  geom_point(position=position_dodge(width=0.5), size=5) + 
  ylab(bquote('Consumed Leaf Area ' (cm^2)))+
  scale_color_manual('Water Regimen', values = c("#0091FF", "#FF6E00")) +
  xlab("Strain Diversity") +
  scale_y_continuous(limits = c(2, 4), breaks = seq(2, 4, 0.25)) + 
  theme(legend.position="none") 

herbivoryint_plot
```


```{r}

####****Figure 1****###
####*
####****Create textbox model results****###

textbox_caterpillargrowth <- ggdraw() + 
  draw_label("
Strain diversity
Drought**
Strain diversity × Drought*
", fontface ="plain",
             x = -0.5, y = 0.8, hjust = 0.5, vjust = 0.5, size = 10)

textbox_herbivory <- ggdraw() + 
  draw_label("
Strain diversity
Drought*
Strain diversity × Drought
", fontface ="plain",
             x = -0.5, y = 0.8, hjust = 0.5, vjust = 0.5, size = 10)


#****Arrange and save caterpillar RGR and herbivory interaction plot (Figure 1)****###
arrange_cat <- ggarrange(caterpillargrowthrateint_plot, textbox_caterpillargrowth, herbivoryint_plot, textbox_herbivory, ncol = 4, nrow = 1, align = "h", labels = c("A", "", "B"), widths = c(10, 2, 10, 2))

arrange_cat

```

```{r Caterpillar RGR CV and emmeans plots}

#Subset data for 1 and 4 diversity levels 
weightdata.sub <- reference %>%
  left_join(read.csv('Feeding Trial/div2023_feedingtrial_clean.csv')) %>%
  filter(culture_id != "control") %>%
  filter(strain_diversity !="2")%>%
  filter(!is.na(relative_growth_rate))

weightdatacv.summary <- weightdata.sub %>%
  group_by(block, strain_diversity) %>%
  dplyr::summarise(RGR=mean(relative_growth_rate), se=sd(relative_growth_rate)/sqrt(n()), cv = sd(relative_growth_rate)/mean(relative_growth_rate))

weightdatacv.summary2 <- subset(weightdatacv.summary,select=c(1,2,5))

wide_data <- weightdatacv.summary2 %>%
  pivot_wider(names_from = strain_diversity, values_from = cv, names_prefix = "cv_")


# Perform paired t-test. p = 0.40. Nothing significant
t.test(wide_data$cv_1, wide_data$cv_4, paired = TRUE)

#Reaction norm plot plotting coefficients of variation in caterpillar RGR by 1 and 4 stain diversity levels for each block.
cvrgr_plot <- ggplot(weightdatacv.summary2, aes(x = factor(strain_diversity, level=c('1', '4')), y = cv, group = block)) +
  geom_line() + theme_bw() +
  labs(x = "Strain Diversity", y = "Coefficients of variation (CV)") + geom_line(size=1) + geom_point(size =1) 

cvrgr_plot


```

```{r plant trait GLMMs}

##Diversity effects not a large driver of univariate plant traits.

#****SLA**** Nothing significant 0.15 watering:herbivory effect. 
sla_model <- lmer(SLA ~ water_regimen + strain_diversity + herb_treat + water_regimen:strain_diversity  +  herb_treat:strain_diversity + water_regimen:herb_treat + (1|block) , data=punchdata)
summary(sla_model)
ranef(sla_model)
drop1(sla_model, test="Chisq")
anova(sla_model)
plot(sla_model)
hist(residuals(sla_model))

#Rerun model with strain diversity as a factor
punchdata$strain_diversity_factor <- factor(punchdata$strain_diversity, 
                                            levels = c(1, 2, 4), 
                                            ordered = TRUE)

#****PWC**** Nothing significant.  
pwc_model <- lmer(PWC ~ water_regimen + strain_diversity + herb_treat + water_regimen:herb_treat + water_regimen:strain_diversity + strain_diversity:herb_treat + water_regimen:herb_treat + (1|block), data = punchdata)
summary(pwc_model)
ranef(pwc_model)
anova(pwc_model)
plot(pwc_model)
hist(residuals(pwc_model))

#Rerun model with strain diversity as a factor
punchdata$strain_diversity_factor <- factor(punchdata$strain_diversity, 
                                            levels = c(1, 2, 4), 
                                            ordered = TRUE)

lwc_model_factor <- lmer(LWC ~ water_regimen + strain_diversity_factor + herb_treat + water_regimen:herb_treat + water_regimen:strain_diversity_factor + strain_diversity_factor:herb_treat + (1|block), data = punchdata)
summary(lwc_model_factor)
ranef(lwc_model_factor)
anova(lwc_model_factor)
plot(lwc_model_factor)
hist(residuals(lwc_model_factor))

#Nothing significant
emm_lwc <- emmeans(lwc_model_factor, pairwise ~ water_regimen|strain_diversity_factor)
emm_lwc

#****Relative Chlorophyll**** Nothing significant. Watering = 0.15
chlorophyll_model <- lmer(relative_chlorophyll ~ water_regimen + strain_diversity + herb_treat + water_regimen:strain_diversity + water_regimen:herb_treat + strain_diversity:herb_treat + (1|block), data = photosynqdata)
summary(chlorophyll_model)
ranef(chlorophyll_model)
anova(chlorophyll_model)
plot(chlorophyll_model)
hist(residuals(chlorophyll_model))

#Rerun model with strain diversity as a factor
photosynqdata$strain_diversity_factor <- factor(photosynqdata$strain_diversity, 
                                                levels = c(1, 2, 4), 
                                                ordered = TRUE)

chlorophyll_model_factor <- lmer(relative_chlorophyll ~ water_regimen + strain_diversity_factor + herb_treat + water_regimen:herb_treat + water_regimen:strain_diversity_factor + strain_diversity_factor:herb_treat + (1|block), data = photosynqdata)
summary(chlorophyll_model_factor)
ranef(chlorophyll_model_factor)
anova(chlorophyll_model_factor)
plot(chlorophyll_model_factor)
hist(residuals(chlorophyll_model_factor))

#Nothing significant
emm_chloro <- emmeans(chlorophyll_model_factor, pairwise ~ water_regimen|strain_diversity_factor)
emm_chloro

#****Leaf Temperature Differential**** Watering significant
LTD_model <- lmer(LTD ~ water_regimen + strain_diversity + herb_treat + water_regimen:strain_diversity + water_regimen:herb_treat + strain_diversity:herb_treat + (1|block), data = photosynqdata)
summary(LTD_model)
ranef(LTD_model)
anova(LTD_model)
plot(LTD_model)
hist(residuals(LTD_model))

#Rerun model with strain diversity as factor
LTD_model_factor <- lmer(LTD ~ water_regimen + strain_diversity_factor + herb_treat + water_regimen:herb_treat + water_regimen:strain_diversity_factor + strain_diversity_factor:herb_treat + (1|block), data = photosynqdata)
summary(LTD_model_factor)
ranef(LTD_model_factor)
anova(LTD_model_factor)
plot(LTD_model_factor)
hist(residuals(LTD_model_factor))

# Difference in watering treatments in all three diversity levels significant
emm_LTD <- emmeans(LTD_model_factor, pairwise ~ water_regimen|strain_diversity_factor)
emm_LTD

#****Phi2****# Watering significant.
phi2_model <- lmer(Phi2 ~ water_regimen + strain_diversity + herb_treat + water_regimen:herb_treat + water_regimen:strain_diversity + strain_diversity:herb_treat + (1|block), data = photosynqdata)
summary(phi2_model)
ranef(phi2_model)
anova(phi2_model)
plot(phi2_model)
hist(residuals(phi2_model))

#Rerun model with strain diversity as factor
Phi2_model_factor <- lmer(LTD ~ water_regimen + strain_diversity_factor + herb_treat + water_regimen:herb_treat + water_regimen:strain_diversity_factor + strain_diversity_factor:herb_treat + (1|block), data = photosynqdata)
summary(Phi2_model_factor)
ranef(Phi2_model_factor)
anova(Phi2_model_factor)
plot(Phi2_model_factor)
hist(residuals(Phi2_model_factor))

# Difference in watering treatments in all three diversity levels significant. 
emm_Phi2 <- emmeans(Phi2_model_factor, pairwise ~ water_regimen|strain_diversity_factor)
emm_Phi2

#****PhiNO**** Watering significant
PhiNO_model <- lmer(PhiNO ~ water_regimen + strain_diversity + herb_treat + water_regimen:herb_treat + water_regimen:strain_diversity + strain_diversity:herb_treat + (1|block), data = photosynqdata)
summary(PhiNO_model)
ranef(PhiNO_model)
anova(PhiNO_model)
plot(PhiNO_model)
hist(residuals(PhiNO_model))

#Rerun model with strain diversity as factor
PhiNO_model_factor <- lmer(PhiNO ~ water_regimen + strain_diversity_factor + herb_treat + water_regimen:herb_treat + water_regimen:strain_diversity_factor + strain_diversity_factor:herb_treat + (1|block), data = photosynqdata)
summary(PhiNO_model_factor)
ranef(PhiNO_model_factor)
anova(PhiNO_model_factor)
plot(PhiNO_model_factor)
hist(residuals(PhiNO_model_factor))

# Difference in watering treatments in all three diversity levels significant. 
emm_PhiNO <- emmeans(PhiNO_model_factor, pairwise ~ water_regimen|strain_diversity_factor)
emm_PhiNO

#**** NPQt ****# Watering significant
npqt_model <- lmer(NPQt ~ water_regimen + strain_diversity + herb_treat + performance_treatment + water_regimen:strain_diversity  + strain_diversity:herb_treat + (1|block), data = photosynqdata)
summary(npqt_model)
ranef(npqt_model)
anova(npqt_model)
plot(npqt_model)
hist(residuals(npqt_model))

#Rerun model with strain diversity as factor
npqt_model_factor <- lmer(NPQt ~ water_regimen + strain_diversity_factor + herb_treat + water_regimen:herb_treat + water_regimen:strain_diversity_factor + strain_diversity_factor:herb_treat + (1|block), data = photosynqdata)
summary(npqt_model_factor)
ranef(npqt_model_factor)
anova(npqt_model_factor)
plot(npqt_model_factor)
hist(residuals(npqt_model_factor))

# Difference in watering treatments is significant between monocultures and bicultures but not polycultures
emm_npqt <- emmeans(npqt_model_factor, pairwise ~ water_regimen|strain_diversity_factor)
emm_npqt

#**** Leaf thickness ****# Watering significant
thickness_model <- lmer(Thickness ~ water_regimen + strain_diversity + herb_treat + water_regimen:strain_diversity + water_regimen:herb_treat + strain_diversity:herb_treat + (1|block), data = photosynqdata)
summary(thickness_model)
ranef(thickness_model)
anova(thickness_model)
plot(thickness_model)
hist(residuals(thickness_model))

#Rerun model with strain diversity as factor
thickness_model_factor <- lmer(Thickness ~ water_regimen + strain_diversity_factor + herb_treat + water_regimen:herb_treat + water_regimen:strain_diversity_factor + strain_diversity_factor:herb_treat + (1|block), data = photosynqdata)
summary(thickness_model_factor)
ranef(thickness_model_factor)
anova(thickness_model_factor)
plot(thickness_model_factor)
hist(residuals(thickness_model_factor))

# Difference in watering treatments is significant in all three diversity treatments
emm_thickness <- emmeans(thickness_model_factor, pairwise ~ water_regimen|strain_diversity_factor)
emm_thickness

#**** Plant height ****# Watering significant, weak evidence for main effect of herbivory and strain diversity, weak evidence for strain diversity by herbivory interaction
plantheight_model <- lmer(height.cm ~ water_regimen + strain_diversity + herb_treat + water_regimen:strain_diversity + water_regimen:herb_treat + strain_diversity:herb_treat + (1|block), data = photosynqdata)
summary(plantheight_model)
ranef(plantheight_model)
anova(plantheight_model)
plot(plantheight_model)
hist(residuals(plantheight_model))

#Rerun model with strain diversity as factor
plantheight_model_factor <- lmer(height.cm ~ water_regimen + strain_diversity_factor + herb_treat + water_regimen:herb_treat + water_regimen:strain_diversity_factor + strain_diversity_factor:herb_treat + (1|block), data = photosynqdata)
summary(plantheight_model_factor)
ranef(plantheight_model_factor)
anova(plantheight_model_factor)
plot(plantheight_model_factor)
hist(residuals(plantheight_model_factor))

# Difference in watering treatments is significant in all three diversity treatments. The difference between the watering treatments decreased in the polycultures, but not by much. 
emm_height <- emmeans(plantheight_model_factor, pairwise ~ water_regimen|strain_diversity_factor)
emm_height

```

```{r plant functional trait univariate plots}

#Phi2 watering effect plot
Phi2data.summary <- photosynqdata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(PhiII=mean(Phi2), se=sd(Phi2)/sqrt(n()))
Phi2data.summary

Phi2data.summary$water_regimen <- factor(Phi2data.summary$water_regimen, levels = c("well-watered", "drought"))

phi2water_plot <- ggplot(Phi2data.summary, aes(x=as.factor(water_regimen), y=PhiII, ymin = PhiII-se, ymax = PhiII+se, color=water_regimen))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) + ylab(bquote('ΦII ' (ΔφF/φF[m]))) +scale_color_manual('Water Regimen', values = c("#0091FF", "#FF6E00")) + theme(legend.position="right")+xlab("Watering Treatment") + theme(legend.position="right")+xlab("Watering Treatment")+ theme_bw() + theme(panel.grid = element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(limits = c(0.53, 0.59), breaks = seq(0.53, 0.59, 0.01)) + theme(legend.position="none")

phi2water_plot

# Phi2 Effect size
((0.5789043-0.5449684)/0.5789043)*100


#NPQt watering effect plot

npqtdata.summary <- photosynqdata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(npqt=mean(NPQt), se=sd(NPQt)/sqrt(n()))
Phi2data.summary

npqtdata.summary$water_regimen <- factor(Phi2data.summary$water_regimen, levels = c("well-watered", "drought"))

npqtwater_plot <- ggplot(npqtdata.summary, aes(x=as.factor(water_regimen), y=npqt, ymin = npqt-se, ymax = npqt+se, color=water_regimen))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) + ylab(bquote('NPQt')) +scale_color_manual('Water Regimen', values = c("#0091FF", "#FF6E00")) + theme(legend.position="right")+xlab("Watering Treatment") + theme_bw() + theme(panel.grid = element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(limits = c(0.7, 1.9), breaks = seq(0.7, 1.9, 0.2)) + theme(legend.position="none")

npqtwater_plot

((1.553747-1.029184)/1.553747)*100

#PhiNO watering effect plot

PhiNOdata.summary <- photosynqdata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(Phino=mean(PhiNO), se=sd(PhiNO)/sqrt(n()))
PhiNOdata.summary

PhiNOdata.summary$water_regimen <- factor(PhiNOdata.summary$water_regimen, levels = c("well-watered", "drought"))

phinowater_plot <- ggplot(PhiNOdata.summary, aes(x=as.factor(water_regimen), y=Phino, ymin = Phino-se, ymax = Phino+se, color=water_regimen))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) + ylab(bquote('ΦNO ' (F[s]/F[m]))) +scale_color_manual('Water Regimen', values = c("#0091FF", "#FF6E00")) + theme(legend.position="right")+xlab("Watering Treatment") + theme_bw() + theme(panel.grid = element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(limits = c(0.19, 0.23), breaks = seq(0.19, 0.23, 0.01)) + theme(legend.position="none")

phinowater_plot

((0.22674470-0.2040070)/0.22674470)*100

#Leaf temperature differntial plot watering
LTDdata.summary <- photosynqdata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(ltd=mean(LTD), se=sd(LTD)/sqrt(n()))
PhiNOdata.summary

LTDdata.summary$water_regimen <- factor(LTDdata.summary$water_regimen, levels = c("well-watered", "drought"))

LTDwater_plot <- ggplot(LTDdata.summary, aes(x=as.factor(water_regimen), y=ltd, ymin = ltd-se, ymax = ltd+se, color=water_regimen))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) + ylab(bquote('Leaf temperature differential (C°) ')) +scale_color_manual('Water Regimen', values = c("#0091FF", "#FF6E00")) + theme(legend.position="right")+xlab("Watering Treatment")+ theme_bw() + theme(panel.grid = element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(limits = c(-3.2, -2), breaks = seq(-3.2, -2, 0.2)) + theme(legend.position="none") 

LTDwater_plot

#Leaf thickness
thicknessdata.summary <- photosynqdata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(thickness=mean(Thickness), se=sd(Thickness)/sqrt(n()))
thicknessdata.summary

thicknessdata.summary$water_regimen <- factor(thicknessdata.summary$water_regimen, levels = c("well-watered", "drought"))

thicknesswater_plot <- ggplot(thicknessdata.summary, aes(x=as.factor(water_regimen), y=thickness, ymin = thickness-se, ymax = thickness+se, color=water_regimen))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) + ylab(bquote('Leaf thickness ' (mm))) +scale_color_manual('Water Regimen', values = c("#0091FF", "#FF6E00")) +xlab("Watering Treatment") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(limits = c(0.45, 0.65), breaks = seq(0.45, 0.65, 0.05)) + theme(legend.position="none") 

thicknesswater_plot

#Plant height strain diversity by watering interaction

heightdata.summary <- photosynqdata %>%
  group_by(water_regimen, strain_diversity) %>%
  dplyr::summarise(height=mean(height.cm), se=sd(height.cm)/sqrt(n()))
heightdata.summary

heightdata.summary$water_regimen <- factor(thicknessdata.summary$water_regimen, levels = c("well-watered", "drought"))

thicknessint_plot <- ggplot(thicknessdata.summary, aes(x=as.numeric(water_regimen), y=thickness, ymin = thickness-se, ymax = thickness+se, color=water_regimen))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) + ylab(bquote('Plant height ' (cm))) +scale_color_manual('Water Regimen', values = c("#0091FF", "#FF6E00")) + theme(legend.position="right")+xlab("Watering Treatment")

thicknesswater_plot + theme_bw() + theme(panel.grid = element_blank()) 

#Plant height strain diversity by herbivory interaction 
heightdata.summary <- photosynqdata %>%
  group_by(strain_diversity, herb_treat) %>%
  dplyr::summarise(height=mean(height.cm), se=sd(height.cm)/sqrt(n()))
heightdata.summary

height_divherbplot <- ggplot(heightdata.summary, aes(x=as.numeric(strain_diversity), y=height, ymin = height-se, ymax = height+se, color=herb_treat))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) +ylab("Plant height (cm)")+ theme(legend.position="right")+xlab("Strain Diversity")+scale_color_manual('Herbivory Treatment', values = c("green4", "purple")) + theme_bw() + theme(panel.grid = element_blank())

height_divherbplot

#Plant height herbivory plot
heightdata.summary <- photosynqdata %>%
  group_by(herb_treat) %>%
  dplyr::summarise(height=mean(height.cm), se=sd(height.cm)/sqrt(n()))
heightdata.summary

heightdata.summary$herb_treat <- factor(heightdata.summary$herb_treat, levels = c("noherbivory", "herbivory"))

height_herbplot <- ggplot(heightdata.summary, aes(x=herb_treat, y=height, ymin = height-se, ymax = height+se, color=herb_treat))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) +ylab("Plant height (cm)")+ theme(legend.position="right")+xlab("Herbivroy Treatment")+scale_color_manual('Herbivory Treatment', values = c("purple", "green4")) + theme_bw() + theme(panel.grid = element_blank())

height_herbplot

```
```{r, fig.width=10, fig.height=10, fig.fullwidth=TRUE}

####****Figure 2****###
####*
####****Create textbox model results****###

textbox_LTD <- ggdraw() + 
  draw_label("
Strain diversity
Drought**
Herbivory
Strain diversity × Drought
Strain diversity × Herbivory
Drought × Herbivory

", fontface ="plain",
             x = -0.5, y = 0.4, hjust = 0.5, vjust = 0.5, size = 8)

textbox_thickness <- ggdraw() + 
  draw_label("
Strain diversity
Drought**
Herbivory
Strain diversity × Drought
Strain diversity × Herbivory
Drought × Herbivory
", fontface ="plain",
             x = -0.5, y = 0.4, hjust = 0.5, vjust = 0.5, size = 8)

textbox_phi2 <- ggdraw() + 
  draw_label("
Strain diversity
Drought**
Herbivory
Strain diversity × Drought
Strain diversity × Herbivory
Drought × Herbivory
", fontface ="plain",
             x = -0.5, y = 0.8, hjust = 0.5, vjust = 0.5, size = 8)

textbox_phiNO <- ggdraw() + 
  draw_label("
Strain diversity
Drought**
Herbivory
Strain diversity × Drought
Strain diversity × Herbivory
Drought × Herbivory
", fontface ="plain",
             x = -0.5, y = 0.8, hjust = 0.5, vjust = 0.5, size = 8)

textbox_NPQt <- ggdraw() + 
  draw_label("
Strain diversity
Drought**
Herbivory
Strain diversity × Drought
Strain diversity × Herbivory
Drought × Herbivory
", fontface ="plain",
             x = -0.5, y = 0.5, hjust = 0.5, vjust = 0.5, size = 8)

#****Arrange and save univariate plant trait plots (Figure 2)****###
arrange_plantuni <- ggarrange(LTDwater_plot, textbox_LTD, thicknesswater_plot, textbox_thickness, phi2water_plot, textbox_phi2, phinowater_plot, textbox_phiNO, npqtwater_plot, textbox_NPQt , ncol = 4, nrow = 3, align = "h", labels = c("A", "", "B", "", "C", "", "D", "", "E"), widths = c(10, 2, 10, 2)) 

arrange_plantuni



```


```{r biomass GLMMs}

##****Shoot biomass**** Significant interaction between strain diversity and herbivory
shoot_model <- lmer(shoot.g ~ water_regimen + strain_diversity + herb_treat + water_regimen:strain_diversity + water_regimen:herb_treat + strain_diversity:herb_treat + (1|block), data = harvestdata)
summary(shoot_model)
ranef(shoot_model)
anova(shoot_model)
plot(shoot_model)
hist(residuals(shoot_model))
summary(shoot_model)

#Rerun model with strain diversity as factor
harvestdata$strain_diversity_factor <- factor(harvestdata$strain_diversity, 
                                              levels = c(1, 2, 4), 
                                              ordered = TRUE)

shootbio_model_factor <- lmer(shoot.g ~ water_regimen + strain_diversity_factor + herb_treat + water_regimen:herb_treat + water_regimen:strain_diversity_factor + strain_diversity_factor:herb_treat + (1|block), data = harvestdata)
summary(shootbio_model_factor)
ranef(shootbio_model_factor)
anova(shootbio_model_factor)
plot(shootbio_model_factor)
hist(residuals(shootbio_model_factor))

# Difference in herbivory treatments is significant in the bi and polycultures.  
emm_shoot <- emmeans(shootbio_model_factor, pairwise ~ herb_treat|strain_diversity_factor)
emm_shoot

##****Root biomass**** Strain diversity weakly significant
root_model <- lmer(root.g ~ water_regimen + strain_diversity + herb_treat + water_regimen:strain_diversity + water_regimen:herb_treat + strain_diversity:herb_treat + (1|block), data = harvestdata)
summary(root_model)
ranef(root_model)
anova(root_model)
plot(root_model)
hist(residuals(root_model))
summary(root_model)

root_model_factor <- lmer(root.g ~ water_regimen + strain_diversity_factor + herb_treat + water_regimen:herb_treat + water_regimen:strain_diversity_factor + strain_diversity_factor:herb_treat + (1|block), data = harvestdata)
summary(root_model_factor)
ranef(root_model_factor)
anova(root_model_factor)
plot(root_model_factor)
hist(residuals(root_model_factor))

# Difference in herbivory treatments is significant in the biculture treatments  
emm_root <- emmeans(root_model_factor, pairwise ~ herb_treat|strain_diversity_factor)
emm_root

##****Total biomass**** Weak strain diversity by herbivory interaction. 
totbio_model <- lmer(totbio.g ~ water_regimen + strain_diversity + herb_treat + water_regimen:strain_diversity + water_regimen:herb_treat + strain_diversity:herb_treat + (1|block), data = harvestdata)
summary(totbio_model)
ranef(totbio_model)
anova(totbio_model)
plot(totbio_model)
hist(residuals(totbio_model))

#3 way interaction
totbio_model_three <- lmer(totbio.g ~ water_regimen*strain_diversity*herb_treat + (1|block), data = harvestdata)
summary(totbio_model_three)
ranef(totbio_model_three)
anova(totbio_model_three)
plot(totbio_model_three)
hist(residuals(totbio_model_three))

# Diversity as a factor for emmeans

#Rerun model with strain diversity as factor
harvestdata$strain_diversity_factor <- factor(harvestdata$strain_diversity, 
                                              levels = c(1, 2, 4), 
                                              ordered = TRUE)

totbio_model_factor <- lmer(totbio.g ~ water_regimen + strain_diversity_factor + herb_treat + water_regimen:strain_diversity_factor + water_regimen:herb_treat + strain_diversity_factor:herb_treat + (1|block), data = harvestdata)
summary(totbio_model_factor)
ranef(totbio_model_factor)
anova(totbio_model_factor)
plot(totbio_model_factor)
hist(residuals(totbio_model_factor))

# Difference in herbivory treatments is significant in the bi and polyculture treatments 
emm_totbio <- emmeans(totbio_model_factor, pairwise ~ herb_treat|strain_diversity_factor)
emm_totbio

#Change contrasts to look at differences between diversity levels
emmeans_totbioherb <- emmeans(totbio_model_factor, ~ strain_diversity_factor * herb_treat)
emmeans_totbioherb
# Run pairwise comparison specifically for herbivory treatment levels only. Nothing significant
pairs(emmeans_totbioherb, by = "herb_treat")


##****Root:shoot ratio**** Watering and strain diversity significant. 
rootshoot_model <- lmer(root_shoot_ratio ~ water_regimen + strain_diversity + herb_treat + water_regimen:strain_diversity + water_regimen:herb_treat + strain_diversity:herb_treat + (1|block), data = harvestdata)
summary(rootshoot_model)
ranef(rootshoot_model)
anova(rootshoot_model)
plot(rootshoot_model)
hist(residuals(rootshoot_model))

##****Nodule count**** Water regimen and herbivory interaction weakly significant (both equal 0.05). Didn't get any boundary fit errors with block as a random effect.  
nodule_model <- lmer(nodule_count ~ water_regimen + strain_diversity + herb_treat + water_regimen:strain_diversity + water_regimen:herb_treat + strain_diversity:herb_treat + (1|block), data = harvestdata)
summary(nodule_model)
ranef(nodule_model)
anova(nodule_model)
plot(nodule_model)
hist(residuals(nodule_model))

nodule_model_factor <- lmer(nodule_count ~ water_regimen + strain_diversity_factor + herb_treat + water_regimen:herb_treat + water_regimen:strain_diversity_factor + strain_diversity_factor:herb_treat + (1|block), data = harvestdata)
summary(nodule_model_factor)
ranef(nodule_model_factor)
anova(nodule_model_factor)
plot(nodule_model_factor)
hist(residuals(nodule_model_factor))

# Difference in herbivory treatments is significant in drought but not watered
emm_nodule <- emmeans(nodule_model_factor, pairwise ~ herb_treat|water_regimen)
emm_nodule
```

```{r plant performance and allocation trait plots}

# Difference in watering treatments is significant in all three diversity treatments, but bicultures is weakly significant
emm_shoot <- emmeans(shootbio_model_factor, pairwise ~ water_regimen|strain_diversity_factor)
emm_shoot

emm_shoot <- emmeans(shootbio_model_factor, pairwise ~ herb_treat|strain_diversity_factor)
emm_shoot

emm_tot <- emmeans(totbio_model_factor, pairwise ~ herb_treat|strain_diversity_factor)
emm_tot

#Shoot biomass strain diversity by herbivory plot
shootdata.summary <- harvestdata %>%
  group_by(strain_diversity, herb_treat) %>%
  dplyr::summarise(bio=mean(shoot.g), se=sd(shoot.g)/sqrt(n()))
shootdata.summary

shootbio_divherbplot <- ggplot(shootdata.summary, aes(x=as.numeric(strain_diversity), y=bio, ymin = bio-se, ymax = bio+se, color=herb_treat))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) +ylab("Shoot biomass (g)")+ theme(legend.position="right")+xlab("Strain Diversity")+scale_color_manual('Herbivory Treatment', values = c("green4", "purple"))

shootbio_divherbplot + theme_bw() + theme(panel.grid = element_blank())

#Shoot biomass strain treatment herbivory plot

shootdata.summary <- harvestdata %>%
  group_by(strain_diversity, herb_treat) %>%
  dplyr::summarise(bio=mean(shoot.g), se=sd(shoot.g)/sqrt(n()))
shootdata.summary

shootbio_strain <- ggplot(shootdata.summary, aes(x=as.factor(strain_treatment), y=bio, ymin = bio-se, ymax = bio+se))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) + scale_x_discrete(limits = c("M1", "M2", "M3", "M4", "M5", "M6", "M7", "M8", "BI1", "BI2", "BI3", "BI4", "BI5", "BI6", "BI7", "BI8", "POLY1", "POLY2", "POLY3", "POLY4", "POLY5", "POLY6", "POLY7", "POLY8"))+ylab("Plant biomass")+ theme(legend.position="right")+xlab("Strain Treatment")

shootbio_strain

#Shoot biomass strain diversity by watering plot


#Root biomass strain diversity by herbivory plot
rootdata.summary <- harvestdata %>%
  group_by(strain_diversity, herb_treat) %>%
  dplyr::summarise(rbio=mean(root.g), se=sd(root.g)/sqrt(n()))
rootdata.summary

rootbio_divherbplot <- ggplot(rootdata.summary, aes(x=as.numeric(strain_diversity), y=rbio, ymin = rbio-se, ymax = rbio+se, color=herb_treat))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) +ylab("Root biomass (g)")+ theme(legend.position="right")+xlab("Strain Diversity")+scale_color_manual('Herbivory Treatment', values = c("green4", "purple"))

rootbio_divherbplot + theme_bw() + theme(panel.grid = element_blank())

#Root biomass strain diversity plot
rootdata.summary <- harvestdata %>%
  group_by(strain_diversity) %>%
  dplyr::summarise(rbio=mean(root.g), se=sd(root.g)/sqrt(n()))
rootdata.summary

rootbio_divherbplot <- ggplot(rootdata.summary, aes(x=as.numeric(strain_diversity), y=rbio, ymin = rbio-se, ymax = rbio+se))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) +ylab("Root biomass (g)")+ theme(legend.position="right")+xlab("Strain Diversity")

rootbio_divherbplot + theme_bw() + theme(panel.grid = element_blank())

#Root biomass by watering plot
rootdata.summary <- harvestdata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(rbio=mean(root.g), se=sd(root.g)/sqrt(n()))
rootdata.summary

rootbio_water <- ggplot(rootdata.summary, aes(x=water_regimen, y=rbio, ymin = rbio-se, ymax = rbio+se, color=water_regimen))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) + ylab(bquote('Root biomass ' (g)))+scale_color_manual('Water Regimen', values = c("#FF6E00", "#0091FF")) + theme(legend.position="right")+xlab("Watering treatment")+ scale_x_discrete(limits = c("well-watered", "drought"))

rootbio_water+ theme_bw() + theme(panel.grid = element_blank())

#Total biomass strain diversity by herbivory plot
totbiodata.summary <- harvestdata %>%
  group_by(strain_diversity, herb_treat) %>%
  dplyr::summarise(totbio=mean(totbio.g), se=sd(totbio.g)/sqrt(n()))
totbiodata.summary

totbio_divherbplot <- ggplot(totbiodata.summary, aes(x=as.numeric(strain_diversity), y=totbio, ymin = totbio-se, ymax = totbio+se, color=herb_treat))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) + ylab(bquote('Total plant biomass ' (g)))+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +scale_color_manual('Water Regimen', values = c("green4", "purple")) + theme(legend.position="right")+xlab("Strain Diversity") + theme(legend.position="none") 

totbio_divherbplot


#Total biomass strain diversity by watering plot

totbiodata.summary <- harvestdata %>%
  group_by(strain_diversity, water_regimen) %>%
  dplyr::summarise(totbio=mean(totbio.g), se=sd(totbio.g)/sqrt(n()))
totbiodata.summary

totbioint_plot <- ggplot(totbiodata.summary, aes(x=as.numeric(strain_diversity), y=totbio, ymin = totbio-se, ymax = totbio+se, color=water_regimen))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) + ylab(bquote('Plant biomass ' (g)))+scale_color_manual('Water Regimen', values = c("#FF6E00", "#0091FF")) + theme(legend.position="right")+xlab("Strain Diversity")

totbioint_plot + theme_bw() + theme(panel.grid = element_blank())

#Plant height

heightdata.summary <- photosynqdata %>%
  group_by(strain_diversity,water_regimen) %>%
  dplyr::summarise(height=mean(height.cm), se=sd(height.cm)/sqrt(n()))
heightdata.summary

heightdata.summary$water_regimen <- factor(heightdata.summary$water_regimen, levels = c("well-watered", "drought"))

heightint_plot <- ggplot(heightdata.summary, aes(x=as.numeric(strain_diversity), y=height, ymin = height-se, ymax = height+se, color=water_regimen))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) + ylab(bquote('Plant height ' (cm))) +scale_color_manual('Water Regimen', values = c("#0091FF", "#FF6E00")) + theme(legend.position="right")+xlab("Watering Treatment") + theme_bw() + theme(panel.grid = element_blank()) 

heightint_plot

#Total biomass watering plot
totbiodata.summary <- harvestdata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(totbio=mean(totbio.g), se=sd(totbio.g)/sqrt(n()))
totbiodata.summary

totbio_water <- ggplot(totbiodata.summary, aes(x=water_regimen, y=totbio, ymin = totbio-se, ymax = totbio+se, color=water_regimen))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) + ylab(bquote('Plant biomass ' (g)))+scale_color_manual('Water Regimen', values = c("#FF6E00", "#0091FF")) + theme(legend.position="right")+xlab("Watering treatment")+ scale_x_discrete(limits = c("well-watered", "drought"))

totbio_water+ theme_bw() + theme(panel.grid = element_blank())


#Leaf area strain diversity plot
leafareadata.summary <- herbivorydata %>%
  group_by(strain_diversity) %>%
  dplyr::summarise(leafarea=mean(totalarea.cm2), se=sd(totalarea.cm2)/sqrt(n()))
leafareadata.summary

leafarea_divplot <- ggplot(leafareadata.summary, aes(x=as.factor(strain_diversity), y=leafarea, ymin = leafarea-se, ymax = leafarea+se))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) + scale_x_discrete(limits = c("1", "2", "4"))+ylab("Total leaf area (cm^2)")+ theme(legend.position="right")+xlab("Strain Diversity")

leafarea_divplot + theme_bw() + theme(panel.grid = element_blank())

#Leaf area watering plot
leafareadata.summary <- herbivorydata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(leafarea=mean(totalarea.cm2), se=sd(totalarea.cm2)/sqrt(n()))
leafareadata.summary

leafarea_watplot <- ggplot(leafareadata.summary, aes(x=as.factor(water_regimen), y=leafarea, ymin = leafarea-se, ymax = leafarea+se))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) + scale_x_discrete(limits = c("well-watered", "drought"))+ylab("Total leaf area (cm^2)")+ theme(legend.position="right")+xlab("Water regimen")

leafarea_watplot + theme_bw() + theme(panel.grid = element_blank())

#Root:shoot ratio diversity plot
rootshootdata.summary <- harvestdata %>%
  group_by(strain_diversity) %>%
  dplyr::summarise(rootshoot=mean(root_shoot_ratio), se=sd(root_shoot_ratio)/sqrt(n()))
rootshootdata.summary

rootshoot_divplot <- ggplot(rootshootdata.summary, aes(x=as.factor(strain_diversity), y=rootshoot, ymin = rootshoot-se, ymax = rootshoot+se))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) + scale_x_discrete(limits = c("1", "2", "4"))+ylab("Root:shoot ratio")+ theme(legend.position="right")+xlab("Strain Diversity")

rootshoot_divplot + theme_bw() + theme(panel.grid = element_blank())

#Root:shoot ratio watering plot
rootshootdata.summary <- harvestdata %>%
  group_by(water_regimen) %>%
  dplyr::summarise(rootshoot=mean(root_shoot_ratio), se=sd(root_shoot_ratio)/sqrt(n()))
rootshootdata.summary

rootshootdata.summary$water_regimen <- factor(rootshootdata.summary$water_regimen, levels = c("well-watered", "drought"))

rootshootwater_plot <- ggplot(rootshootdata.summary, aes(x=as.factor(water_regimen), y=rootshoot, ymin = rootshoot-se, ymax = rootshoot+se, color=water_regimen))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) + ylab(bquote('Root:shoot ratio ')) +scale_color_manual('Water Regimen', values = c("#0091FF", "#FF6E00")) + theme(legend.position="right")+xlab("Watering Treatment")

rootshootwater_plot + theme_bw() + theme(panel.grid = element_blank()) 

#Root:shoot strain diversity by watering interaction plot

rootshootdata.summary <- harvestdata %>%
  group_by(strain_diversity, water_regimen) %>%
  dplyr::summarise(rootshoot=mean(root_shoot_ratio), se=sd(root_shoot_ratio)/sqrt(n()))
rootshootdata.summary

rootshootdata.summary$water_regimen <- factor(rootshootdata.summary$water_regimen, levels = c("well-watered", "drought"))

rootshootint_plot <- ggplot(rootshootdata.summary, aes(x=as.numeric(strain_diversity), y=rootshoot, ymin = rootshoot-se, ymax = rootshoot+se, color=water_regimen))+geom_pointrange(position=position_dodge(width=0.5), width= 0.2, alpha=4) + geom_point(position=position_dodge(width=0.5), size=5) + ylab(bquote('Root:shoot ratio')) +scale_color_manual('Water Regimen', values = c("#0091FF", "#FF6E00")) +xlab("Strain Diversity") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(limits = c(0.65, 0.90), breaks = seq(0.65, 0.90, 0.05)) + theme(legend.position="none") 
rootshootint_plot

#Nodule watering x herbivory interaction

harvestnoduledata <- na.omit(harvestdata)

noduledata.summary <- harvestnoduledata %>%
  group_by(water_regimen, herb_treat) %>%
  dplyr::summarise(
    nodules = mean(nodule_count), se = sd(nodule_count)/sqrt(n()), cv = sd(nodule_count)/mean(nodule_count)*100, sd=sd(nodule_count)
  )
noduledata.summary

noduledata.summary$herb_treat <- factor(noduledata.summary$herb_treat, levels = c("noherbivory", "herbivory"))
nodulelabs <- c("non leaf area removed", "leaf area removed")

nodulewaterherb_plot <- ggplot(noduledata.summary, aes(x = herb_treat, y = nodules, group = water_regimen, color=water_regimen, ymax = nodules+se, ymin= nodules-se)) + geom_line(size=1) + geom_point(size=3) + geom_errorbar(size = 1, width = 0.1) + theme_bw() + scale_color_manual('Watering Treatment', values = c("darkorange", "#0091FF")) + xlab('Herbivory Treatment') + ylab(bquote('Nodule count')) + geom_line(size=1) + geom_point(size =1) + theme(panel.grid = element_blank())+
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.y = element_text(color="black", size=rel(1)),
        axis.text.x = element_text(color="black", size=rel(0.8),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(0.8)),
        strip.text = element_text(size=rel(1),face="bold"),
        axis.ticks = element_line(color = "black"),
        legend.position ="none") + scale_x_discrete(labels = nodulelabs, expand = c(0.2, 0))

nodulewaterherb_plot



```
```{r}
####****Figure 3****###
####*
####****Create textbox model results****###

textbox_biomass<- ggdraw() + 
  draw_label("
Strain diversity
Drought
Herbivory
Strain diversity × Drought
Strain diversity × Herbivory*
Drought × Herbivory
", fontface ="plain",
             x = -0.5, y = 0.45, hjust = 0.5, vjust = 0.5, size = 8)

textbox_rootshoot <- ggdraw() + 
  draw_label("
Strain diversity**
Drought**
Herbivory
Strain diversity × Drought
Strain diversity × Herbivory
Drought × Herbivory
", fontface ="plain",
             x = -0.5, y = 0.8, hjust = 0.5, vjust = 0.5, size = 8)

#****Arrange and save biomass and root:shoot ratio plots (Figure 3)****###
arrange_fig3 <- ggarrange(totbio_divherbplot, textbox_biomass, rootshootint_plot, textbox_rootshoot, ncol = 4, nrow = 1, align = "h", labels = c("A", "", "B"), widths = c(10, 2, 10, 2))

arrange_fig3

#totbio_divherbplot, rootshootint_plot
```


```{r multivariate plant trait analysis load in dataframes}

reference <- read.csv('Reference/div2023_reference.csv') #%>%
#filter(culture_id != "control")

#****Primary plant trait dataframe****#
traitdata_all <- reference %>%
  left_join(read.csv('Multivariate Analysis/beanDIP_GH2023_traitall_multivariate.csv')) %>%
  filter(culture_id != "control") %>%
  filter(!is.na(trichome_density)) %>%
  filter(!is.na(SLA)) %>%
  filter(!is.na(Phi2)) %>%
  filter(!is.na(nodule_count))

traitdata_all <- na.exclude(traitdata_all)


```

```{r multivariate plant trait analysis create and standardardize environmental and species matrices}

#Create environmental matrix
traitdata_clean.env <- select(traitdata_all, pot_id, block, sub_block, water_regimen, rhiz_treat, herb_treat, strain_treatment, strain_diversity, performance_treatment, performance_treatment_diversity, trait_selected, culture_id)

#Convert fixed effects to factors and characters for environmental matrix with control plants
traitdata_clean.env$block <- as.factor(traitdata_clean.env$block)
traitdata_clean.env$sub_block <- as.factor(traitdata_clean.env$sub_block)
traitdata_clean.env$water_regimen <- as.factor(traitdata_clean.env$water_regimen)
traitdata_clean.env$herb_treat <- as.factor(traitdata_clean.env$herb_treat)
traitdata_clean.env$strain_treatment <- as.factor(traitdata_clean.env$strain_treatment)
traitdata_clean.env$strain_diversity <- as.numeric(traitdata_clean.env$strain_diversity)
traitdata_clean.env$performance_treatment <- as.factor(traitdata_clean.env$performance_treatment)
traitdata_clean.env$culture_id <- as.factor(traitdata_clean.env$culture_id)

#Create species matrix
traitdata_all_responses <- select(traitdata_all, LTD, NPQt, Phi2, PhiNO, chlorophyll, SLA, LWC, shoot.g, root.g, nodule_count, RS.ratio, height.cm, thickness, trichome_density)

```

```{r multivariate plant trait analysis running the RDA model}

#Standardize variance of the dataframe
###Set RDA scale to 3. Scale = 3 means that that the scores of both species (traits) and and site (plants) are scaled to have unit variance, essentially emphasizing the correlations between the two groups equally. 
scale <- 3

#Standardize variance of the dataframe
traitdata_all_responses <- scale(traitdata_all_responses)

#RDA model- Watering, strain diversity, and herbivory all significant. No interactions significant.   
rda <- rda(traitdata_all_responses ~ water_regimen + strain_diversity + herb_treat + water_regimen:herb_treat + water_regimen:strain_diversity  + Condition(block), data=traitdata_clean.env, scale=TRUE)
rda
anova(rda)
anova(rda, by="term")
anova(rda, by="axis")
anova(rda, by="margin")
summary(rda)

plot(rda)

```

```{r RDA plots}

# Place fig.width = 14, fig.height = 14 in top of code chunk

#Framework of plot using weighted averages of each diversity treatment.
rdadivplot <- ordiplot(rda, display=c("wa"), cex= 0.5, cex.axis=0.8, cex.lab=0.9, tck=0.0, mgp=c(2, 0.5, 0), type="none", scaling=scale, xlim=c(-1, 1), ylim=c(-1, 1), xlab= "RDA 1", ylab = "RDA 2")

#Strain diversity treatments. Seem to be in different spots, polycultures seem to be located in the middle of the trait space, not being pulled in any direction. Qualitatively it looks like the ellipses is a bit larger than the rest, and it looks like monocultures have more constrained responses.  
#with(traitdata_all, ordiellipse(rda, display = "wa", traitdata_all$strain_diversity, kind = "se", conf=0.95, lty = 1, lwd = 2, col = "blue", show.groups ="1"))
#with(traitdata_all, ordiellipse(rda, display = "wa", traitdata_all$strain_diversity, kind = "se", conf=0.95, lty = 1, lwd = 2, col = "red", show.groups ="2"))
#with(traitdata_all, ordiellipse(rda, display = "wa", traitdata_all$strain_diversity, kind = "se", conf=0.95, lty = 1, lwd = 2, col = "green4", show.groups ="4"))

#Watering treatments 
with(traitdata_all, ordiellipse(rda, display = "wa", traitdata_all$water_regimen, kind = "se", conf=0.95, draw="polygon", lty = 1, lwd = 2, col = "#0091FF", show.groups ="well-watered", label=TRUE))
with(traitdata_all, ordiellipse(rda, display = "wa", traitdata_all$water_regimen, kind = "se", conf=0.95, draw= "polygon", lty = 1, lwd = 2, col = "darkorange", show.groups ="drought", label=TRUE))
with(traitdata_all, ordiellipse(rda, display = "wa", traitdata_all$herb_treat, kind = "se", conf=0.95, draw="polygon", lty = 1, lwd = 2, col = "pink", show.groups ="noherbivory"))
with(traitdata_all, ordiellipse(rda, display = "wa", traitdata_all$herb_treat, kind = "se", conf=0.95, draw= "polygon", lty = 1, lwd = 2, col = "green", show.groups ="herbivory"))
plot(envfit(rda~strain_diversity, traitdata_all, display="wa"), arrow.mul=2,add = TRUE, col="brown", cex=1.3)
#legend("topleft", legend=c("well-watered", "drought"), col="black", pt.bg=c("#0091FF", "darkorange"), pch=c(21), pt.cex=1.5)

#Add plant trait response vectors
plot(envfit(rda ~ SLA, traitdata_all, display = "wa"), arrow.mul = 2, add = TRUE, col="black")
plot(envfit(rda ~ nodule_count, traitdata_all, display = "wa"), arrow.mul = 2, add = TRUE, col="black")
plot(envfit(rda ~ shoot.g, traitdata_all, display = "wa"), arrow.mul = 2, add = TRUE, col="black")
plot(envfit(rda ~ root.g, traitdata_all, display = "wa"), arrow.mul = 2, add = TRUE, col="black")
plot(envfit(rda ~ RS.ratio, traitdata_all, display = "wa"), arrow.mul = 2, add = TRUE, col="black")
plot(envfit(rda ~ LTD, traitdata_all, display = "wa"), arrow.mul = 2, add = TRUE, col="black")
plot(envfit(rda ~ NPQt, traitdata_all, display = "wa"), arrow.mul = 2, add = TRUE, col="black")
plot(envfit(rda ~ Phi2, traitdata_all, display = "wa"), arrow.mul = 2, add = TRUE, col="black")
plot(envfit(rda ~ PhiNO, traitdata_all, display = "wa"), arrow.mul = 2, add = TRUE, col="black")
plot(envfit(rda ~ chlorophyll, traitdata_all, display = "wa"), arrow.mul = 2, add = TRUE, col="black")
plot(envfit(rda ~ thickness, traitdata_all, display = "wa"), arrow.mul = 2, add = TRUE, col="black")
plot(envfit(rda ~ height.cm, traitdata_all, display = "wa"), arrow.mul = 2, add = TRUE, col="black")
plot(envfit(rda ~ trichome_density, traitdata_all, display = "wa"), arrow.mul = 2, add = TRUE, col="black")
plot(envfit(rda ~ LWC, traitdata_all, display = "wa"), arrow.mul = 2, add = TRUE, col="black")

#Save RDA figure

```

```{r trichome analysis}

trichomes_div2023_mod <- lmer(trichome_density ~ water_regimen + strain_diversity + herb_treat  + water_regimen:strain_diversity + water_regimen:herb_treat + strain_diversity:herb_treat + (1|sub_block), data=trichomedata)
summary(trichomes_div2023_mod)
drop1(trichomes_div2023_mod, test="Chisq")
anova(trichomes_div2023_mod)


```



```{r trichome plots}

trichome_avg <- trichomedata %>%
  group_by(strain_diversity, water_regimen) %>%
  dplyr::summarise(trichomes = mean(trichome_density), se=sd(trichome_density)/sqrt(n()))

tric <- ggplot(trichome_avg, aes(x=as.numeric(strain_diversity), y=trichomes, ymin= trichomes-se, ymax= trichomes + se, shape=water_regimen)) + geom_point(size=3, position=position_dodge(0.9)) + 
  geom_errorbar(aes(ymin=trichomes-se, ymax=trichomes+se), width=0.2, position=position_dodge(0.9)) +
  ylab('Trichome Density (#/leaf punch) ') + xlab('Strain Diversity') + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) + scale_y_continuous(limits = c(300, 450), breaks = seq(300, 450, 25))

tric



```












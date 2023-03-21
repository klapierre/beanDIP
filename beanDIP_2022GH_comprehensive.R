## beanDIP 2022 Greenhouse
#10/25/2022
#Brendan Randall

library(lme4)
library(ggplot2)
theme_set(
  theme_classic() +
    theme(legend.position = "top")
)
library(RColorBrewer)
library(dplyr)
library(lmerTest)
library (ggpubr)
library(extrafont) 
font_import()
loadfonts(device = "win", y)






# Calling data frames
setwd("C:/Users/brend/Dropbox/PC/Downloads/beanDIP 2022 Greenhouse Datasheets/csv")
getwd()

reference <- read.csv('beanDIP_GH2022_Reference_Analysis.csv')

# Predictor variable data frames 
caterpillar_weights <- read.csv('Feeding Trial_GH2022_caterpillar_weights_clean.csv')
herbivory <- read.csv('Feeding Trial_Herbivory Data_LeafByte_clean.csv')
trichomes <- read.csv('LeafTraits_trichome data_clean.csv')
sla <- read.csv('Leaf Traits_Leaf punch data_raw.csv')
photosynq <- read.csv('Leaf Traits_photosynq_clean.csv')
heightleaves <- read.csv('LeafTraits_heightleaves.csv')
harvest <- read.csv('Harvest Summer 2022 Greenhouse.csv')


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

heightleavesdata <- merge(reference, heightleaves, by ="pot_id")
heightleavesdata <- na.omit(heightleavesdata)

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
  group_by(strain_num) %>%
  dplyr::summarise(
    RGR = mean(relative_growth_rate), se =sd(relative_growth_rate)/sqrt(n()), na.rm = TRUE)
   
weightdata.summary

caterpillargrowthrate <- ggplot(weightdata.summary, aes(x=strain_num, y=RGR, ymin = RGR-se, ymax = RGR+se))
caterpillargrowthrate + geom_pointrange() + geom_errorbar(width = 0.8) + geom_point(size = 5)+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size=rel(3)),
        axis.text.x = element_text(size=rel(3),angle=0,hjust=0.5,vjust=1),
        axis.title = element_text(size=rel(3)),
        strip.text = element_text(size=rel(3),face="bold"))+scale_x_discrete(limits = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24"))+
  ylab("Relative Growth Rate")+ xlab("Strain ID")



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



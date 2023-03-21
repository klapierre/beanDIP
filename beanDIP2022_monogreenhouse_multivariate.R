#beanDIP mponoculture greenhouse 

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
loadfonts(device = "win", y)library(fossil)
library(devtools)
library(vegan)
library(ggvegan)
library(lme4)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggpubr)
library(iNEXT)
library(devtools)
library(readxl)
library(tidyverse)
library(lmerTest)
library(gridExtra)
library(corrplot)
library(scales)


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
traitsall <- read.csv('Leaf_traits_all.csv')


#Merging data frames

photosynqdata <- merge(reference, photosynq, by="pot_id")
photosynqdata <- na.omit(photosynqdata)

harvestdata <- merge(reference, harvest, by="pot_id")

traitdata <- merge(reference, traitsall, by="pot_id")


#traitdata <- merge(traitdata, caterpillar_weights, by="pot_id")

#traitdata <- merge(traitdata, herbivory, by ="pot_id")

traitdata <- na.exclude(traitdata)

#Transformations
#A few ways to go about transformations. Log transformations can only be used on variables that are greater than zero (SLA, LDMC, thickness). Square root transformations can be used for variables which are greater than or equal to zero (nodule count, chlorophyll). Arcsine transformations are used primarily for proportion data (caterpillar growth rates, percent herbivory)

#Residuals of variables we're most interested in look good. Talked with Karin and showed residuals, don't need to do transformations!

#Log transformations
#hist(traitdata$specific_leaf_area)
#traitdata$specific_leaf_area <- log(traitdata$specific_leaf_area)
#hist(traitdata$specific_leaf_area)

#hist(traitdata$leaf_dry_matter_content)
#traitdata$leaf_dry_matter_content <- log(traitdata$leaf_dry_matter_content)
#hist(traitdata$leaf_dry_matter_content)

#hist(traitdata$Thickness)
#traitdata$Thickness <- log(traitdata$Thickness)
#hist(traitdata$Thickness)
# Histogram for thickness still has 1 outlier that is pulling it to the right. Ask Karin if this is ok? 

#hist(traitdata$root_shoot_ratio)
#Residuals skewed right
#traitdata$root_shoot_ratio <- log(traitdata$root_shoot_ratio)
#hist(traitdata$root_shoot_ratio)

#hist(traitdata$trifoliate_leaflet_num)
#traitdata$trifoliate_leaflet_num <- sqrt_trans(traitdata$trifoliate_leaflet_num)
#hist(traitdata$trifoliate_leaflet_num)


# These traits have some skewdness in their histograms pre-transformation. Log transformation to normalize distributions- not necessary. Just my eyes seeing things that aren't there.  

#Square root transformations

#harvestdata$nodule_count <- sqrt(harvestdata$nodule_count) Histogram of the residuals for nodule count seems to be very normal, don't think it needs to be transformed

#hist(traitdata$trifoliate_leaflet_num)
#traitdata$trifoliate_leaflet_num <- sqrt_trans(traitdata$trifoliate_leaflet_num)
#hist(traitdata$trifoliate_leaflet_num)

#traitdata$Relative.Chlorophyll <- sqrt(traitdata$Relative.Chlorophyll) # Histogram of residuals very left skewed, square root transformation. When looking up common transformations for photosynQ data, the only transformations they recommend is square root transforming the PAR data.

#traitdata$PAR <- sqrt(traitdata$PAR)
#hist(traitdata$PAR)

#hist(traitdata$shoot_weight_g)
#Residuals normally distributed

#hist(traitdata$root_weight_g)
#Residuals normally distributed

#hist(traitdata$tot_biomass_g)
#Residuals normally distributed

#hist(traitdata$trifoliate_leaflet_num)
#Residuals skewed
#traitdata$trifoliate_leaflet_num <- log(traitdata$trifoliate_leaflet_num)
#hist(traitdata$trifoliate_leaflet_num)

#Arcsine transformations
#hist(traitdata$relative_growth_rate)

#traitdata$relative_growth_rate <- asin(sqrt(traitdata$relative_growth_rate))# When arc sine transformed, created NaNs. Histogram displayed high normality. 

#traitdata$percent_consumed <- asin(sqrt(traitdata$percent_consumed))
#When arc sine transformed, created NaNs.

#pivot wider 

pivot_trait <- pivot_wider(traitdata, id_cols = pot_id,  values_fill=0)


traitdata_clean.env <- select(traitdata, pot_id, bench, block, sub_block, water_regimen, herb_treat, water_herb, tray_number, tray_position, strain_pres, strain_num, culture_id, selected_from, species)
traitdata_clean.env$pot_id <- as.factor(traitdata_clean.env$pot_id)
traitdata_clean.env$bench <- as.factor(traitdata_clean.env$bench)



traitdata_clean.env$bench <- as.character(traitdata_clean.env$bench)
traitdata_clean.env$block <- as.factor(traitdata_clean.env$block)
traitdata_clean.env$sub_block <- as.factor(traitdata_clean.env$sub_block)
traitdata_clean.env$tray_number <- as.factor(traitdata_clean.env$tray_number)
traitdata_clean.env$tray_position <- as.factor(traitdata_clean.env$tray_position)
traitdata_clean.env$water_regimen <- as.factor(traitdata_clean.env$water_regimen)
traitdata_clean.env$strain_num <- as.factor(traitdata_clean.env$strain_num)
traitdata_clean.env$species <- as.factor(traitdata_clean.env$species)
traitdata_clean.env$herb.treat <- as.factor(traitdata_clean.env$herb_treat)
traitdata_clean.env$water_herb <- as.factor(traitdata_clean.env$water_herb)

#Create trait space for PCA
traitdata_sp <- select(traitdata, Leaf.Angle, Leaf.Temperature.Differential, NPQt, PhiNO, Relative.Chlorophyll, Thickness, specific_leaf_area, leaf_dry_matter_content)

#Running PCA ##
# Not sure if I need to specify scale?
#scl <- 2 ## scaling == 2

pca.trait <- rda(traitdata_sp, scale=TRUE)
summary(pca.trait)

#Extract eigenvalues
ev <- pca.trait$CA$eig

#plot eigenvalues


par(mfrow = c(1,1))
barplot(ev, main= "Eigenvalues", col="bisque", las=2)
abline(h=mean(ev), col="red")


# The plot is not aligning correctly, but we can see from the graph that the first 4 axes reach the average eigenvalue (most important variables). Also don't have all of the variables in the PCA


# PCA with all of the variables

#traitdata_all <- select(traitdata, PAR, Leaf.Angle, Leaf.Temperature.Differential, NPQt, PhiNO, Relative.Chlorophyll, Thickness, specific_leaf_area, leaf_dry_matter_content, nodule_count, shoot_weight_g, root_weight_g, tot_biomass_g, root_shoot_ratio, trifoliate_leaflet_num, relative_growth_rate, total_leaf_area, consumed_leaf_area, percent_consumed)

# without caterpillar data 
traitdata_all <- select(traitdata, Leaf.Angle, Leaf.Temperature.Differential, NPQt, PhiNO, Phi2, Relative.Chlorophyll, Thickness, specific_leaf_area, leaf_dry_matter_content, nodule_count, shoot_weight_g, root_weight_g, tot_biomass_g, root_shoot_ratio)

pca.trait.all <- rda(traitdata_all, scale=TRUE)
summary(pca.trait.all)


#Extract eigenvalues
ev_a <- pca.trait.all$CA$eig

par(mfrow = c(1,1))
barplot(ev_a, main= "Eigenvalues", col="bisque", las=2)
abline(h=mean(ev_a), col="red")

#First two axes seem to be driving variation the most. Also first 6 variables are above the mean (most important variables). PCA doesn't explain which variables are most important, rather which variable specifically is driving variation.

#Also, this data doesn't include anything from the feeding trial, would it be best to have leaf traits and caterpillar data plotted on different multivariate ordinations? 
# Yes, data collected at different levels, can't include herbivore data with the plant traits. 

#Plot PCA. In the plot, each individual replicate is represented as a site object. Need to narrow down to strain as our primary sites. should we be plotting things different. 

biplot(pca.trait.all, scaling=1, type=c("text", "points"))

#Need to get plot to show only strains and species as dots. Should I run the model with only. 


ef <- envfit(pca.trait.all ~ water_regimen, data= traitdata_clean.env)
ef
plot(ef)

############Full RDA################

scale <- 3

#Watering regimen and strain number interaction
#rda <- rda(traitdata_all ~ water_regimen + strain_num  +  water_regimen:strain_num + Condition(sub_block), data= traitdata_clean.env, scale=TRUE)
#rda
#anova(rda, by="term")

#Watering regimen and strain number interaction and watering regimen and herbivory interaction
rda <- rda(traitdata_all ~ water_regimen + strain_num  + herb.treat + water_regimen:herb.treat+ water_regimen:strain_num + Condition(sub_block), data= traitdata_clean.env, scale=TRUE)
rda
anova(rda, by="term")

#Extract % explained by the first two axes--- is this right? 

perc <- round(100*(summary(rda)$cont$importance[2, 1:2]),2)



#Three way interaction
#rda <- rda(traitdata_all ~ water_regimen*herb.treat*strain_num + Condition(sub_block), data= traitdata_clean.env, scale=TRUE)
#rda
#anova(rda, by="term")

#rda <- rda(traitdata_all ~ water_regimen*herb.treat*strain_num + Condition(sub_block), data= traitdata_clean.env, scale=TRUE)
#rda
#anova(rda, by="margin")



rda <- rda(traitdata_all ~ water_regimen + strain_num  + herb.treat +selected_from + water_regimen:herb.treat+ water_regimen:strain_num  + Condition(sub_block), data= traitdata_clean.env, scale=TRUE)
rda
anova(rda, by="term")

# RDA with strain identity pulled out  

#rda <- rda(traitdata_all ~ water_regimen + Condition(sub_block) + Condition(strain_num), data= traitdata_clean.env, scale=TRUE)
#rda

#Scale of 3, run ANOVA. 
scale <- 3
anova(rda)
anova(rda, by="margin")

# Build framework of the plot- plot all strains as ellipses
rdaplot <- ordiplot(rda, display=c("wa"), cex=0.5, cex.axis=0.8, cex.lab=0.9, tck=0.0, mgp=c(1.7, 0.5, 0), type="none", scaling=scale, xlim=c(-1.5,1.5), ylim=c(-2, 2), xlab="RDA 1", ylab="RDA 2")

# Add site information to plot (individual plants listed as numbers)
points(rda, "sites", col="gray", cex= 0.3, pch=17)
text(rda, "sites", col="gray", scaling=scale, cex=0.7, pch=17)
#text(rda, dis="cn", col= "brown")


#Add ellipses/bubbles for each strain
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="azure4", show.groups="1") )

# I think we have the beginning of a plot!! Finally an ellipse is on a plot! Repeat for the other 23 strains

# Line type 1- selected from ambient environments , line type 2 selected from drought environments 
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="aquamarine2", show.groups="2"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="yellow3", show.groups="3"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="darkorange", show.groups="4"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="lightskyblue", show.groups="5"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="magenta", show.groups="6"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="blue", show.groups="7"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="darkorchid4", show.groups="8"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="red", show.groups="9"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="darkseagreen", show.groups="10"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="darkgoldenrod3", show.groups="11"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="coral1", show.groups="12"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="burlywood4", show.groups="13"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="deeppink3", show.groups="14"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="chartreuse3", show.groups="15"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="cadetblue", show.groups="16"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="chocolate4", show.groups="17"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="black", show.groups="18"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="blue4", show.groups="19"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="darkkhaki", show.groups="20"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="darkgreen", show.groups="21"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="cornsilk4", show.groups="22"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="cyan3", show.groups="23"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="darksalmon", show.groups="24"))

#text(rda, "species", col="black", scaling = 3, cex=0.9)

#Lot of colors to choose from. There is definitely a lot going on in this plot but I am really excited looking at it! Interactions playing out on the plot even if it looks crazy. Woohoo the experiment really worked and is telling us cool stuff!! 

# Still, what is going on? 

# add plant traits vectors for variables on top of plot
#plot(envfit(rda ~ PAR, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ Leaf.Temperature.Differential, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ NPQt, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ Phi2, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ PhiNO, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ Relative.Chlorophyll, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ Thickness, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ specific_leaf_area, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ leaf_dry_matter_content, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ nodule_count, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ shoot_weight_g, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ root_weight_g, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ tot_biomass_g, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ root_shoot_ratio, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
#plot(envfit(rda ~ trifoliate_leaflet_num, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
#plot(envfit(rda ~ relative_growth_rate, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
#plot(envfit(rda ~ total_leaf_area, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
#plot(envfit(rda ~ consumed_leaf_area, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
#plot(envfit(rda ~ percent_consumed, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
#plot(envfit(rda ~ water_regimen, traitdata, display="wa"), arrow.mul = 2, add =TRUE, col="red", cex = 0.8)


###########RDA plot- Drought vs well watered ellipses, traits visualized as species. Arrows represent treatments


rdaplot <- ordiplot(rda, display=c("wa"), cex=0.5, cex.axis=0.8, cex.lab=0.9, tck=0.0, mgp=c(1.7, 0.5, 0), type="none", scaling=scale, xlim=c(-2,2), ylim=c(-2, 2), xlab="RDA 1", ylab="RDA 2")

legend("bottomleft", legend = c("watered_herbivory","watered_no_herbivory", "drought_herbivory", "drought_no_herbivory"), pch = 1:4,
       col = c("blue4","cyan3", "darkorange", "yellow3"))

text(rda, "species", col="black", scaling=scale, cex=0.7, pch=17)
#text(rda, dis="cn", col="brown")

# Add drought ellipses- grouped with herbivore treats
with(traitdata, ordiellipse(rda, display="wa", traitdata$water_regimen, kind="se", conf=0.95, lty = 1, lwd =2, col ="darkorange", show.groups="Drought"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$water_regimen, kind="se", conf=0.95, lty = 1, lwd =2, col ="blue", show.groups="Watered"))

# Add herbivory ellipses- grouped with drought treats 

with(traitdata, ordiellipse(rda, display="wa", traitdata$herb_treat, kind="se", conf=0.95, lty= 1, lwd=1, col ="green", show.groups="herbivory"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$herb_treat, kind="se", conf=0.95, lty= 1, lwd=1, col ="purple", show.groups="no_herbivory"))

#Add herbivory_Water treatment ellipses for aggregation plot
with(traitdata, ordiellipse(rda, display="wa", traitdata$water_herb, kind="se", conf=0.95, lty= 1, lwd=1, col ="blue4", show.groups="watered_herbivory"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$water_herb, kind="se", conf=0.95, lty= 2, lwd=2, col ="cyan3", show.groups="watered_no_herbivory"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$water_herb, kind="se", conf=0.95, lty= 1, lwd=1, col ="darkorange", show.groups="drought_herbivory"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$water_herb, kind="se", conf=0.95, lty= 2, lwd=2, col ="yellow3", show.groups="drought_no_herbivory"))

#legend("bottomleft", legend = c("watered_herbivory","watered_no_herbivory", "drought_herbivory", "drought_no_herbivory"), pch = 1:4,
       #col = c("blue4","cyan3", "darkorange", "yellow3"))
# Why are legends so awful to add? 
# Label= TRUE to add labels to plot

#Not using PAR
#plot(envfit(rda ~ PAR, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ Leaf.Temperature.Differential, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ NPQt, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ PhiNO, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ Relative.Chlorophyll, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ Thickness, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ specific_leaf_area, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ leaf_dry_matter_content, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ nodule_count, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ shoot_weight_g, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ root_weight_g, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ tot_biomass_g, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ root_shoot_ratio, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
#plot(envfit(rda ~ trifoliate_leaflet_num, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
#plot(envfit(rda ~ relative_growth_rate, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
#plot(envfit(rda ~ total_leaf_area, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
#plot(envfit(rda ~ consumed_leaf_area, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
#plot(envfit(rda ~ percent_consumed, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )

######Two bubble plots- watered and droughted plants separated to show water_regimen:strain_num significant interaction

#Well-watered plot
rdaplot <- ordiplot(rda, display=c("wa"), cex=0.5, cex.axis=0.8, cex.lab=0.9, tck=0.0, mgp=c(1.7, 0.5, 0), type="none", scaling=scale, xlim=c(-2,2), ylim=c(-2, 2), xlab="RDA 1", ylab="RDA 2")


#Filter doesn't work
#water <- filter(traitdata, water_regimen=='Watered')
#drought <- filter(traitdata, water_regimen=="Drought")
#water



with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="aquamarine2", show.groups="2"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="yellow3", show.groups="3"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="darkorange", show.groups="4"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="lightskyblue", show.groups="5"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="magenta", show.groups="6"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="blue", show.groups="7"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="darkorchid4", show.groups="8"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="red", show.groups="9"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="darkseagreen", show.groups="10"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="darkgoldenrod3", show.groups="11"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="coral1", show.groups="12"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="burlywood4", show.groups="13"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="deeppink3", show.groups="14"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="chartreuse3", show.groups="15"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="cadetblue", show.groups="16"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="chocolate4", show.groups="17"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="black", show.groups="18"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="blue4", show.groups="19"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="darkkhaki", show.groups="20"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="darkgreen", show.groups="21"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="cornsilk4", show.groups="22"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="cyan3", show.groups="23"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="darksalmon", show.groups="24"))

plot(envfit(rda ~ Leaf.Temperature.Differential, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ NPQt, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ PhiNO, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ Relative.Chlorophyll, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ Thickness, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ specific_leaf_area, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ leaf_dry_matter_content, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ nodule_count, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ shoot_weight_g, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ root_weight_g, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ tot_biomass_g, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ root_shoot_ratio, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )


#Drought plot
rdaplot <- ordiplot(rda, display=c("wa"), cex=0.5, cex.axis=0.8, cex.lab=0.9, tck=0.0, mgp=c(1.7, 0.5, 0), type="none", scaling=scale, xlim=c(-2,2), ylim=c(-2, 2), xlab="RDA 1", ylab="RDA 2")


water <- filter(traitdata, water_regimen=="Watered")

with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="aquamarine2", show.groups="2"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="yellow3", show.groups="3"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="darkorange", show.groups="4"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="lightskyblue", show.groups="5"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="magenta", show.groups="6"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="blue", show.groups="7"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="darkorchid4", show.groups="8"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="red", show.groups="9"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="darkseagreen", show.groups="10"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="darkgoldenrod3", show.groups="11"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="coral1", show.groups="12"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="burlywood4", show.groups="13"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="deeppink3", show.groups="14"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="chartreuse3", show.groups="15"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="cadetblue", show.groups="16"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="chocolate4", show.groups="17"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="black", show.groups="18"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="blue4", show.groups="19"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 1, lwd =2, col ="darkkhaki", show.groups="20"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="darkgreen", show.groups="21"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="cornsilk4", show.groups="22"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="cyan3", show.groups="23"))
with(traitdata, ordiellipse(rda, display="wa", traitdata$strain_num, kind="se", conf=0.95, lty = 2, lwd =2, col ="darksalmon", show.groups="24"))

plot(envfit(rda ~ Leaf.Temperature.Differential, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ NPQt, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ PhiNO, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ Relative.Chlorophyll, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ Thickness, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ specific_leaf_area, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ leaf_dry_matter_content, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ nodule_count, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ shoot_weight_g, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ root_weight_g, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ tot_biomass_g, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )
plot(envfit(rda ~ root_shoot_ratio, traitdata, display="wa"), arrow.mul = 2, add = TRUE, col="red", cex = 0.8 )

# If the number of variables increases, corrections for multiple testing becomes desirable
#efit.adj <- efit
#pvals.adj <- p.adjust(efit$vectors$pvals, method = 'bonferroni')
#efit.adj$vectors$pvals <- pvals.adj
#efit.adj


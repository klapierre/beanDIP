################################################################################
##  fieldBacterialCommunityAnalysis.R: Sequence data from UMD field trials
##
##  Author: Kimberly Komatsu
##  Date created: April 28, 2025
################################################################################


#install and load package
# if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
# devtools::install_github("jbisanz/qiime2R")

library(qiime2R)
library(openxlsx)
library(performance)
library(PerformanceAnalytics)
library(nlme)
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

#designating colors for drought and grazing treatments
droughtColor <- c('#6baed6', '#fdbe85', '#fd8d3c', '#e6550d', '#a63603') #from 0 to 99 ##change this to be blue at 0 to red at 99
grazingColor <- c('#ABDEFF', '#469BEC', '#6D882B') #from HHMMM to MMMMM to MLLMM


#set working directory
setwd('C:\\Users\\kjkomatsu\\Box\\RESEARCH - kjkomatsu - Komatsu Lab Research\\beanDIP\\20250708_beanDIP_fieldAmplicon_kim version')


##### Import Data #####
#sequence variants
soils <- read_qza("20250710_beanDIP_fieldSoils\\deblur_output\\deblur_table_final.qza")
nodules <- read_qza("20250710_beanDIP_fieldNodules\\deblur_output\\deblur_table_final.qza")

#metadata
exptMetaData <- read.csv('20250710_beanDIP_fieldNodules\\combined_metadata.csv') %>% 
  mutate(type=ifelse(plate < 12, 'soil', 'nodule'))

#taxonomy
staxonomy <- read_qza("20250710_beanDIP_fieldSoils\\taxa\\classification.qza")
staxonomy<-parse_taxonomy(staxonomy$data)
staxonomy2 <- staxonomy %>% 
  rownames_to_column("Feature.ID")

ntaxonomy <- read_qza("20250710_beanDIP_fieldNodules\\taxa\\classification.qza")
ntaxonomy<-parse_taxonomy(ntaxonomy$data)
ntaxonomy2 <- ntaxonomy %>% 
  rownames_to_column("Feature.ID")

# ##### Shannon's Entropy #####
# bshannon <- read_qza("KomatsuV6V8\\diversity\\shannon_vector.qza")$data %>% 
#   mutate(type='bacteria') %>% 
#   rownames_to_column("SampleID")  #creates sample name column to merge with metadata
# fshannon <- read_qza("KomatsuITS2\\diversity\\shannon_vector.qza")$data %>% 
#   mutate(type='fungi') %>% 
#   rownames_to_column("SampleID")  #creates sample name column to merge with metadata
# 
# shannonTrt <- rbind(bshannon, fshannon) %>%  
#   left_join(metadata) %>% 
#   rename(year=Year) %>% 
#   filter(year!='ctrl')
# 
# # write.csv(shannonTrt, 'GMDR_all_shannons.csv', row.names=F)
# 
# shannonTrt$grazing_category=factor(shannonTrt$grazing_category,levels=c('heavy', 'stable', 'destock'))
# 
# anova_t3(IndVars=c('year', 'rainfall_reduction', 'grazing_treatment'), 
#          DepVar='shannon_entropy', 
#          RndForm='~1|block/paddock/plot', 
#          Data=subset(shannonTrt, site=='FK' & type=='bacteria' & year>2018))
# 
# anova_t3(IndVars=c('year', 'rainfall_reduction', 'grazing_treatment'), 
#          DepVar='shannon_entropy', 
#          RndForm='~1|block/paddock/plot', 
#          Data=subset(shannonTrt, site=='FK' & type=='fungi' & year>2018))
# 
# anova_t3(IndVars=c('year', 'rainfall_reduction', 'grazing_treatment'), 
#          DepVar='shannon_entropy', 
#          RndForm='~1|block/paddock/plot', 
#          Data=subset(shannonTrt, site=='TB' & type=='bacteria' & year>2018))
# 
# anova_t3(IndVars=c('year', 'rainfall_reduction', 'grazing_treatment'), 
#          DepVar='shannon_entropy', 
#          RndForm='~1|block/paddock/plot', 
#          Data=subset(shannonTrt, site=='TB' & type=='fungi' & year>2018))


# ##### Faith's PD #####
# spd <- read_qza("20250710_beanDIP_fieldSoils\\diversity\\faith_pd_vector.qza")$data %>% 
#   mutate(type='soils')
# npd <- read_qza("20250710_beanDIP_fieldNodules\\diversity\\faith_pd_vector.qza")$data %>% 
#   mutate(type='nodules')
# 
# pdTrt <- rbind(bpd, fpd) %>% 
#   rename(SampleID=V1, faith_pd=V2) %>%  #creates sample name column to merge with metadata
#   left_join(metadata) %>% 
#   rename(year=Year) %>% 
#   filter(year!='ctrl') %>% 
#   mutate(site2=ifelse(site=='TB', 'WY', 'MT'))
# 
# # write.csv(pdTrt, 'C:\\Users\\kjkomatsu\\OneDrive - UNCG\\GMDR\\soils ms derived data\\GMDR_all_faith pd.csv', row.names=F)

# pdTrt$grazing_category=factor(pdTrt$grazing_category,levels=c('heavy', 'stable', 'destock'))
# 
# anova_t3(IndVars=c('year', 'rainfall_reduction', 'grazing_treatment'), 
#          DepVar='faith_pd', 
#          RndForm='~1|block/paddock/plot', 
#          Data=subset(pdTrt, site=='FK' & type=='bacteria' & year>2018))
# 
# anova_t3(IndVars=c('year', 'rainfall_reduction', 'grazing_treatment'), 
#          DepVar='faith_pd', 
#          RndForm='~1|block/paddock/plot', 
#          Data=subset(pdTrt, site=='FK' & type=='fungi' & year>2018))
# 
# anova_t3(IndVars=c('year', 'rainfall_reduction', 'grazing_treatment'), 
#          DepVar='faith_pd', 
#          RndForm='~1|block/paddock/plot', 
#          Data=subset(pdTrt, site=='TB' & type=='bacteria' & year>2018))
# 
# anova_t3(IndVars=c('year', 'rainfall_reduction', 'grazing_treatment'), 
#          DepVar='faith_pd', 
#          RndForm='~1|block/paddock/plot', 
#          Data=subset(pdTrt, site=='TB' & type=='fungi' & year>2018))
# 
# pdTrt$grazing_treatment = factor(pdTrt$grazing_treatment, levels=c('destock', 'stable', 'heavy'))
# 
# ggplot(data=barGraphStats(data=subset(pdTrt, type=='fungi' & year>2018), variable="faith_pd", byFactorNames=c("grazing_treatment", "site2", "year")), 
#        aes(x=grazing_treatment, y=mean, color=as.factor(grazing_treatment))) +
#   xlab('Grazing Treatment') + ylab("Fungal Faith's PD") +
#   facet_grid(cols=vars(year), rows=vars(site2), scales='free') +
#   geom_point(size=5) +
#   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=2) +
#   scale_color_manual(values=grazingColor, name='Grazing\nTreatment') +
#   theme(legend.position='none')
# 
# # ggsave('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_GMDR_soils\\figures\\Fig6_FKTB_microbial_FaithPD_grazing.png', width=12, height=6, units='in', dpi=300, bg='white')
# 
# 
# # ggplot(data=barGraphStats(data=subset(pdTrt, site=='TB' & year>2018 & type=='fungi'), variable="faith_pd", byFactorNames=c("grazing_treatment", "site2", "type", "year")), 
# #        aes(x=grazing_treatment, y=mean, color=as.factor(grazing_treatment))) +
# #   xlab('Grazing Treatment') + ylab("Faith's PD") +
# #   facet_grid(cols=vars(year), rows=vars(type), scales='free') +
# #   geom_point(size=5) +
# #   geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2, size=2) +
# #   scale_color_manual(values=grazingColor, name='Grazing\nTreatment') +
# #   theme(legend.position='none')
# # 
# # # ggsave('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_GMDR_soils\\figures\\Fig_TB_fungal_FaithPD_grazing.png', width=12, height=5, units='in', dpi=300, bg='white')




##### Multi-Variate Analysis #####
#Distance matrix
sunifraqDistance <- as.matrix(read_qza('20250710_beanDIP_fieldSoils\\diversity\\weighted_unifrac_distance_matrix.qza')$data)
nunifraqDistance <- as.matrix(read_qza('20250710_beanDIP_fieldNodules\\diversity\\weighted_unifrac_distance_matrix.qza')$data)

sunifracDistanceTrt <- as.data.frame(sunifraqDistance) %>%
  rownames_to_column("SampleID") %>%
  mutate(type='soils') %>%
  left_join(metadata) %>%
  filter(!grepl('ctrl',SampleID))

sunifraqDistanceMatrix <- sunifracDistanceTrt %>%
  select(-type:-livestock_util_2021) %>%
  column_to_rownames(var="SampleID") %>%
  as.matrix()

sunifraqEnv <- sunifracDistanceTrt %>%
  select(SampleID, type:livestock_util_2021)

# #plotting PCOA
unifracPCOA<-read_qza("20250710_beanDIP_fieldSoils\\diversity\\weighted_unifrac_pcoa_results.qza")

unifracPCOATrt<-unifracPCOA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  left_join(pdTrt) %>%
  mutate(site2=ifelse(site=='TB', 'WY',
                      ifelse(site=='FK', 'MT', 'NA')))

# ggplot(subset(unifracPCOATrt, !is.na(site2)), aes(x=PC1, y=PC2, color=as.factor(rainfall_reduction), shape=grazing_treatment)) +
#   facet_wrap(~site2, scales='free') +
#   geom_point(alpha=0.5) + #alpha controls transparency and helps when points are overlapping
#   theme_q2r() +
#   scale_size_continuous(name="Faith's PD") +
#   scale_color_manual(name="Rainfall Reduction", values=droughtColor) +
#   scale_shape_discrete(name="Grazing Treatment")
# 
# # ggsave('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_GMDR_soils\\figures\\Fig_FKTB_bacterial_PCOA.png', width=14, height=7, units='in', dpi=300, bg='white')


##### PERMANOVA #####
bSVrelative <- as.data.frame(apply(bSVs$data, 2, function(x) x/sum(x)*100)) %>% #convert to percent
  rownames_to_column("Feature.ID") %>%
  gather(-Feature.ID, key="SampleID", value="Abundance") %>%
  left_join(btaxonomy2) %>%
  group_by(SampleID, Kingdom, Phylum) %>%
  summarise(Abundance=sum(Abundance)) %>%
  ungroup()  %>%
  filter(!grepl('ctrl',SampleID),
         !(Phylum %in% c('Campilobacterota', 'FCPU426', 'MBNT15'))) %>%
  left_join(metadata, multiple='all') %>%
  select(-type) %>%
  unique() %>%
  mutate(taxa=paste(Kingdom, Phylum, sep='::')) %>%
  select(-Kingdom, -Phylum) %>%
  filter(Year>2018) %>%
  spread(key=taxa, value=Abundance, fill=0)
# write.csv(bSVrelative, 'SVrelative_bacteria_20240717.csv', row.names=F)
chart.Correlation(bSVrelative[,15:41], histogram=T, pch=19)

fSVrelative <- as.data.frame(apply(fSVs$data, 2, function(x) x/sum(x)*100)) %>% #convert to percent
  rownames_to_column("Feature.ID") %>%
  gather(-Feature.ID, key="SampleID", value="Abundance") %>%
  left_join(ftaxonomy2) %>%
  group_by(SampleID, Kingdom, Phylum) %>%
  summarise(Abundance=sum(Abundance)) %>%
  ungroup()  %>%
  filter(!grepl('ctrl',SampleID),
         !(Phylum %in% c('Aphelidiomycota'))) %>%
  left_join(metadata, multiple='all') %>%
  select(-type) %>%
  unique() %>%
  mutate(taxa=paste(Kingdom, Phylum, sep='::')) %>%
  select(-Kingdom, -Phylum) %>%
  filter(Year>2018) %>%
  spread(key=taxa, value=Abundance, fill=0)
# write.csv(fSVrelative, 'SVrelative_fungi_20240717.csv', row.names=F)
chart.Correlation(fSVrelative[,15:27], histogram=T, pch=19)

# #run first 4 lines of each of above taxa code
# SVrelative <- rbind(bSVrelative, fSVrelative) %>%
#   group_by(SampleID, Kingdom, Phylum) %>%
#   summarise(Abundance=sum(Abundance)) %>%
#   ungroup()  %>%
#   filter(!grepl('ctrl',SampleID),
#          !(Phylum %in% c('Campilobacterota', 'FCPU426', 'MBNT15', 'Aphelidiomycota'))) %>%
#   left_join(metadata, multiple='all') %>%
#   select(-type) %>%
#   unique() %>%
#   mutate(taxa=paste(Kingdom, Phylum, sep='::')) %>%
#   select(-Kingdom, -Phylum) %>%
#   filter(Year>2018) %>%
#   spread(key=taxa, value=Abundance, fill=0)
# # write.csv(SVrelative, 'SVrelative_alltaxa_20240717.csv', row.names=F)

set.seed(123)

#FK permanova - bacterial
fkSVmatrix <- bSVrelative %>%
  filter(site=='FK') %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

fkSVenv <- bSVrelative %>%
  filter(site=='FK') %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaFK <- adonis2(formula = fkSVmatrix ~ Year*rainfall_reduction*grazing_treatment,
                       data=fkSVenv,
                       permutations = 999, method = "bray")
print(permanovaFK) #sig rainfall effect (marginal interaction with year)

# ggplot(subset(unifracPCOATrt), aes(x=PC1, y=PC2, color=as.factor(rainfall_reduction), shape=grazing_treatment)) +
#   facet_wrap(~site2, scales='free') +
#   geom_point() + #alpha controls transparency and helps when points are overlapping
#   theme_q2r() +
#   scale_size_continuous(name="Faith's PD") +
#   scale_color_manual(name="Rainfall Reduction", values=droughtColor) +
#   scale_shape_discrete(name="Grazing Treatment")

#FK permanova - fungal
fkSVmatrix <- fSVrelative %>%
  filter(site=='FK') %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

fkSVenv <- fSVrelative %>%
  filter(site=='FK') %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaFK <- adonis2(formula = fkSVmatrix ~ Year*rainfall_reduction*grazing_treatment,
                       data=fkSVenv,
                       permutations = 999, method = "bray")
print(permanovaFK) #sig rainfall*year effect


#FK permanova - fungal 2019
fkSVmatrix <- fSVrelative %>%
  filter(site=='FK', Year==2019) %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

fkSVenv <- fSVrelative %>%
  filter(site=='FK', Year==2019) %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaFK <- adonis2(formula = fkSVmatrix ~ rainfall_reduction,
                       data=fkSVenv,
                       permutations = 999, method = "bray")
print(permanovaFK) #sig rainfall effect

#FK permanova - fungal 2020
fkSVmatrix <- fSVrelative %>%
  filter(site=='FK', Year==2020) %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

fkSVenv <- fSVrelative %>%
  filter(site=='FK', Year==2020) %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaFK <- adonis2(formula = fkSVmatrix ~ rainfall_reduction,
                       data=fkSVenv,
                       permutations = 999, method = "bray")
print(permanovaFK)

#FK permanova - fungal 2021
fkSVmatrix <- fSVrelative %>%
  filter(site=='FK', Year==2021) %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

fkSVenv <- fSVrelative %>%
  filter(site=='FK', Year==2021) %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaFK <- adonis2(formula = fkSVmatrix ~ rainfall_reduction,
                       data=fkSVenv,
                       permutations = 999, method = "bray")
print(permanovaFK) 

#FK permanova - fungal 2022
fkSVmatrix <- fSVrelative %>%
  filter(site=='FK', Year==2022) %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

fkSVenv <- fSVrelative %>%
  filter(site=='FK', Year==2022) %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaFK <- adonis2(formula = fkSVmatrix ~ rainfall_reduction,
                       data=fkSVenv,
                       permutations = 999, method = "bray")
print(permanovaFK) #sig rainfall effect 

ggplot(subset(unifracPCOATrt, type=='fungi' & site=='FK'), aes(x=PC1, y=PC2, color=as.factor(rainfall_reduction))) +
  facet_wrap(~year, scales='free') +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  scale_size_continuous(name="Faith's PD") +
  scale_color_manual(name="Rainfall Reduction", values=droughtColor) +
  scale_shape_discrete(name="Grazing Treatment")



#TB permanova - bacterial
tbSVmatrix <- bSVrelative %>%
  filter(site=='TB') %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

tbSVenv <- bSVrelative %>%
  filter(site=='TB') %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaTB <- adonis2(formula = tbSVmatrix ~ Year*rainfall_reduction*grazing_treatment,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) #significant grazing * year effect

#TB permanova - bacterial 2019
tbSVmatrix <- bSVrelative %>%
  filter(site=='TB', Year==2019) %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

tbSVenv <- bSVrelative %>%
  filter(site=='TB', Year==2019) %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaTB <- adonis2(formula = tbSVmatrix ~ grazing_treatment,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) #significant grazing

#TB permanova - bacterial 2020
tbSVmatrix <- bSVrelative %>%
  filter(site=='TB', Year==2020) %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

tbSVenv <- bSVrelative %>%
  filter(site=='TB', Year==2020) %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaTB <- adonis2(formula = tbSVmatrix ~ grazing_treatment,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) 

#TB permanova - bacterial 2021
tbSVmatrix <- bSVrelative %>%
  filter(site=='TB', Year==2021) %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

tbSVenv <- bSVrelative %>%
  filter(site=='TB', Year==2021) %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaTB <- adonis2(formula = tbSVmatrix ~ grazing_treatment,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) #significant grazing

#TB permanova - bacterial 2022
tbSVmatrix <- bSVrelative %>%
  filter(site=='TB', Year==2022) %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

tbSVenv <- bSVrelative %>%
  filter(site=='TB', Year==2022) %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaTB <- adonis2(formula = tbSVmatrix ~ grazing_treatment,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) 


#TB permanova - fungal
tbSVmatrix <- fSVrelative %>%
  filter(site=='TB') %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

tbSVenv <- fSVrelative %>%
  filter(site=='TB') %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaTB <- adonis2(formula = tbSVmatrix ~ Year*rainfall_reduction*grazing_treatment,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) #significant grazing * year effect


#TB permanova - fungal 2019
tbSVmatrix <- fSVrelative %>%
  filter(site=='TB', Year==2019) %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

tbSVenv <- fSVrelative %>%
  filter(site=='TB', Year==2019) %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaTB <- adonis2(formula = tbSVmatrix ~ rainfall_reduction,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) #sig rainfall effect 

permanovaTB <- adonis2(formula = tbSVmatrix ~ grazing_treatment,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) #no effect

#TB permanova - fungal 2020
tbSVmatrix <- fSVrelative %>%
  filter(site=='TB', Year==2020) %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

tbSVenv <- fSVrelative %>%
  filter(site=='TB', Year==2020) %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaTB <- adonis2(formula = tbSVmatrix ~ rainfall_reduction,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) 

permanovaTB <- adonis2(formula = tbSVmatrix ~ grazing_treatment,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) 

#TB permanova - fungal 2021
tbSVmatrix <- fSVrelative %>%
  filter(site=='TB', Year==2021) %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

tbSVenv <- fSVrelative %>%
  filter(site=='TB', Year==2021) %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaTB <- adonis2(formula = tbSVmatrix ~ rainfall_reduction,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) 

permanovaTB <- adonis2(formula = tbSVmatrix ~ grazing_treatment,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) 

#TB permanova - fungal 2022
tbSVmatrix <- fSVrelative %>%
  filter(site=='TB', Year==2022) %>%
  select(-SampleID:-livestock_util_2021) %>%
  as.matrix()

tbSVenv <- fSVrelative %>%
  filter(site=='TB', Year==2022) %>%
  select(SampleID:livestock_util_2021) %>%
  unique()

permanovaTB <- adonis2(formula = tbSVmatrix ~ rainfall_reduction,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) 

permanovaTB <- adonis2(formula = tbSVmatrix ~ grazing_treatment,
                       data=tbSVenv,
                       permutations = 999, method = "bray")
print(permanovaTB) 

ggplot(subset(unifracPCOATrt, type=='fungi' & site=='TB' & year>2018), aes(x=PC1, y=PC2, color=as.factor(rainfall_reduction))) +
  facet_wrap(~year, scales='free') +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  scale_size_continuous(name="Faith's PD") +
  scale_color_manual(name="Rainfall Reduction", values=droughtColor) +
  scale_shape_discrete(name="Grazing Treatment")

ggplot(subset(unifracPCOATrt, type=='fungi' & site=='TB' & year>2018), aes(x=PC1, y=PC2, color=as.factor(grazing_treatment))) +
  facet_wrap(~year, scales='free') +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  scale_size_continuous(name="Faith's PD") +
  scale_color_manual(name="Grazing Intensity", values=grazingColor) +
  scale_shape_discrete(name="Grazing Treatment")

ggplot(subset(unifracPCOATrt, type=='fungi' & site=='TB' & year>2019), aes(x=PC1, y=PC2, color=as.factor(grazing_treatment))) +
  # facet_wrap(~year, scales='free') +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  scale_size_continuous(name="Faith's PD") +
  scale_color_manual(name="Grazing Intensity", values=grazingColor) +
  scale_shape_discrete(name="Grazing Treatment")

ggplot(subset(unifracPCOATrt, type=='fungi' & site=='TB' & year>2018), aes(x=PC1, y=PC2, color=as.factor(rainfall_reduction))) +
  # facet_wrap(~year, scales='free') +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  scale_size_continuous(name="Faith's PD") +
  scale_color_manual(name="Rainfall Reduction", values=droughtColor) +
  scale_shape_discrete(name="Rainfall Reduction")



### fungal taxa specific responses
fSVsToPlot <- fSVrelative %>% #find the average abundance of a SV
  rename('Fungi::Other'='Fungi::NA') %>% 
  pivot_longer(cols=c('Fungi::Ascomycota':'Fungi::Olpidiomycota'), names_to='taxa', values_to='abundance') %>% 
  group_by(site, plot, Year, taxa) %>% 
  summarize(abundance=sum(abundance)) %>% 
  ungroup() %>% 
  group_by(taxa) %>% 
  summarize(mean_abundance=mean(abundance)) %>% 
  ungroup() %>% 
  arrange(desc(mean_abundance)) %>%
  top_n(20, mean_abundance) %>%
  pull(taxa) #extract only the names from the table

metadata2 <- metadata %>% 
  select(-type) %>% 
  unique()

fSVphylaHeatMap <- fSVrelative %>%
  rename('Fungi::Other'='Fungi::NA') %>%
  pivot_longer(cols=c('Fungi::Ascomycota':'Fungi::Olpidiomycota'), names_to='taxa', values_to='abundance') %>% 
  group_by(site, plot, Year, taxa) %>% 
  summarize(abundance=sum(abundance)) %>% 
  ungroup() %>% 
  # filter(taxa %in% SVsToPlot) %>% 
  left_join(metadata2) %>%
  mutate(NormAbundance=log10(abundance+0.0001)) %>%  # do a log10 transformation after adding a 0.01% pseudocount. Could also add 1 read before transformation to percent
  mutate(taxa2=fct_reorder(taxa, NormAbundance, .desc=T)) 
# mutate(Feature=paste(Feature.ID, Taxon)) %>%
# mutate(Feature=gsub("[kpcofgs]__", "", Feature)) # trim out leading text from taxonomy string

# write.csv(SVphylaHeatMap, 'C:\\Users\\kjkomatsu\\OneDrive - UNCG\\GMDR\\soils ms derived data\\GMDR_microbial_abundance_2019-2022.csv', row.names=F)


fSVphylaHeatMap$grazing_treatment=factor(fSVphylaHeatMap$grazing_treatment,levels=c('destock', 'stable', 'heavy'))
fSVphylaHeatMap <- fSVphylaHeatMap[order(fSVphylaHeatMap$NormAbundance),]

#fungal by drought TB
ggplot(data=subset(fSVphylaHeatMap, site=='TB' & grepl('Fungi',taxa2)), aes(x=as.factor(rainfall_reduction), y=taxa2, fill=NormAbundance)) +
  geom_tile() +
  # facet_grid(~`Year`, scales="free_x") +
  # theme_q2r() +
  # theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_viridis_c(name="log10(% Abundance)") +
  scale_y_discrete(limits=rev) +
  # scale_x_discrete(labels=c('destock', 'supplement', 'maintain')) +
  xlab('Rainfall Reduction (%)') + ylab('Taxa')
# ggsave('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_GMDR_soils\\figures\\Fig7_TB_fungalRainfall_heatmap.png', width=12, height=8, units='in', dpi=300, bg='white')

#fungal by grazing TB
ggplot(data=subset(fSVphylaHeatMap, site=='TB' & grepl('Fungi',taxa2)), aes(x=as.factor(grazing_treatment), y=taxa2, fill=NormAbundance)) +
  geom_tile() +
  # facet_grid(~`Year`, scales="free_x") +
  # theme_q2r() +
  # theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_viridis_c(name="log10(% Abundance)") +
  scale_y_discrete(limits=rev) +
  scale_x_discrete(labels=c('destock', 'stable', 'heavy')) +
  xlab('Grazing Treatment') + ylab('Taxa')
# ggsave('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_GMDR_soils\\figures\\Fig7_TB_fungalGrazing_heatmap.png', width=12, height=8, units='in', dpi=300, bg='white')


#fungal by rainfall FK
ggplot(data=subset(fSVphylaHeatMap, site=='FK' & grepl('Fungi',taxa2)), aes(x=as.factor(rainfall_reduction), y=taxa2, fill=NormAbundance)) +
  geom_tile() +
  facet_grid(~`Year`) +
  # theme_q2r() +
  # theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_viridis_c(name="log10(% Abundance)") +
  scale_y_discrete(limits=rev) +
  # scale_x_discrete(labels=c('destock', 'supplement', 'maintain')) +
  xlab('Rainfall Reduction (%)') + ylab('Taxa')
# ggsave('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_GMDR_soils\\figures\\FigS15_FK_fungalRainfall_heatmap.png', width=12, height=8, units='in', dpi=300, bg='white')



### bacterial taxa specific responses
bSVsToPlot <- bSVrelative %>% #find the average abundance of a SV
  # rename('Fungi::Other'='Fungi::NA') %>% 
  pivot_longer(cols=c('d__Bacteria::Abditibacteriota':'d__Bacteria::WPS-2'), names_to='taxa', values_to='abundance') %>% 
  group_by(site, plot, Year, taxa) %>% 
  summarize(abundance=sum(abundance)) %>% 
  ungroup() %>% 
  group_by(taxa) %>% 
  summarize(mean_abundance=mean(abundance)) %>% 
  ungroup() %>% 
  arrange(desc(mean_abundance)) %>%
  top_n(20, mean_abundance) %>%
  pull(taxa) #extract only the names from the table

metadata2 <- metadata %>% 
  select(-type) %>% 
  unique()

bSVphylaHeatMap <- bSVrelative %>%
  # rename('Fungi::Other'='Fungi::NA') %>%
  pivot_longer(cols=c('d__Bacteria::Abditibacteriota':'d__Bacteria::WPS-2'), names_to='taxa', values_to='abundance') %>% 
  group_by(site, plot, Year, taxa) %>% 
  summarize(abundance=sum(abundance)) %>% 
  ungroup() %>% 
  # filter(taxa %in% SVsToPlot) %>% 
  left_join(metadata2) %>%
  mutate(NormAbundance=log10(abundance+0.0001)) %>%  # do a log10 transformation after adding a 0.01% pseudocount. Could also add 1 read before transformation to percent
  mutate(taxa2=fct_reorder(taxa, NormAbundance, .desc=T)) 
# mutate(Feature=paste(Feature.ID, Taxon)) %>%
# mutate(Feature=gsub("[kpcofgs]__", "", Feature)) # trim out leading text from taxonomy string

# write.csv(SVphylaHeatMap, 'C:\\Users\\kjkomatsu\\OneDrive - UNCG\\GMDR\\soils ms derived data\\GMDR_microbial_abundance_2019-2022.csv', row.names=F)


bSVphylaHeatMap$grazing_treatment=factor(bSVphylaHeatMap$grazing_treatment,levels=c('destock', 'stable', 'heavy'))

#fungal by grazing TB
ggplot(data=subset(bSVphylaHeatMap, site=='TB' & Year>2019), aes(x=as.factor(grazing_treatment), y=taxa2, fill=NormAbundance)) +
  geom_tile() +
  facet_grid(~`Year`) +
  # theme_q2r() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  scale_fill_viridis_c(name="log10(% Abundance)") +
  scale_y_discrete(limits=rev) +
  scale_x_discrete(labels=c('destock', 'stable', 'heavy')) +
  xlab('Rainfall Reduction (%)') + ylab('Taxa')
# ggsave('C:\\Users\\kjkomatsu\\OneDrive - UNCG\\manuscripts\\first author\\2024_GMDR_soils\\figures\\FigS13_TB_bacterialGrazing_heatmap.png', width=12, height=8, units='in', dpi=300, bg='white')
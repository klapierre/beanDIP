#####################################################################
##  beanDIP_greenhouse_phylogeneticDiversity.R
##  Phylogenetic diversity for all treatment combinations for Randall greenhouse
##  and Parker common garden
##
##  Author: Kimberly Komatsu
#####################################################################

library(phytools)
library(picante)
library(grid)
library(vegan)
# library(BiodiversityR)
library(tidyverse)

setwd('C:\\Users\\kjkomatsu\\Smithsonian Dropbox\\Kimberly Komatsu\\bean_dip_2018-2024\\greenhouse trials\\phylogeneticDistance_greenhouseCommonGarden')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))

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

#####################
# tree management
#####################

#mr bayes tree: 40 million runs
treeNexus<-read.nexus(file='012023_final_sequences_aligned.nex.con.tre')

#makes tree into a .phylo file instead of .newick or .tre
tree<-as.phylo(treeNexus) 

#root tree with Otu152
treeR<-root(tree, "Otu152", r=T)

#plot trees
plot.phylo(treeR, use.edge.length=T)
plot.phylo(treeR, use.edge.length=F)



#####################
# Randall Diversity Greenhouse Experiment
#####################

#select strains from experiment
treeRandall<-keep.tip(treeR, c('Otu006', 'Otu008', 'Otu043', 'Otu045', 'Otu066', 'Otu090', 'Otu119', 'Otu123'))
plot.phylo(treeRandall, show.node.label=T)


#making a tree with only Brady japonicum strains, excluding Australian Ulex
BjITStree<-drop.tip(concITStree, c('ITS_026', 'ITS_028', 'ITS_029', 'ITS_031', 'ITS_032', 'ITS_033', 'ITS_011', 'ITS_021', 'Rhizobium_leguminosarum_X01z_ITS'))
plot.phylo(BjITStree, use.edge.length=T)
nodelabels(frame='none', cex=0.8)


#making a tree with only Brady japonicum strains, excluding Australian Ulex and ITS_030 (long branch length)
BjITStree<-drop.tip(concITStree, c('ITS_026', 'ITS_028', 'ITS_029', 'ITS_031', 'ITS_032', 'ITS_033', 'ITS_011', 'ITS_021', 'Rhizobium_leguminosarum_X01z_ITS', 'ITS_030'))
plot.phylo(BjITStree, use.edge.length=T)
nodelabels(frame='none', cex=0.8)


#making a tree with only Brady japonicum strains from the Bay Area
BjBayAreaITStree<-drop.tip(BjITStree, c('ITS_011', 'ITS_021', 'ITS_027')) 
plot.phylo(BjBayAreaITStree, use.edge.length=T)
nodelabels(frame='none', cex=0.8)



#####################
# strain data management
#####################
#read in nodule data
ITSnodules <- read.csv('strain data\\La Pierre_invasion molecular manuscript_strain information_01252016.csv')%>%
  select(plant_species, plant_status, nodule_ID, ITS_OTU_97, concatenated_OTU_98)

#create an interaction matrix of strains for each plant species
ITSinteractionMatrix <- ITSnodules%>%
  select(plant_species, plant_status, nodule_ID, ITS_OTU_97, concatenated_OTU_98)%>%
  filter(ITS_OTU_97!='', concatenated_OTU_98!='')%>%
  select(-concatenated_OTU_98)%>%
  mutate(interaction=1)%>%
  spread(key=ITS_OTU_97, value=interaction, fill=0)

#subset out only Bay Area and brady japonicum strains
ITSbjBayAreaInteractionMatrix <- ITSinteractionMatrix%>%
  select(-ITS_032, -ITS_011, -ITS_021, -ITS_027)%>%
  filter(nodule_ID!='X01z', plant_species!='Ulex europaeus - Australia', plant_species!='Spartium junceum - Italy', plant_species!='Ulex europaeus - Portugul')%>%
  #get summary interaction matrix (sum of interactions by species)
  gather(key=ITS_OTU_98, value=interaction, -plant_species, -plant_status, -nodule_ID)%>%
  group_by(plant_status, plant_species, ITS_OTU_98)%>%
  summarise(interaction=sum(interaction))%>%
  spread(key=ITS_OTU_98, value=interaction)





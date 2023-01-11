
setwd("C:/Users/Sarah Alley/Dropbox (Smithsonian)/bean_dip_2018-2024/lab work/Sanger Sequencing/sequence data/all_samples_consolidated/all_nifH_as_of_12_2022")

OTUs<-read.csv("12292022_OTUs.csv")

library(tidyverse)

install.packages("vegan")
library(vegan)

#Rank abundance curve

#Making data long form 
#Useful code for separating out Site/Plot/Drought Treatment/Nodule number into separate columns

cleanOTUs<-OTUs%>%pivot_longer(!OTU, names_to="delete", values_to="sample")%>%
  select(-delete)%>%
  filter(sample!="")%>%
  separate(col=sample, into=c("plant", "nodule"), sep="-")%>%
  separate(col=plant, into=c("site_plot", "drought"), sep="(?<=\\d)(?=[a-z]?)")%>%
  mutate(plot=readr::parse_number(site_plot), site=str_remove(site_plot, "\\d+"))%>%
  select(-site_plot)%>%
  mutate(drt_trt=ifelse(drought%in%c("a", "b", "c", "d", "e", "f", "g", "h", "i"),"ambient", "drought"))%>%
  mutate(id=paste(site, plot, drought, nodule, sep="_"))


#Need to put data back into wide form after separating samples into different columns 
widecleanOTUs <-cleanOTUs%>%
  mutate(presence=1)%>%
  pivot_wider(names_from=OTU, values_from=presence, values_fill=0)

species_accumulation<- specaccum(widecleanOTUs[, 7:188])

plot(species_accumulation)

#Want all ambient, all droughted, then by four different sites (ambient and droughted)

#Ambient samples

ambientcleanOTUs <-filter(widecleanOTUs, drt_trt == "ambient")
species_accumulation_ambient<- specaccum(ambientcleanOTUs[, 7:188])
plot(species_accumulation_ambient)
  
#Droughted samples 

droughtedcleanOTUs <-filter(widecleanOTUs, drt_trt == "drought")
species_accumulation_droughted<- specaccum(droughtedcleanOTUs[, 7:188])
plot(species_accumulation_droughted)

#Clarksville samples

CvillecleanOTUs <-filter(widecleanOTUs, site == "C")
species_accumulation_Cville<- specaccum(CvillecleanOTUs[, 7:188])
plot(species_accumulation_Cville)

#Keedysville samples 

KvillecleanOTUs <-filter(widecleanOTUs, site == "K")
species_accumulation_Kville<- specaccum(KvillecleanOTUs[, 7:188])
plot(species_accumulation_Kville)

#Wye samples

WyecleanOTUs <-filter(widecleanOTUs, site == "W")
species_accumulation_Wye<- specaccum(WyecleanOTUs[, 7:188])
plot(species_accumulation_Wye)

#Poplar Hill samples

PopHillcleanOTUs <-filter(widecleanOTUs, site == "P")
species_accumulation_PopHill<- specaccum(PopHillcleanOTUs[, 7:188])
plot(species_accumulation_PopHill)

#Within sites: everything once, little overlap
#Among sites: some overlap (OTU 1)





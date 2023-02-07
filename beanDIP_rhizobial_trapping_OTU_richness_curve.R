
setwd("C:/Users/Sarah Alley/Dropbox (Smithsonian)/bean_dip_2018-2024/lab work/Sanger Sequencing/sequence data/consensus_mothur_consolidated/OTUs")

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
species_accumulation_ambient<- specaccum(ambientcleanOTUs[, 7:188], method="collector")
plot(species_accumulation_ambient)
  
#Droughted samples 

droughtedcleanOTUs <-filter(widecleanOTUs, drt_trt == "drought")
species_accumulation_droughted<- specaccum(droughtedcleanOTUs[, 7:188])
plot(species_accumulation_droughted)

#Clarksville samples

CvillecleanOTUs <-filter(widecleanOTUs, site == "C")
species_accumulation_Cville<- specaccum(CvillecleanOTUs[, 7:188], method="collector")
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


#Look at new data
#Make these curves for each site separated by droughted by undroughted
#Want about the same number of samples for each of these curves
#If we need more samples, target specifically ones that will even things out
#Finish prepping for sequencing and send off ones I have already done in the fridge
#Let the group know details of the results 

#02/03/2023, repeating analyses with all 599 (or 598?) samples we have sequences for 









#In short form, not all the OTUs show up along the top, but works when converted to long form in next step 
OTUs2<-read.csv("OTUs_01312023.csv")

cleanOTUs2<-OTUs2%>%pivot_longer(!OTU, names_to="delete", values_to="sample")%>%
  select(-delete)%>%
  filter(sample!="")%>%
  separate(col=sample, into=c("plant", "nodule"), sep="-")%>%
  separate(col=plant, into=c("site_plot", "drought"), sep="(?<=\\d)(?=[a-z]?)")%>%
  mutate(plot=readr::parse_number(site_plot), site=str_remove(site_plot, "\\d+"))%>%
  select(-site_plot)%>%
  mutate(drt_trt=ifelse(drought%in%c("a", "b", "c", "d", "e", "f", "g", "h", "i"),"ambient", "drought"))%>%
  mutate(id=paste(site, plot, drought, nodule, sep="_"))

#Between this step and the last step, 4 observations are lost... 
widecleanOTUs2 <-cleanOTUs2%>%
  mutate(presence=1)%>%
  pivot_wider(names_from=OTU, values_from=presence, values_fill=0)


#Ask Kim what the 7:188 means 
species_accumulation2<- specaccum(widecleanOTUs2[, 7:188])

plot(species_accumulation2)

#Ambient samples

#Using filter doesn't work for some reason, had to switch to subset
#ambientcleanOTUs2 <-filter(widecleanOTUs2, drt_trt == "ambient")
#R reports the correct number of samples for ambient 

ambientcleanOTUs2 <- subset(widecleanOTUs2, drt_trt == "ambient")
species_accumulation_ambient2<- specaccum(ambientcleanOTUs2[, 7:188], method="collector")
plot(species_accumulation_ambient2)

#*Have reached a plateau with ambient samples!341 samples are ambient  

#Droughted samples 
#R reports the correct number of samples for drought

droughtedcleanOTUs2 <-subset(widecleanOTUs2, drt_trt == "drought")
species_accumulation_droughted2<- specaccum(droughtedcleanOTUs2[, 7:188])
plot(species_accumulation_droughted2)

#Clarksville samples

#R does not report the correct number of samples, there are 169 samples from Clarksville total, 68 droughted, 101 ambient
#Some of the sites have a space before or after the letter, for example, Clarksville could be C space instead of just C
#Be sure new dataframes give number of samples you're expecting 
  
#All Clarksville 
CvillecleanOTUs2 <- widecleanOTUs2 %>% filter(site %in% c("C", " C"))
species_accumulation_Cville2<- specaccum(CvillecleanOTUs2[, 7:188], method="collector")
plot(species_accumulation_Cville2)

#Clarksville ambient 
ambientCvillecleanOTUs2 <- widecleanOTUs2 %>% filter(site %in% c("C", " C") & drt_trt == "ambient")
ambient_species_accumulation_Cville2<- specaccum(ambientCvillecleanOTUs2[, 7:188], method="collector")
plot(ambient_species_accumulation_Cville2)

#Clarksville droughted
droughtedCvillecleanOTUs2 <- widecleanOTUs2 %>% filter(site %in% c("C", " C") & drt_trt == "drought")
drought_species_accumulation_Cville2<- specaccum(ambientCvillecleanOTUs2[, 7:188], method="collector")
plot(drought_species_accumulation_Cville2)

#Keedysville samples 
KvillecleanOTUs2 <- widecleanOTUs2 %>% filter(site %in% c("K", " K"))
species_accumulation_Kville2<- specaccum(KvillecleanOTUs2[, 7:188])
plot(species_accumulation_Kville2)

#Keedysville ambient 
ambientKvillecleanOTUs2 <- widecleanOTUs2 %>% filter(site %in% c("K", " K") & drt_trt == "ambient")
ambient_species_accumulation_Kville2<- specaccum(ambientKvillecleanOTUs2[, 7:188])
plot(ambient_species_accumulation_Kville2)

#Keedysville droughted
droughtKvillecleanOTUs2 <- widecleanOTUs2 %>% filter(site %in% c("K", " K", "K ") & drt_trt == "drought")
drought_species_accumulation_Kville2<- specaccum(droughtKvillecleanOTUs2[, 7:188])
plot(drought_species_accumulation_Kville2)

#Wye samples

WyecleanOTUs2 <-filter(widecleanOTUs2, site %in% c("W", " W", "W "))
species_accumulation_Wye2<- specaccum(WyecleanOTUs2[, 7:188])
plot(species_accumulation_Wye2)

#Wye ambient 
ambientWyecleanOTUs2 <-filter(widecleanOTUs2, site %in% c("W", " W", "W ") & drt_trt == "ambient")
ambient_species_accumulation_Wye2<- specaccum(ambientWyecleanOTUs2[, 7:188])
plot(ambient_species_accumulation_Wye2)

#Wye droughted
droughtWyecleanOTUs2 <-filter(widecleanOTUs2, site %in% c("W", " W", "W ") & drt_trt == "drought")
drought_species_accumulation_Wye2<- specaccum(droughtWyecleanOTUs2[, 7:188])
plot(drought_species_accumulation_Wye2)

#Poplar Hill samples
PopHillcleanOTUs2 <-filter(widecleanOTUs2, site %in% c("P", " P", "P "))
species_accumulation_PopHill2<- specaccum(PopHillcleanOTUs2[, 7:188])
plot(species_accumulation_PopHill2)

#Poplar Hill ambient
ambientPopHillcleanOTUs2 <-filter(widecleanOTUs2, site %in% c("P", " P", "P ")& drt_trt == "ambient")
ambient_species_accumulation_PopHill2<- specaccum(ambientPopHillcleanOTUs2[, 7:188])
plot(ambient_species_accumulation_PopHill2)

#Poplar Hill droughted
droughtPopHillcleanOTUs2 <-filter(widecleanOTUs2, site %in% c("P", " P", "P ")& drt_trt == "drought")
drought_species_accumulation_PopHill2<- specaccum(droughtPopHillcleanOTUs2[, 7:188])
plot(drought_species_accumulation_PopHill2)

#######################################################################################
##  field_site_map.R: Creates map of Maryland with 4 soy field sites 
##
##  Author: Kelsey McGurrin
##
#######################################################################################

####setup####
library(tidyverse)
library(readxl)
library(sf) #data frames with spatial features
library(tigris) #us census data downloads


# working directory path (add yours if different)
setwd("~/Dropbox/bean_dip_2018-2024/field trials/plot maps")

# read in site coordinates as shapefile
# use NAD83 as CRS because that's what tigris package will use for state outline
# https://www.nceas.ucsb.edu/sites/default/files/2020-04/OverviewCoordinateReferenceSystems.pdf
sites <- read_excel("site_coordinates.xlsx") %>%
  filter(year==2019) %>%
  select(name,long_dd,lat_dd) %>%
  st_as_sf(coords = c("long_dd", "lat_dd"), crs = 4269)
sites

# get outline of maryland
md <- counties(state = "maryland", cb = TRUE, resolution = "500k", year = 2022) 

# basic plot
ggplot() +
  geom_sf(data=md, fill="white") +
  geom_sf(data = sites, size = 4, color="darkgreen")+
  theme_void()

# with labels for site names
ggplot() +
  geom_sf(data=md) +
  geom_sf_label(data=sites,aes(label = sites$name),nudge_y = 0.15)+
  geom_sf(data = sites, size = 4, color="darkgreen")+
  theme_bw()+
  theme(axis.title=element_blank())

ggsave("site_map_labels.png",width = 6,height = 4,units=c("in"))  

# regional map
midA_states<-counties(cb=T) %>% filter(STATEFP %in% c("24","42","54","51","10","11"))
dc<-counties(cb=T) %>% filter(STATEFP =="11")

ggplot()+
  geom_sf(data=midA_states,fill="#f0f0f0")+
  geom_sf(data=dc,fill="#f0f0f0")+
  geom_sf(data=md,fill="#bdbdbd")+
  geom_sf_label(data=sites,aes(label = sites$name),nudge_y = 0.15)+
  geom_sf(data=sites, size = 4, fill="darkgreen",shape = 21, colour = "black")+
  coord_sf(xlim=c(-80,-74.5), ylim=c(37.5, 40.5), expand=F)+
  theme_bw()+
  theme(axis.title=element_blank(),panel.background = element_rect(fill = "#deebf7"))

ggsave("regional_site_map_labels.png",width = 6,height = 4,units=c("in"))  


# all of contig US
st <- states(cb=T) 
st_contig <- st %>% filter(!(NAME %in% c("Alaska", "Hawaii","Commonwealth of the Northern Mariana Islands",
                                         "Puerto Rico","Guam", "American Samoa","United States Virgin Islands")))

ggplot()+
  geom_sf(data=st_contig)+
  theme_bw()

ggsave("contig_us.png",width = 6,height = 4,units=c("in"))  


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

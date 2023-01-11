################################################################################
##  clean_photosynq.R: Cleaning photosynq data.
##
##  Author: Kelsey McGurrin
################################################################################


####setup####
library(dplyr)
library(lubridate)
setwd("~/Dropbox/bean_dip_2018-2024/field trials/data")

####2019####
data <- read.csv("raw_data/2019_Leaf_multispec.csv")

#make changes according to "notes" column
#remove row from dataset
#1315404 1313013
datab<-data[!(data$ID %in% c(1315404,1313013)),] 

#change indiv to 5
#1308606 1308604 1308598 1308595 1308584 1308581 1308579 1308573 1308567	1308565	1308563	1308558	1308571	1308569	1308589 
datab$Individual[datab$ID %in% c(1308606,1308604, 1308598, 1308595, 1308584, 1308581, 1308579, 1308573, 1308567, 1308565, 1308563,	1308558,	1308571,	1308569,	1308589)]<-5

#change indiv to 6
#1308559	1308553	1308597 1308591 1308583 1308575 1386911
datab$Individual[datab$ID %in% c(1308559,	1308553,	1308597, 1308591, 1308583, 1308575, 1386911)]<-6

#change indiv to 4
#1308592 
datab$Individual[datab$ID==1308592]<-4

#change indiv to 3
#1386912 
datab$Individual[datab$ID==1386912]<-3

#change indiv to 1
#1312993
datab$Individual[datab$ID==1312993]<-1

#change variety to 83
#1308590 1308570	1308572 1308571	1308569	1308589 
datab$Variety[datab$ID %in% c(1308590, 1308570,	1308572, 1308571,	1308569,	1308589)]<-83

#change plot 3 to 1 because extra 67-3-6/missing 67-1-6
datab$Plot[datab$ID==1310003]<-1

#67-3-1 was measured twice. no note but keep second
datab<-datab[!(data$ID==1310163),] 

#remove column from dataset because no difference between variety and variety 2
#variety 2
datab$Variety.1<-NULL

#output cleaned file
write.csv(datab,file="clean_data/clean_photosynq_2019.csv",row.names = F)

####2020####
data <- read.csv("raw_data/2020_Leaf_multispec.csv")

#create shortcut view of issues
notes<-tibble(data$ID,data$note)

#make changes according to "notes" column
datab<-data %>%
  select(-c(note,Issues))

#1790006 actually variety 70 at clarksville
datab$Variety[datab$ID==1790006]<-70
datab$Site[datab$ID==1790006]<-"Clarksville"

#1786113 actually variety 5
datab$Variety[datab$ID==1786113]<-5

#1790022 actually indiv 6
datab$Individual[datab$ID==1790022]<-6

#1790023 actually indiv 4
datab$Individual[datab$ID==1790023]<-4
  
#1787843 and 1786114 actually indiv 3
datab$Individual[datab$ID==1787843]<-3
datab$Individual[datab$ID==1786114]<-3

#1786102 actually indiv 2
datab$Individual[datab$ID==1786102]<-2


#output cleaned file
write.csv(datab,file="clean_data/clean_photosynq_2020.csv",row.names = F)

####2021####
data <- read.csv("raw_data/2021_Leaf_multispec.csv")

#drop extra columns from output
datab<-data %>%
  select(-c(note,Issues))

#make changes according to "notes" column
data_notes<-data %>%
  select(ID,time,Individual,Plot,Site,Variety,note) %>%
  filter(note!="")

#2255129 actually variety 58
datab$Variety[datab$ID==2255129]<-58
  
#2247503 actually indiv 3
datab$Individual[datab$ID==2247503]<-3
  
#2247501 actually indiv 1
datab$Individual[datab$ID==2247501]<-1

#2247490 actually site clarksville
datab$Site[datab$ID==2247490]<-"Clarksville"

#2247479 actually indiv 1
datab$Individual[datab$ID==2247479]<-1

#some plants at keedysville dying/dead - exclude?
#2258556 & 2258555 very sad
#2258554 & 2258553 sad

# wye has entries on two days: 8/16 and 8/18. everything from the 18th should be PH
datab$day<-parse_date_time(as.character(datab$time),"mdy_HM")
datab$Site[date(datab$day)=="2021-08-18"]<-"Poplar Hill"

#output cleaned file
write.csv(datab,file="clean_data/clean_photosynq_2021.csv",row.names = F)


library(dplyr)
setwd("~/Dropbox/bean_dip_2018-2024/data")
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

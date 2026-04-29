

setwd("C:/Users/nh1087/OneDrive - USNH/Documents/NECC/Diet Data/")


library(tidyverse)
library(ggplot2)
library(readr)



#Editing the table from Steimle et al. 1985, never have to do again
#preyEDs<-read.csv("Steimle1985_ED.csv")%>%
#  separate(col="Taxa",into=c("Taxa","N_combustions","Ash_perc","H2O_perc",
#                             "Dry_KJg_mean","Dry_KJg_SD","Ash_KJg","Wet_KJg","Shell_perc"),
#           sep=" ")
#
#write.csv(file="Steimle1985_ED_2.0.csv",preyEDs)


preyEDs<-read.csv("Steimle1985_ED_2.0.csv")

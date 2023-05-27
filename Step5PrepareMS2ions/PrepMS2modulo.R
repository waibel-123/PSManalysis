# Prepare MS2 ion m/z list as csv. data for Mass Defect (Modulo 1) analysis

#setwd
setwd("x")

#library
library(ggplot2)
library(dplyr)
library(ggExtra)


#load prepared Rda dataframe (csv too big)
load("Water_3_DF_2.Rda")

DF_final <- DF_2
rm(DF_2)

length(unique(DF_final$mz))

# data prep for those spectra that are different ----
"filter full set below to DF_size"

df <- filter(DF_final, Category =="Assigned_QC" | Category =="Assigned_contams")
length(unique(df$mz))
df <- df[,c(9,10)]
write.csv(df, file="Water_3_Assigned.csv", row.names = FALSE)

df <- filter(DF_final, Category =="Unassigned")
length(unique(df$mz))
df <- df[,c(9,10)]
write.csv(df, file="Water_3_Unassigned.csv", row.names = FALSE)

df <- filter(DF_final, Category =="Unmatched")
length(unique(df$mz))
df <- df[,c(9,10)]
write.csv(df, file="Water_3_Unmatched.csv", row.names = FALSE)
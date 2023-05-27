# Prepare data tables

# working dir
setwd("x")

#library
library(tidyr)
library(dplyr)
library(matrixStats)

# Read in data tables ----
list.files(path = "x", pattern = ".csv")

DF <- read.csv("Water_3.csv", sep = ";", header = TRUE)


# Check correct naming of headers and categories ----
DF$Category <- as.factor(DF$Category)

str(DF)
summary(DF$Category)
summary(DF)

# Split and remove doubly assigned, where Assigned have precedent over Anchored and Anchored over Unassigned ----

DF_Assigned <- filter(DF, Category=="Assigned_QC" | Category=="Assigned_contams")
summary(DF_Assigned)

DF_Anchored <- filter(DF, Category=="Anchored_QC" | Category=="Anchored_contams")
DF_Anchored_new <- DF_Anchored[!(DF_Anchored$pep_query %in% DF_Assigned$pep_query),]

DF_Unmatched <- filter(DF, Category=="Unmatched")
DF_Unmatched_new <- DF_Unmatched[!(DF_Unmatched$pep_query %in% DF_Assigned$pep_query), ]

DF1 <- rbind(DF_Assigned, DF_Anchored_new, DF_Unmatched_new)

DF2 <- filter(DF, Category=="Unassigned")
DF2$Category <- droplevels(DF2$Category)
summary(DF2)

DF2_new <- DF2[!(DF2$pep_query %in% DF1$pep_query), ]


A <- Reduce(base::union, list(DF2_new$pep_query,
                              DF1$pep_query))

DF_new <- rbind(DF1, DF2_new)

# check summary ----
length(DF_new$pep_query)
length(unique(DF_new$pep_query))

summary(DF_new)

# export as prep table ----

write.csv(DF_new, file = "Water_3_prep.csv",
          row.names = FALSE)

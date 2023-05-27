# Extract MS2 ions from string

#wd
setwd("x")

#libraries

library(tidyr)
library(dplyr)
library(stringi)
library(ggplot2)
library(matrixStats)
library(ggExtra)

#library(interp)
#library(gghighlight)

# File
Water_3 <- read.csv("Water_3_prep.csv", header = TRUE)
DF_input <- Water_3
rm(Water_3)


length(unique(DF_input$pep_query))


# Get CID fragment ions from strings ----
# Split column with ion species, can be hundreds of new columns 
DList <- strsplit(DF_input$StringIons1, ",")

DF <- stri_list2matrix(DList, byrow = TRUE)
#DF [1,] # check length of select spectra
dim(DF)

# put new columns of ions together with prev DF & remove pep scan title if not needed for analysing scan
DF_Water <- cbind(DF_input[,c(1,2,3,4,6,7,8)], DF)

str(DF_Water)

# Pivot columns to long format for further string split
A <- dim(DF_Water)
A[2]

DF_Water_long <- pivot_longer(DF_Water, cols = 8:A[2], names_to = "Fragment", values_to = "String", values_drop_na = TRUE)
str(DF_Water_long)

# final string split into mz and intensity value for each fragment ion
V_string_split <- strsplit(DF_Water_long$String, ":")
DF_string_split <- stri_list2matrix(V_string_split, byrow = TRUE)
colnames(DF_string_split) <-c("mz","Intensity")

DF_final <- cbind(DF_Water_long[,1:8], DF_string_split)
str(DF_final)

DF_final$Category <- as.factor(DF_final$Category)
DF_final$Fragment <- as.integer(DF_final$Fragment)
DF_final$mz <- round(as.numeric(DF_final$mz), digits = 4) # precision on this instrument not more than 4 digits, accuracy to 0.05 Da
DF_final$Intensity <- as.numeric(DF_final$Intensity)
DF_final$charge <- gsub('[+]',"", DF_final$charge)
DF_final$charge <- as.integer(DF_final$charge)

dim(DF_final)
str(DF_final)
summary(DF_final)
length(unique(DF_final$pep_query))

DF_final$Category <- factor(DF_final$Category, levels = c("Assigned_QC",  
                                                          "Anchored_QC",
                                                          "Unassigned", 
                                                    "Unmatched",
                                                    "Assigned_contams",
                                                    "Anchored_contams"
                                                    ))

DF_final$Retention.time.range <- DF_final$Retention.time.range/60

rm(DF)
rm(DList)
rm(DF_string_split)
rm(V_string_split)
rm(Water_3)
rm(DF_Water)
rm(DF_Water_long)
rm(A)

# Seperate data into LC gradient phases ----


DF_final$Phase <- ifelse(DF_final$Retention.time.range>0 & DF_final$Retention.time.range <=30, "Phase 1A",
                         ifelse(DF_final$Retention.time.range>30 & DF_final$Retention.time.range <=60, "Phase 1B",
                                ifelse(DF_final$Retention.time.range>60 & DF_final$Retention.time.range <=67.5, "Phase 2A",
                                       ifelse(DF_final$Retention.time.range>67.5 & DF_final$Retention.time.range <=75, "Phase 2B", 
                                              ifelse(DF_final$Retention.time.range>75 & DF_final$Retention.time.range <=80, "Phase 3",
                                                     ifelse(DF_final$Retention.time.range>80, "Phase 4", "Other"))))))




# # Calculate parent ion (= pep_query) specific mass differences & plot density diagram ----


DF_2 <- DF_final %>% group_by(pep_query) %>%  arrange( mz, .by_group = TRUE) %>% 
  mutate(mass_diff = mz-lag(mz), .keep="all") %>% ungroup()

str(DF_2)

summary(subset(DF_2, !is.na(mass_diff))$mass_diff)
hist(log(subset(DF_2, !is.na(mass_diff))$mass_diff))

print(length(which(DF_2$mass_diff == 2.0157)))

max(subset(DF_2, !is.na(mass_diff))$mass_diff)


# Save files 
save(DF_2, file= "Water_3_DF_2.Rda")


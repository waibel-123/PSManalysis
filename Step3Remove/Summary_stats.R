# LC-ESI-QqTOF MS/MS spectral analysis

# working directory
setwd("x")

#library
library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)

# Import Data from stackoverflow and Claus Wilke hp 
file_names <- list.files(path = "x", pattern = "prep.csv")
# DF <- do.call(rbind, sapply(file_names, read.csv(header=TRUE)))
data_path <- "x"

data <- tibble(Sample=file_names) %>% 
                          mutate(file_contents = map(Sample, ~read_csv(file.path(data_path, .)))
                                )                   
  
DF <- unnest(data, cols = c(1,2))

DF$Sample <- gsub(".csv","", DF$Sample)

DF$Sample <- gsub("_prep","", DF$Sample)

DF$Sample <- as.factor(DF$Sample)

DF$Category <- as.factor(DF$Category)

DF$Retention.time.range <- DF$Retention.time.range/60

DF$charge <- gsub("[+]","", DF$charge) %>% as.integer()

summary(DF)

# Prepare dataframe for plotting spectral counts panel ----
df <- DF[ , c(1,2,7)]

Summary <- df %>% group_by(Sample,Category) %>% summarise(n(), .groups="keep")
write.csv(Summary, file = "Spectral_counts_summary.csv")

## Plot samples main run (only NCBI FDR 1%) remove 4x first run samples ----
unique(df$Sample)

DF_fil <- filter(df, Sample!="NaPPi_1_old", Sample!="EDTA_2_old", Sample!="Water_3_old", Sample!="NaPPi_control_7_old")


DF_fil <- droplevels(DF_fil)
unique(DF_fil$Sample)

#DF_fil$Sample <- factor(DF_fil$Sample, levels = c(
#                                         "EDTA_control_8", "Water_control_9", "NaPPi_1",
#                                          "EDTA_2", "Water_3", "NaPPi_control_7", "NaPPi_10",
#                                          "EDTA_11", "Water_12", "NaPPi_19", "EDTA_20",
#                                          "Water_21"))

DF_fil$Sample <- factor(DF_fil$Sample, levels = c("Water_control_9", "Water_3", "Water_12", "Water_21",
                                                  "NaPPi_control_7", "NaPPi_1", "NaPPi_10", "NaPPi_19",
                                                  "EDTA_control_8", "EDTA_2", "EDTA_11", "EDTA_20"))
                                                  
                                                 
DF_fil$Category <- as.character(DF_fil$Category)
DF_fil <- DF_fil %>%
  mutate(., Category = with(., case_when(
    (Category=="Assigned_QC") ~ "Other",
    (Category=="Anchored_QC") ~ "Other",
    (Category=="Assigned_contams") ~ "Other",
    (Category=="Anchored_contams") ~ "Other",
    TRUE                      ~ Category
  )
  )
  )

unique(DF_fil$Category)

DF_fil$Category <- factor(DF_fil$Category, levels = c("Unassigned", 
                                                      "Unmatched", 
                                                      "Other"))

V_label <- c( 
             "Water control", "Water 1", "Water 2", "Water 3",
             "NaPPi control", "NaPPi 1", "NaPPi 2", "NaPPi 3",
             "EDTA control", "EDTA 1", "EDTA 2", "EDTA 3"
              )

names(V_label) <- c("Water_control_9", "Water_3", "Water_12", "Water_21",
                    "NaPPi_control_7", "NaPPi_1", "NaPPi_10", "NaPPi_19",
                    "EDTA_control_8", "EDTA_2", "EDTA_11", "EDTA_20")


PlotHist_all <-   ggplot(DF_fil, aes(x=Retention.time.range, colour=Category)) + 
                  geom_histogram(aes(y=..count.., fill=Category), binwidth = 2.5, position = "dodge")+
                  geom_density(aes(y=2.5*..count..), show.legend = TRUE, linetype="solid", trim=TRUE)+
                  facet_wrap(~Sample, scales = "fixed", ncol=4, labeller = labeller(Sample=V_label))+ 
                  theme_bw()+
                  #xlim(0,100)+
                  #scale_x_continuous(breaks = seq(0, 100, by = 10), minor_breaks = seq(0, 100,2.5), expand = c(0, 0))+
                  scale_x_continuous(limits=c(0,100), breaks = c(0,10,20,30,40,50,60,70,80,90,100), labels = c(0,10,20,30,40,50,60,70,80,90,100))+
                  #geom_vline(xintercept = c(60,75,80,85,87.5,100), linetype='dashed', color='black')+
                  xlab("Retention time (min)")+
                  ylab("Spectral counts \n")+
                  #scale_color_brewer(palette="Set1")+
                  theme(legend.position = "none") +
                  #theme(legend.title = element_blank())+
                  scale_fill_manual(values=c("#56B4E9", "#9966CC", "1", "#999999", "#FF00FF"))+
                  scale_color_manual(values=c("#56B4E9", "#9966CC", "1", "#999999", "#FF00FF"))+
                  theme(axis.text.x = element_text(size = 7))

## Plot density plot of samples remove 4x first run samples ----
unique(df$Sample)

DF_fil <- filter(df, Sample!="NaPPi_1_old", Sample!="EDTA_2_old", Sample!="Water_3_old", Sample!="NaPPi_control_7_old")


DF_fil <- droplevels(DF_fil)
unique(DF_fil$Sample)

#DF_fil$Sample <- factor(DF_fil$Sample, levels = c(
#                                         "EDTA_control_8", "Water_control_9", "NaPPi_1",
#                                          "EDTA_2", "Water_3", "NaPPi_control_7", "NaPPi_10",
#                                          "EDTA_11", "Water_12", "NaPPi_19", "EDTA_20",
#                                          "Water_21"))

DF_fil$Sample <- factor(DF_fil$Sample, levels = c("Water_control_9", "Water_3", "Water_12", "Water_21",
                                                  "NaPPi_control_7", "NaPPi_1", "NaPPi_10", "NaPPi_19",
                                                  "EDTA_control_8", "EDTA_2", "EDTA_11", "EDTA_20"))


unique(DF_fil$Category)

DF_fil$Category <- factor(DF_fil$Category, levels = c("Assigned_QC", "Anchored_QC", 
                                                      "Unassigned", "Unmatched", "Assigned_contams", "Anchored_contams"))

V_label <- c( 
  "Water control", "Water 1", "Water 2", "Water 3",
  "NaPPi control", "NaPPi 1", "NaPPi 2", "NaPPi 3",
  "EDTA control", "EDTA 1", "EDTA 2", "EDTA 3"
)

names(V_label) <- c("Water_control_9", "Water_3", "Water_12", "Water_21",
                    "NaPPi_control_7", "NaPPi_1", "NaPPi_10", "NaPPi_19",
                    "EDTA_control_8", "EDTA_2", "EDTA_11", "EDTA_20")


PlotDens <-   ggplot(DF_fil, aes(x=Retention.time.range, colour=Category)) + 
  #geom_histogram(aes(y=..count.., fill=Category), binwidth = 2.5, position = "dodge")+
  geom_density(show.legend = TRUE, linetype="solid", trim=TRUE)+
  facet_wrap(~Sample, scales = "free_y", ncol=4, labeller = labeller(Sample=V_label))+
  theme_bw()+
  #xlim(0,100)+
  #scale_x_continuous(breaks = seq(0, 100, by = 10), minor_breaks = seq(0, 100,2.5), expand = c(0, 0))+
  scale_x_continuous(limits=c(0,100), breaks = c(0,10,20,30,40,50,60,70,80,90,100), labels = c(0,10,20,30,40,50,60,70,80,90,100))+
  #geom_vline(xintercept = c(60,75,80,85,87.5,100), linetype='dashed', color='black')+
  xlab("Retention time (min)")+
  ylab("Normalised frequency \n")+
  #scale_color_brewer(palette="Set1")+
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c("#FF6633", "#66CC00", "#56B4E9", "#9966CC", "1", "#999999", "#FF00FF"))+
  scale_color_manual(values=c("#FF6633", "#66CC00", "#56B4E9", "#9966CC", "1", "#999999", "#FF00FF"))+
  theme(axis.text.x = element_text(size = 7), legend.text = element_text(size = 7),
        legend.position = "bottom")


### Plot 4x first run samples as run on MS instrument (only NCBI FDR 1%)  ----
unique(DF$Sample)

DF_test <- filter(DF, Sample=="NaPPi_1_old"| Sample=="EDTA_2_old"| Sample=="Water_3_old"| Sample=="NaPPi_control_7_old")


DF_test <- droplevels(DF_test)
unique(DF_test$Sample)

DF_test$Sample <- factor(DF_test$Sample, levels = c("NaPPi_1_old",
                                                   "EDTA_2_old",
                                                   "Water_3_old",
                                                   "NaPPi_control_7_old"
  ))

DF_test$Category <- as.character(DF_test$Category)
DF_test <- DF_test %>%
  mutate(., Category = with(., case_when(
    (Category=="Assigned_QC") ~ "Other",
    (Category=="Anchored_QC") ~ "Other",
    (Category=="Assigned_contams") ~ "Other",
    (Category=="Anchored_contams") ~ "Other",
    TRUE                      ~ Category
  )
  )
  )

unique(DF_test$Category)

DF_test$Category <- factor(DF_test$Category, levels = c("Unassigned", 
                                                      "Unmatched", 
                                                      "Other"))

V_label <- c("NaPPi 1 test",
             "EDTA 1 test",
             "Water 1 test",
             "NaPPi control test"
)

names(V_label) <- c("NaPPi_1_old",
                    "EDTA_2_old",
                    "Water_3_old",
                    "NaPPi_control_7_old")


PlotHistDens_test <-   ggplot(DF_test, aes(x=Retention.time.range, colour=Category)) + 
  geom_histogram(aes(y=..count.., fill=Category), binwidth = 2.5, position = "dodge")+
  geom_density(aes(y=2.5*..count..), show.legend = TRUE, linetype="solid", trim=TRUE)+
  facet_wrap(~Sample, scales = "fixed", ncol=4, labeller = labeller(Sample=V_label))+ 
  theme_bw()+
  #xlim(0,100)+
  #scale_x_continuous(breaks = seq(0, 100, by = 10), minor_breaks = seq(0, 100,2.5), expand = c(0, 0))+
  scale_x_continuous(limits=c(0,100), breaks = c(0,10,20,30,40,50,60,70,80,90,100), labels = c(0,10,20,30,40,50,60,70,80,90,100))+
  geom_vline(xintercept = c(60,75,80,85,87.5,100), linetype='dashed', color='black')+
  xlab("Retention time (min)")+
  ylab("Spectral counts \n")+
  #scale_color_brewer(palette="Set1")+
  theme(legend.title = element_blank())+
  scale_fill_manual(values=c("#56B4E9", "#9966CC", "1", "#999999", "#FF00FF"))+
  scale_color_manual(values=c("#56B4E9", "#9966CC", "1", "#999999", "#FF00FF"))

## Plot all samples for supplementary for revision 2 EJSS ----
DF_fil <- filter(df, Sample!="NaPPi_1_old", Sample!="EDTA_2_old", Sample!="Water_3_old", Sample!="NaPPi_control_7_old")

DF_fil <- droplevels(DF_fil)
unique(DF_fil$Sample)

DF_fil$Category <- as.character(DF_fil$Category)
DF_fil <- DF_fil %>%
  mutate(., Category = with(., case_when(
    (Category=="Assigned_QC") ~ "Other",
    (Category=="Anchored_QC") ~ "Other",
    (Category=="Assigned_contams") ~ "Other",
    (Category=="Anchored_contams") ~ "Other",
    TRUE                      ~ Category
  )
  )
  )

unique(DF_fil$Category)

DF_fil$Category <- factor(DF_fil$Category, levels = c("Unassigned", 
                                                      "Unmatched", 
                                                      "Other"))

DF_fil <- DF_fil %>% mutate(Interaction=paste(Sample, Category, sep="_"))
unique(DF_fil$Interaction)

DF_fil$Interaction <- factor(DF_fil$Interaction, levels = c(
                                                      "Water_control_9_Other", "Water_control_9_Unassigned", "Water_control_9_Unmatched",
                                                      "Water_3_Other", "Water_3_Unassigned", "Water_3_Unmatched",
                                                      "Water_12_Other", "Water_12_Unassigned", "Water_12_Unmatched",
                                                      "Water_21_Other", "Water_21_Unassigned", "Water_21_Unmatched",
                                                      "NaPPi_control_7_Other", "NaPPi_control_7_Unassigned", "NaPPi_control_7_Unmatched",
                                                      "NaPPi_1_Other", "NaPPi_1_Unassigned", "NaPPi_1_Unmatched",
                                                      "NaPPi_10_Other", "NaPPi_10_Unassigned", "NaPPi_10_Unmatched",
                                                      "NaPPi_19_Other", "NaPPi_19_Unassigned", "NaPPi_19_Unmatched",
                                                      "EDTA_control_8_Other", "EDTA_control_8_Unassigned", "EDTA_control_8_Unmatched",
                                                      "EDTA_2_Other", "EDTA_2_Unassigned", "EDTA_2_Unmatched",
                                                      "EDTA_11_Other", "EDTA_11_Unassigned", "EDTA_11_Unmatched",
                                                      "EDTA_20_Other", "EDTA_20_Unassigned", "EDTA_20_Unmatched"
                                                      ))



V_label <- c(
  "Water control Other", "Water control Unassigned", "Water control Unmatched",
  "Water 1 Other", "Water 1 Unassigned", "Water 1 Unmatched",
  "Water 2 Other", "Water 2 Unassigned", "Water 2 Unmatched",
  "Water 3 Other", "Water 3 Unassigned", "Water 3 Unmatched",
  "NaPPi control Other", "NaPPi control Unassigned", "NaPPi control Unmatched",
  "NaPPi 1 Other", "NaPPi 1 Unassigned", "NaPPi 1 Unmatched",
  "NaPPi 2 Other", "NaPPi 2 Unassigned", "NaPPi 2 Unmatched",
  "NaPPi 3 Other", "NaPPi 3 Unassigned", "NaPPi 3 Unmatched",
  "EDTA control Other", "EDTA control Unassigned", "EDTA control Unmatched",
  "EDTA 1 Other", "EDTA 1 Unassigned", "EDTA 1 Unmatched",
  "EDTA 2 Other", "EDTA 2 Unassigned", "EDTA 2 Unmatched",
  "EDTA 3 Other", "EDTA 3 Unassigned", "EDTA 3 Unmatched"
)


names(V_label) <- c(
  "Water_control_9_Other", "Water_control_9_Unassigned", "Water_control_9_Unmatched",
  "Water_3_Other", "Water_3_Unassigned", "Water_3_Unmatched",
  "Water_12_Other", "Water_12_Unassigned", "Water_12_Unmatched",
  "Water_21_Other", "Water_21_Unassigned", "Water_21_Unmatched",
  "NaPPi_control_7_Other", "NaPPi_control_7_Unassigned", "NaPPi_control_7_Unmatched",
  "NaPPi_1_Other", "NaPPi_1_Unassigned", "NaPPi_1_Unmatched",
  "NaPPi_10_Other", "NaPPi_10_Unassigned", "NaPPi_10_Unmatched",
  "NaPPi_19_Other", "NaPPi_19_Unassigned", "NaPPi_19_Unmatched",
  "EDTA_control_8_Other", "EDTA_control_8_Unassigned", "EDTA_control_8_Unmatched",
  "EDTA_2_Other", "EDTA_2_Unassigned", "EDTA_2_Unmatched",
  "EDTA_11_Other", "EDTA_11_Unassigned", "EDTA_11_Unmatched",
  "EDTA_20_Other", "EDTA_20_Unassigned", "EDTA_20_Unmatched")


PlotHistDens_all <-   ggplot(DF_fil, aes(x=Retention.time.range, colour=Category)) + 
  geom_histogram(aes(y=..count.., fill=Category), binwidth = 2.5, position = "dodge")+
  geom_density(aes(y=2.5*..count..), show.legend = TRUE, linetype="solid", trim=TRUE)+
  facet_wrap(~Interaction, scales = "free_y", ncol = 3, labeller = labeller(Interaction=V_label))+ 
  theme_bw()+
  #xlim(0,100)+
  #scale_x_continuous(breaks = seq(0, 100, by = 10), minor_breaks = seq(0, 100,2.5), expand = c(0, 0))+
  scale_x_continuous(limits=c(0,100), breaks = c(0,10,20,30,40,50,60,70,80,90,100), labels = c(0,10,20,30,40,50,60,70,80,90,100))+
  #geom_vline(xintercept = c(60,75,80,85,87.5,100), linetype='dashed', color='black')+
  xlab("Retention time (min)")+
  ylab("Spectral counts \n")+
  #scale_color_brewer(palette="Set1")+
  theme(legend.position = "none") +
  #theme(legend.title = element_blank())+
  scale_fill_manual(values=c("#56B4E9", "#9966CC", "1", "#999999", "#FF00FF"))+
  scale_color_manual(values=c("#56B4E9", "#9966CC", "1", "#999999", "#FF00FF"))

PlotHistDens_all

### Plot water_3 (Water rep 1) only ----
DF_fil <- filter(DF, Sample=="Water_3")

DF_fil <- droplevels(DF_fil) 

DF_fil$Category <- as.character(DF_fil$Category)
DF_fil <- DF_fil %>%
  mutate(., Category = with(., case_when(
    (Category=="Assigned_QC") ~ "Other",
    (Category=="Anchored_QC") ~ "Other",
    (Category=="Assigned_contams") ~ "Other",
    (Category=="Anchored_contams") ~ "Other",
    TRUE                      ~ Category
  )
  )
  )

unique(DF_fil$Category)

DF_fil$Category <- factor(DF_fil$Category, levels = c("Unassigned", 
                                                      "Unmatched", 
                                                      "Other"))

DF_fil <- DF_fil %>% mutate(Interaction=paste(Sample, Category, sep="_"))
unique(DF_fil$Interaction)

DF_fil$Interaction <- factor(DF_fil$Interaction, levels = c(
  
  "Water_3_Other", "Water_3_Unassigned", "Water_3_Unmatched"
  ))



V_label <- c(
  
  "Water 1 Other", "Water 1 Unassigned", "Water 1 Unmatched"
  )


names(V_label) <- c(
                    "Water_3_Other", "Water_3_Unassigned", "Water_3_Unmatched"
                    )


PlotHistDens <-   ggplot(DF_fil, aes(x=Retention.time.range, colour=Category)) + 
  geom_histogram(aes(y=..count.., fill=Category), binwidth = 2.5, position = "dodge")+
  geom_density(aes(y=2.5*..count..), show.legend = TRUE, linetype="solid", trim=TRUE)+
  facet_wrap(~Interaction, scales = "free_y", ncol = 1, labeller = labeller(Interaction=V_label))+ 
  theme_bw()+
  #xlim(0,100)+
  #scale_x_continuous(breaks = seq(0, 100, by = 10), minor_breaks = seq(0, 100,2.5), expand = c(0, 0))+
  scale_x_continuous(limits=c(0,100), breaks = c(0,10,20,30,40,50,60,70,80,90,100), labels = c(0,10,20,30,40,50,60,70,80,90,100))+
  geom_vline(xintercept = c(60,75,80,85,87.5,100), linetype='dashed', color='black')+
  xlab("Retention time (min)")+
  ylab("Spectral counts \n")+
  #scale_color_brewer(palette="Set1")+
  theme(legend.position = "none") +
  #theme(legend.title = element_blank())+
  scale_fill_manual(values=c("#56B4E9", "#9966CC", "1", "#999999", "#FF00FF"))+
  scale_color_manual(values=c("#56B4E9", "#9966CC", "1", "#999999", "#FF00FF"))

PlotHistDens

# MS2 fragments normalised to parent ion mass ----
df <- DF[,-c(6,10)]

max(df$charge)
boxplot(df$charge)
hist(df$charge)

df$parent_mz <- df$moverz*df$charge
max(df$parent_mz)
min(df$parent_mz)
boxplot(df$parent_mz)
hist(log(df$parent_mz))

df$fragment_norm <- df$NumVals/df$parent_mz
max(df$fragment_norm)
min(df$fragment_norm)
boxplot(df$fragment_norm)
hist(df$fragment_norm)


df$Category <- as.factor(df$Category)

df$Category <- factor(df$Category, levels = c("Assigned_QC", 
                                              "Anchored_QC",
                                              "Unassigned", "Unmatched", 
                                              "Assigned_contams", "Anchored_contams"
))  

df$Retention.time.range <- as.numeric(df$Retention.time.range)
df$fragment_norm <- as.numeric(df$fragment_norm)

summary(df)
summary(df$Category)


DF_fil <- filter(df, Sample!="NaPPi_1_old", Sample!="EDTA_2_old", Sample!="Water_3_old", Sample!="NaPPi_control_7_old")


DF_fil <- droplevels(DF_fil)
unique(DF_fil$Sample)

DF_fil$Sample <- factor(DF_fil$Sample, levels = c("Water_control_9", "Water_3", "Water_12", "Water_21",
                                                  "NaPPi_control_7", "NaPPi_1", "NaPPi_10", "NaPPi_19",
                                                  "EDTA_control_8", "EDTA_2", "EDTA_11", "EDTA_20"))

V_label <- c( 
  "Water control", "Water 1", "Water 2", "Water 3",
  "NaPPi control", "NaPPi 1", "NaPPi 2", "NaPPi 3",
  "EDTA control", "EDTA 1", "EDTA 2", "EDTA 3"
)

names(V_label) <- c("Water_control_9", "Water_3", "Water_12", "Water_21",
                    "NaPPi_control_7", "NaPPi_1", "NaPPi_10", "NaPPi_19",
                    "EDTA_control_8", "EDTA_2", "EDTA_11", "EDTA_20")


PlotNumVals <-   ggplot(DF_fil%>%filter(fragment_norm<=0.2), aes(x=Retention.time.range, y= fragment_norm, colour=Category)) + 
                  geom_point(size=0.0002,alpha=0.5)+ 
                  geom_smooth(method = "auto", aes(fill=Category))+
                  facet_wrap(~Sample, scales = "fixed", ncol=4, labeller = labeller(Sample=V_label))+
                  scale_x_continuous(limits=c(0,100), breaks = c(0,10,20,30,40,50,60,70,80,90,100), labels = c(0,10,20,30,40,50,60,70,80,90,100))+
                  #geom_vline(xintercept = c(60,75,80,85,87.5,100), linetype='dashed', color='black')+
                  theme_bw()+
                  xlab("Retention time (min)")+
                  ylab("Fragment ions normalised to parent ion mass \n")+
                  #scale_color_brewer(palette="Set1")+
                  theme(legend.title = element_blank(), legend.position = "bottom")+
                  scale_fill_manual(values=c("#FF6633", "#66CC00", "#56B4E9", "#9966CC", "1", "#999999", "#FF00FF"))+
                  scale_color_manual(values=c("#FF6633", "#66CC00", "#56B4E9", "#9966CC", "1", "#999999", "#FF00FF"))

PlotNumVals

## Pick specific for larger plot ----
unique(df$Sample)

DF_fil <- filter(df, Sample=="Water_3")

DF_fil <- droplevels(DF_fil)
unique(DF_fil$Sample)

PlotNumVals <-   ggplot(DF_fil%>%filter(fragment_norm<=0.2), aes(x=Retention.time.range, y= fragment_norm, colour=Category)) + 
                  geom_point(size=0.0002,alpha=0.5)+ 
                  geom_smooth(method = "auto", aes(fill=Category))+
                  #facet_wrap(~Sample, scales = "fixed", ncol=4, labeller = labeller(Sample=V_label))+
                  scale_x_continuous(limits=c(0,100), breaks = c(0,10,20,30,40,50,60,70,80,90,100), labels = c(0,10,20,30,40,50,60,70,80,90,100))+
                  #geom_vline(xintercept = c(60,75,80,85,87.5,100), linetype='dashed', color='black')+
                  theme_bw()+
                  xlab("Retention time (min)")+
                  ylab("Fragment ions normalised to parent ion mass \n")+
                  #scale_color_brewer(palette="Set1")+
                  theme(legend.title = element_blank(), legend.position = "bottom")+
                  scale_fill_manual(values=c("#FF6633", "#66CC00", "#56B4E9", "#9966CC", "1", "#999999", "#FF00FF"))+
                  scale_color_manual(values=c("#FF6633", "#66CC00", "#56B4E9", "#9966CC", "1", "#999999", "#FF00FF"))

  
  scale_fill_manual(values=c("#56B4E9", "#9966CC", "1"))+
  scale_color_manual(values=c("#56B4E9", "#9966CC", "1"))
  
  
  scale_fill_manual(values=c("#56B4E9", "#9966CC", "1", "#999999", "#FF00FF"))+
  scale_color_manual(values=c("#56B4E9", "#9966CC", "1", "#999999", "#FF00FF"))

library(ggExtra)
ggMarginal(PlotNumVals, type = "density", margins = "both", groupColour = TRUE, groupFill = FALSE, xparams = list(trim=TRUE), yparams = list(trim=TRUE))


# Extra if no regression due too few points
DF_fil$Category <- as.character(DF_fil$Category)
DF_fil <- DF_fil %>%
  mutate(., Category = with(., case_when(
    (Category=="Assigned_QC") ~ "Other",
    (Category=="Anchored_QC") ~ "Other",
    (Category=="Assigned_contams") ~ "Other",
    (Category=="Anchored_contams") ~ "Other",
    TRUE                      ~ Category
  )
  )
  )

unique(DF_fil$Category)

DF_fil$Category <- factor(DF_fil$Category, levels = c("Unassigned", 
                                                      "Unmatched", 
                                                      "Other"))

#############################################
library(data.table)
library(dplyr)

#input
myfile1="../../0_data/R2_result/result_R2_20PRS_to_ObsMetS_ATTRaCT.txt"
myfile2="../../0_data/R2_result/result_R2_20PRS_to_ObsMetS_GUSTOchildren.txt"
myfile3="../../0_data/R2_result/result_R2_20PRS_to_ObsMetS_UKBB.txt"
myfile4="../../0_data/R2_result/result_R2_20PRS_to_otherTraits_ATTRaCT.txt"
myfile5="../../0_data/R2_result/result_R2_20PRS_to_otherTraits_GUSTOchildren.txt"
myfile6="../../0_data/R2_result/result_R2_20PRS_to_otherTraits_UKBB.txt"

#output
myoutput="R2_result.txt"

#data <- fread("your_file.csv", colClasses = c(R2 = "character"))
mytmp1 <- fread(myfile1); mytmp1[1:3, ]
mytmp2 <- fread(myfile2); mytmp2[1:3, ]
mytmp3 <- fread(myfile3); mytmp3[1:3, ]

print(colnames(mytmp1))
print(colnames(mytmp2))
print(colnames(mytmp3))

mydata1 <- rbind(mytmp1, mytmp2, mytmp3) 
mydata1 <- mydata1 %>% dplyr::filter(predictor == "PRS_GenMetS")
mydata1[1:3, ]

#GenMetS_to_relatedTraits
################2. combination to output2
attract_trait<- c("BMI", "CHOLESTROL", "DBP", "FastingGlucose", "HDL", "LDL", "SBP", "TG", "ObsMetS")
mytmp1 <- fread(myfile4)  
idx<-which(mytmp1$target=="FastingGlu"); mytmp1$target[idx]<- "FastingGlucose"
mytmp1 <- mytmp1 %>% dplyr::filter(target %in% attract_trait)

#gusto children
children_trait <- c("Abdominal Circumference", "BMI", "Fatty Liver Index", "HDL",
                    "HOMA_IR", "Insulin", "LDL", "logHSCRP", "ObsMetS", "TG")

mytmp2 <- fread(myfile5)
mytmp2$cohort <- "GUSTOchildren6yr"
mytmp2$target <- gsub("|yr6|child|mmol_L|mUL", "", mytmp2$target)
mytmp2$target <- gsub("_", "", mytmp2$target)
idx<-which(mytmp2$target=="FLI"); mytmp2$target[idx]<- "Fatty Liver Index"
idx<-which(mytmp2$target=="cabdominal"); mytmp2$target[idx]<- "Abdominal Circumference"
idx<-which(mytmp2$target=="HOMAIR"); mytmp2$target[idx]<- "HOMA_IR"
idx<-which(mytmp2$target=="INSmUL"); mytmp2$target[idx]<- "Insulin"
mytmp2<- mytmp2 %>% dplyr::filter(target %in% children_trait)
mytmp2 <- mytmp2 %>% dplyr::select(-stage)

#ukb
ukb_trait <- c("BMI", "DBP", "HbA1c", "HDL", "ObsMetS","SBP", "TG","WC")
mytmp3 <- fread(myfile6)  
mytmp3 <- mytmp3 %>% dplyr::filter(target %in% ukb_trait)
unique(mytmp3$target)

mydata2<- rbind(mytmp1, mytmp2, mytmp3)

mydata <- rbind(mydata2, mydata1) %>% data.frame()
colnames(mydata)[which(colnames(mydata)== "X95.CI")]<- "95%CI"
colnames(mydata)[which(colnames(mydata)== "p.value")]<- "p-value"

idx<-which(mydata$target == "WC"); mydata$target[idx]<- "Waist Circumference"

mydata[duplicated(mydata),]$predictor
# remove duplicate row
mydata <- mydata[duplicated(mydata) == "FALSE",]
dim(mydata);  

mydata$predictor <- "GenMetS"

write.table(mydata, myoutput, row.names=F, sep="\t", quote=F)


mydata %>% dplyr::filter(cohort == "GUSTOchildren6yr") %>% dplyr::select(target) %>% unique()
mydata %>% dplyr::filter(cohort == "ATTRaCT") %>% dplyr::select(target) %>% unique()
mydata %>% dplyr::filter(cohort %in% c( "UKBB Asian", "UKBB European")) %>% dplyr::select(target) %>% unique()


cat("output = ", myoutput , "\n")

##############

rm(list=ls())
library(data.table)
library(dplyr)
library(forestplot)

#setwd(mywd)
data_all <- fread("../../0_data/R2_result/R2_result.txt")
dim(data_all); head(data_all)

colnames(data_all)[10]<- "p.value"
# check if any duplicate row
duplicated(data_all)
#data_all[duplicated(data_all),]

# remove duplicate row
data <- data_all[duplicated(data_all) == "FALSE",]
dim(data); head(data)

# filter for only predictor == "GenMetS"
data <- data[data$predictor == "GenMetS",]
dim(data); head(data)

# change Abdominal circumference to AC*
# and waist circumference to WC*
data$target[data$target == "Abdominal Circumference"] <- "AC"
data$target[data$target == "Waist Circumference"] <- "WC"

# for UKBB European
dat1 <- data[cohort=="UKBB European"]
dat1 <- dat1[order(target, sex),]

# Find the row with "ObsMetS"
obsMetS_row <- dat1[target == "ObsMetS", ]

# Remove "ObsMetS" from dat1 to avoid duplication
dat1 <- dat1[target != "ObsMetS", ]

# sort by R2 of girls
sorted_girls <- dat1[sex == "women"][order(R2, decreasing = TRUE)]

# Extract the sorted target names
sorted_targets <- sorted_girls$target

# Reorder the original data frame based on the sorted target names
dat1 <- dat1[target %in% sorted_targets]
dat1 <- dat1[order(match(target, sorted_targets), sex),]

# Add "ObsMetS" row on top
dat1 <- rbind(obsMetS_row, dat1)


# for UKBB Asian
dat2 <- data[cohort=="UKBB Asian"]
dat2 <- dat2[order(target,sex),]

# Find the row with "ObsMetS"
obsMetS_row <- dat2[target == "ObsMetS", ]

# Remove "ObsMetS" from dat1 to avoid duplication
dat2 <- dat2[target != "ObsMetS", ]

# sort by R2 of girls
sorted_girls <- dat2[sex == "women"][order(R2, decreasing = TRUE)]

# Extract the sorted target names
sorted_targets <- sorted_girls$target

# Reorder the original data frame based on the sorted target names
dat2 <- dat2[target %in% sorted_targets]
dat2 <- dat2[order(match(target, sorted_targets), sex),]

# Add "ObsMetS" row on top
dat2 <- rbind(obsMetS_row, dat2)
dat1 <- rbind(dat1, dat2)

# ATTRaCT
dat3 <- data[cohort=="ATTRaCT"]
dat3 <- dat3[order(target,sex),]

# find and replace some words
dat3$target <- gsub("CHOLESTROL", "Cholestrol", dat3$target)
dat3$target <- gsub("FastingGlucose", "Fasting Glucose", dat3$target)

# Find the row with "ObsMetS"
obsMetS_row <- dat3[target == "ObsMetS", ]

# Remove "ObsMetS" from dat1 to avoid duplication
dat3 <- dat3[target != "ObsMetS", ]

# sort by R2 of girls
sorted_girls <- dat3[sex == "women"][order(R2, decreasing = TRUE)]

# Extract the sorted target names
sorted_targets <- sorted_girls$target

# Reorder the original data frame based on the sorted target names
dat3 <- dat3[target %in% sorted_targets]
dat3 <- dat3[order(match(target, sorted_targets), sex),]

# Add "ObsMetS" row on top
dat3 <- rbind(obsMetS_row, dat3)

dat1 <- rbind(dat1, dat3)
dim(dat1); head(dat1); tail(dat1)

# GUSTO
dat4 <- data[(cohort=="GUSTOchildren" | cohort == "GUSTOchildren6yr")]
dat4 <- dat4[order(target,sex),]

dat4$cohort <- "GUSTO"

# remove fatty liver index
dat4 <- dat4[-which(dat4$target == "Fatty Liver Index"),]

# Find the row with "ObsMetS"
obsMetS_row <- dat4[target == "ObsMetS", ]

# Remove "ObsMetS" from dat1 to avoid duplication
dat4 <- dat4[target != "ObsMetS", ]

# sort by R2 of girls
sorted_girls <- dat4[sex == "girl"][order(R2, decreasing = TRUE)]

# Extract the sorted target names
sorted_targets <- sorted_girls$target

# Reorder the original data frame based on the sorted target names
dat4 <- dat4[target %in% sorted_targets]
dat4 <- dat4[order(match(target, sorted_targets), sex),]

# Add "ObsMetS" row on top
dat4 <- rbind(obsMetS_row, dat4)

dat1 <- rbind(dat1, dat4)
dat1$sex[dat1$sex == "boy"] <- "men"
dat1$sex[dat1$sex == "girl"] <- "women"
dim(dat1); head(dat1); tail(dat1)
table(dat1$cohort)

#dat1$R2 <- format(signif(dat1$R2,4), scientific = T, digits = 3)
#dat1$low <- format(signif(dat1$low,4), scientific = T, digits = 3)
#dat1$high <- format(signif(dat1$high,4), scientific = T, digits = 3)
dat1$p.value <- format(signif(dat1$p.value,4), scientific = T, digits = 3)
#dat1$CI2 <- paste0("(", dat1$low, ",", dat1$high, ")")

# Plot
dat1$target[dat1$target == "ObsMetS"] <- "MetS"
write.table(dat1, "data_table.txt", row.names=F, sep="\t", quote=F) 



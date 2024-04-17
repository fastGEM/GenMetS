#############################################
library(data.table)
library(dplyr)

#input
myfile1="R2_result/result_R2_20PRS_to_ObsMetS_ATTRaCT.txt"
myfile2="R2_result/result_R2_20PRS_to_ObsMetS_GUSTOchildren.txt"
myfile3="R2_result/result_R2_20PRS_to_ObsMetS_UKBB.txt"
myfile4="R2_result/result_R2_20PRS_to_otherTraits_ATTRaCT.txt"
myfile5="R2_result/result_R2_20PRS_to_otherTraits_GUSTOchildren.txt"
myfile6="R2_result/result_R2_20PRS_to_otherTraits_UKBB.txt"

#output
myoutput="R2_result.txt"

###############
mytmp1<-fread(myfile1); mytmp2=fread(myfile2); mytmp3=fread(myfile3)
idx<-which(mytmp1$predictor=="MetS_PRS")
if(length(idx)>0 ) mytmp1$predictor[idx]= "GenMetS"

idx<-which(mytmp2$predictor=="predicted_mets") ; 
if(length(idx)>0 )  mytmp2$predictor[idx]= "GenMetS"

colnames(mytmp2)<- colnames(mytmp1)
colnames(mytmp3)<- colnames(mytmp1)
mydata1 <- rbind(mytmp1, mytmp2, mytmp3)


#GenMetS_to_relatedTraits
################2. combination to output2
myref <- fread("old_coef.txt")[, c(1:10)]

#attract
attract_trait <- myref %>% dplyr::filter(cohort=="ATTRaCT") %>% dplyr::filter(!is.na(Trait)) 
attract_trait <- unique(attract_trait$Trait)
print(attract_trait)
mytmp1 <- fread(myfile4) %>% data.frame() 
idx<-which(mytmp1$target=="FastingGlu"); mytmp1$target[idx]<- "FastingGlucose"
mytmp1 <- mytmp1 %>% dplyr::filter(target %in% attract_trait)
idx<-which(mytmp1$predictor=="MetS_PRS"); mytmp1$predictor[idx]<- "GenMetS"

#gusto children
children_trait <- myref %>% dplyr::filter(cohort=="GUSTOchildren") %>% dplyr::filter(!is.na(Trait)) 
children_trait <- unique(children_trait$Trait)
print(children_trait)
mytmp2 <- fread(myfile5) %>% data.frame() 
mytmp2$target <- gsub("|yr6|child|mmol_L|mUL", "", mytmp2$target)
mytmp2$target <- gsub("_", "", mytmp2$target)
unique(mytmp2$target)
idx<-which(mytmp2$target=="FLI"); mytmp2$target[idx]<- "Fatty Liver Index"
idx<-which(mytmp2$target=="cabdominal"); mytmp2$target[idx]<- "Abdominal Circumference"
idx<-which(mytmp2$target=="HOMAIR"); mytmp2$target[idx]<- "HOMA_IR"
idx<-which(mytmp2$target=="INSmUL"); mytmp2$target[idx]<- "Insulin"
mytmp2<- mytmp2 %>% dplyr::filter(target %in% children_trait)
mytmp2 <- mytmp2 %>% dplyr::select(-stage)
idx<-which(mytmp2$predictor == "predicted_mets"); 
if(length(idx) > 0) mytmp2$predictor[idx]="GenMetS"

#ukb
#gusto children
ukb_trait <- myref %>% dplyr::filter(cohort=="UKB") 
ukb_trait <- unique(ukb_trait$Trait)
print(ukb_trait)
ukb_trait[3]<- "HbA1c"
mytmp3 <- fread(myfile6) %>% data.frame() 
unique(mytmp3$target)
idx<-which(mytmp3$target=="Glucose"); mytmp3$target[idx]<- "FastingGlucose"
mytmp3 <- mytmp3 %>% dplyr::filter(target %in% ukb_trait)
unique(mytmp3$target)

mydata2<- rbind(mytmp1, mytmp2, mytmp3)

mydata <- rbind(mydata2, mydata1) %>% data.frame()
colnames(mydata)[which(colnames(mydata)== "X95.CI")]<- "95%CI"
colnames(mydata)[which(colnames(mydata)== "p.value")]<- "p-value"

idx<-which(mydata$target == "WC"); mydata$target[idx]<- "Waist Circumference"

mydata <- mydata %>% dplyr::filter(sex !="all")
# check if any duplicate row

mydata[duplicated(mydata),]$predictor
# remove duplicate row
mydata <- mydata[duplicated(mydata) == "FALSE",]
dim(mydata);  
write.table(mydata, myoutput, row.names=F, sep="\t", quote=F)


cat("output = ", myoutput , "\n")
#mydata1 <- mydata %>% dplyr::filter(predictor=="GenMetS" &  sex != "all")
#write.table(mydata1, myoutput, row.names=F, sep="\t", quote=F)
#myoutput="Coefficient_determination_PRS_GenMetS_to_continoustraits_CI_modified.txt"
#mydata1 <- mydata %>% dplyr::filter(predictor=="PRS_GenMetS" &  sex != "all")
#write.table(mydata1, myoutput, row.names=F, sep="\t", quote=F)

#GenMetS, Lind2019, Walree2022 in explaining
#ObsMetS in ATTRaCT, GUSTO children, UKB Asian
#mydata <- fread("r2_result_all.txt")
#mydata1 <- mydata %>% dplyr::filter(target=="ObsMetS") %>%
#  dplyr::filter(predictor %in% c("GenMetS", "PRS_GenMetS", "PRS_Lind2019", "PRS_Walree2022")) %>%
# dplyr::filter(cohort != "UKBB European")

#mydata1$low <-  format(round(mydata1$low, 2), nsmall = 2)
#mydata1$high <- format(round(mydata1$high, 2), nsmall = 2)
#mydata1$X95.CI <- paste0("(", mydata1$low, ",", mydata1$high, ")")
#write.table(mydata1, "table_GenMetS_Lind2019_Walree2022_vs_ObsMetS.txt", row.names=F, sep="\t", quote=F)


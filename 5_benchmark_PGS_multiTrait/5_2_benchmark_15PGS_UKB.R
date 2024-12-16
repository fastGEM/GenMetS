######################
######################

rm(list=ls())
library(tidyverse)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(ggplot2)
library(scales) #for scientific number
library(epitools) #for odds ratio
library(questionr) #odds.ratio
library(gridExtra)
library(forcats)
#####################
set.seed(1234)

####################

OddsRatio_for_diseases <- function(myinputdata, mycohort){
#7 diseases x 2 (case-control, early-late) x 2 (all and matched) x 2 (PredMetS and ObsMetS)
#matched case and control
myresult <- matrix(rep(NA, 20*7*11), ncol=11) 
colnames(myresult)<- c("disease", "predictor" , "samples", "test", "case", "control", "OddsRatio", "low", "high", "p-value", "cohort")
library(MatchIt)
library(oddsratio)
K=1
for( I in c(1:7)){
	diseases <- c("t2dm", "hypertension", "HF", "CAD", "stroke", "NAFLD", "MI")
  	disease=diseases[I]; 
  	mydata <- myinputdata 
  	mydata <- mydata %>% dplyr::filter(!is.na(!!sym(disease)))

   	#match data by Age and Education 
   	fit <-  as.formula(paste0("factor(", disease, ") ~ Age + Education "))
   	matched <- matchit(fit, mydata, method="nearest", ratio=1,  na.rm=TRUE)
   	matched_data <- get_matches(matched); dim(matched_data)
   	
   	 	
   	#########matched data, 20PRS, case-control
   	for(prs in pgslist$PGS) {
     		fit <- as.formula(paste0("factor(", disease, ") ~ ", prs))
     		gfit <- glm(fit, data = matched_data, family = "binomial")
     		pvalue <- scientific(summary(gfit)$coef[2, 4], digits = 3)
     		t <- or_glm(data=matched_data, model=gfit, incr=setNames(list(1), prs))
     		print(t)
     		tmp <- matched_data %>% dplyr::select(!!sym(disease))
     		N <- unlist(table(tmp)); rm(tmp)
           	myresult[K, ] <- c(disease, prs, "matched", "case-control", 
				      N[2], N[1], t$oddsratio, t$`ci_low (2.5)`, t$`ci_high (97.5)`, pvalue, mycohort)
     		K <- K + 1
   	} #end for 20 PGSs
 } #for 7 diseases 

return(myresult)
} #close function

##########################################################MAIN
cat("reading data\n")

mydata <- fread("../../0_data/UKB_Asianwomenmen_s8792_mets_diseases_others_updated.txt") %>% data.frame()
##################GenMetS <- GenMetS + 19PRS = 20PRS
mydata$ID <-as.character(mydata$ID)
myPRS19 <- fread("../../0_data/studysample_20PGS.txt")

mydata2 <- merge(mydata, myPRS19, by.x="ID", by.y="sampleID")
mydata <- mydata2; rm(mydata2)
################################

cat("ensure health control\n")
tmp <- mydata %>% dplyr::filter(t2dm==0 & HF ==0 & hypertension==0 & stroke==0 & CAD==0 & MI==0) 
dim(tmp)
healthcontrol <- which(mydata$ID %in% tmp$ID)
##in R: mydata[[col]] equals mydata$col
set.seed(1234) 
diseases <- c("t2dm", "hypertension", "HF", "CAD", "stroke", "NAFLD", "MI")
diseasenames <-c("T2DM", "Hypertension", "Heart Failure", "CAD", "Stroke", "NAFLD", 
                 "Myocardial infarction")
                  
for (col in diseases) {
  idx <- which(mydata[[col]] == 1)
  mydata[[col]][idx] <- 1
  mydata[[col]][-idx] <- NA
  mydata[[col]][healthcontrol] <- 0
  print(table(mydata[[col]], useNA = "ifany"))
}

cat("prepare for or_glm\n")
mydata$t2dm <- factor(mydata$t2dm, levels=c(0, 1), labels=c("control", "case"))
mydata$HF <- factor(mydata$HF, levels=c(0, 1), labels=c("control", "case"))
mydata$hypertension <- factor(mydata$hypertension, levels=c(0, 1), labels=c("control", "case"))
mydata$NAFLD <- factor(mydata$NAFLD, levels=c(0, 1), labels=c("control", "case"))
mydata$stroke <- factor(mydata$stroke, levels=c(0, 1), labels=c("control", "case"))
mydata$CAD <- factor(mydata$CAD, levels=c(0, 1), labels=c("control", "case"))
mydata$MI <- factor(mydata$MI, levels=c(0,1), labels=c("control", "case"))

pgslist<-fread("PGSlist.txt")

cat("UKB ASN women\n")
myUKB_ASN_women <-  mydata %>% dplyr::filter(sex=="Female")
dim(myUKB_ASN_women)
colidx <- which(colnames(myUKB_ASN_women) %in% pgslist$PGS) 
Nrow <- length(colidx)
cat("we are studying ", Nrow, " PGSs\n")
myUKB_ASN_women[, c(colidx) ] <- lapply(myUKB_ASN_women[, c(colidx) ], function(x) c(scale(x)))
mycohort1 <- "UKB ASN women"
myresult1 <- OddsRatio_for_diseases(myUKB_ASN_women, mycohort1) 

cat("UKB ASN men\n")
myUKB_ASN_men  <- mydata %>% dplyr::filter(sex=="Male") ################
dim(myUKB_ASN_men)
colidx <- which(colnames(myUKB_ASN_men) %in% pgslist$PGS) 
Nrow <- length(colidx)
cat("we are studying ", Nrow, " PGSs\n")
myUKB_ASN_men[, c(colidx) ] <- lapply(myUKB_ASN_men[, c(colidx) ], function(x) c(scale(x)))
mycohort2 <- "UKB ASN men"
myresult2 <- OddsRatio_for_diseases(myUKB_ASN_men, mycohort2) 

myresult_all <- rbind(myresult1, myresult2)
cat("result table = OddsRatio_diseases_UKB.txt \n")
write.table(myresult_all, "../../0_data/OddsRatio_result/OddsRatio_diseases_UKB.txt", row.names=F, sep="\t", quote=F)


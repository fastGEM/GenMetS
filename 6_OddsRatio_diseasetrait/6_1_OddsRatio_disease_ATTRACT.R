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
#####################

OddsRatio_for_diseases <- function(myinputdata, mycohort){
#9 diseases x 2 (case-control, early-late) x 2 (all and matched) x 2 (PredMetS and ObsMetS)
#matched case and control
diseases <- c("T2DM", "STROKE", "PVD", "HF", "CAD", "MI", "HTN", "COPD", "DEPRESSION")
myresult <- matrix(rep(NA, 20*9*11), ncol=11) 
colnames(myresult)<- c("disease", "predictor" , "samples", "test", "case", "control", "OddsRatio", "low", "high", "p-value", "cohort")
library(MatchIt)
library(oddsratio)
K=1
for( I in c(1:9)){
  	disease=diseases[I]; 
  	mydata <- myinputdata 
  	mydata <- mydata %>% dplyr::filter(!is.na(!!sym(disease)))

   	#match data by Age and Education 
   	fit <-  as.formula(paste0("factor(", disease, ") ~ AGE "))
   	matched <- matchit(fit, mydata, method="nearest", ratio=1,  na.rm=TRUE)
   	matched_data <- get_matches(matched); dim(matched_data)
   	
   	#########matched data, 20PRS, case-control
   	for(prs in pgs20$PGS20) {
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



#####################################################
cat("reading data\n")
set.seed(1234)

myraw<- fread("study_data_ATTRACT.csv") %>% data.frame()
idx<-which(myraw$GENDER==1);myraw$GENDER[idx]<- "Male"
idx<-which(myraw$GENDER==2); myraw$GENDER[idx]<- "Female"
table(myraw$GENDER)

myPRS20 <- fread("../1_make_20PGS/studysample_20PRS.txt")
mydata2 <- merge(myraw, myPRS20, by.x="LID", by.y="sampleID")
myraw <- mydata2; rm(mydata2)



mydata<- myraw
diseases <- c("T2DM", "STROKE", "PVD", "HF", "CAD", "MI", "HTN", "COPD", "DEPRESSION")
cat("prepare for or_glm\n")
mydata$T2DM <- factor(mydata$T2DM, levels=c(0, 1), labels=c("control", "case"))
mydata$STROKE <- factor(mydata$STROKE, levels=c(0, 1), labels=c("control", "case"))
mydata$PVD <- factor(mydata$PVD, levels=c(0, 1), labels=c("control", "case"))
mydata$CAD <- factor(mydata$CAD, levels=c(0, 1), labels=c("control", "case"))
mydata$MI <- factor(mydata$MI, levels=c(0, 1), labels=c("control", "case"))
mydata$HTN <- factor(mydata$HTN, levels=c(0, 1), labels=c("control", "case"))
mydata$COPD <- factor(mydata$COPD, levels=c(0,1), labels=c("control", "case"))
mydata$DEPRESSION <- factor(mydata$DEPRESSION, levels=c(0,1), labels=c("control", "case"))

pgs20<-fread("../PGS20.txt")

myraw <- mydata; #backup myraw
cat("OddsRatio for ATTRACT women\n")
myATTRACT_women <- myraw %>% dplyr::filter(GENDER=="Female")
dim(myATTRACT_women)
colidx <- which(colnames(myATTRACT_women) %in% pgs20$PGS20)
Nrow <- length(colidx)
cat("we are studying ", Nrow, " PGSs\n")
colnames(myATTRACT_women)[colidx]
myATTRACT_women[, c(colidx) ] <- lapply(myATTRACT_women[, c(colidx) ], function(x) c(scale(x)))
mycohort1 = "ATTRaCT women"
myresult1 <- OddsRatio_for_diseases(myATTRACT_women, mycohort1)

cat("OddsRatio for ATTRACT men\n")
myATTRACT_men <- myraw %>% dplyr::filter(GENDER=="Male")
dim(myATTRACT_men)
colidx <- which(colnames(myATTRACT_men) %in% pgs20$PGS20)
Nrow <- length(colidx)
cat("we are studying ", Nrow, " PGSs\n")
colnames(myATTRACT_men)[colidx]
myATTRACT_men[, c(colidx) ] <- lapply(myATTRACT_men[, c(colidx) ], function(x) c(scale(x)))
mycohort2 = "ATTRaCT men"
myresult2 <- OddsRatio_for_diseases(myATTRACT_men, mycohort2)


myresult_all <- rbind(myresult1, myresult2)
write.table(myresult_all, "OddsRatio_diseases_ATTRaCT.txt", row.names=F, sep="\t", quote=F)




 




 



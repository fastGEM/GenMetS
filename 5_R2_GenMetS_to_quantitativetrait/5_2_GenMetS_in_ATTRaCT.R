########################

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
#data
#####################

myraw<- fread("attract_mets_data.txt") %>% data.frame()
myraw$Ethnicity <- factor(myraw$Ethnicity)

myPRS20 <- fread("../1_make_20PGS/studysample_20PRS.txt")
mydata2 <- merge(myraw, myPRS20, by.x="IID", by.y="sampleID")
myraw <- mydata2; rm(mydata2)

myfemale <- myraw %>% filter(GENDER=="Female"); dim(myfemale)
mymale <- myraw %>% filter(GENDER=="Male"); dim(mymale)

###################################ATTRACT women
#matrix(rep(NA, 12*6), ncol=6)
#model, r2, ci, low, high p-value
#mode1: PRS, PRS+age, PRS+ethnicty
###################################
#######################################
#r2_cofunder
#######################################
myr2_95CI <- function(mydata, mymodel){
  set.seed(1234)
  library(boot)
  r2 <- function(mydata, idx){
    fit<- lm( mymodel, data=mydata[idx, ])
    summary(fit)$r.square
  }
  myfun <- boot(mydata,r2, R=1000) #permutation, 100
  
  fit<- lm( mymodel, data=mydata)
  tmpr2 <- summary(fit)$r.square
 
  if(tmpr2 > 0.0001) { tmpr2_formated<- round(tmpr2,4)}
	else {tmpr2_formated <- scientific(tmpr2, digits=4) } 

  tmppv<- scientific(summary(fit)$coef[2,4], digits=4)
  
  tmpsd <- sd(myfun$t)
  tmplow <- tmpr2 - 1.96 * tmpsd
  tmphigh <- tmpr2 + 1.96 * tmpsd
  if(tmplow > 0.0001) { tmplow <- round(tmplow,4)}
	else {tmplow <- scientific(tmplow, digits=4) } 

  if(tmphigh > 0.0001) { tmphigh <- round(tmphigh,4)}
	else {tmphigh <- scientific(tmphigh, digits=4) } 

  tmpci<-paste0("(", tmplow, ",", tmphigh, ")")
  
  tmpresult <-c(tmpr2_formated, tmpci, tmplow, tmphigh, tmppv)
  return(tmpresult)
}


###############################################
#ObsMetS~otherPRS(include GenMets)
#ObsMetS~GenMetS+age
#ObsMetS~GenMetS+age+ethnicity
#alltraits~GenMetS
#################################################
#all, female, male
pgs20<-fread("PGS20.txt")
colidx <- which(colnames(myraw) %in% pgs20$PGS20)
Nrow <- length(colidx)
cat("we are studying ", Nrow, " PGSs\n")

tmpsize=(3*Nrow*20+10)*11
r2_result <- matrix(rep(NA, tmpsize ), ncol=11)
colnames(r2_result)<- c("target", "name", "N", "sex", "predictor", "R2", "95%CI", "low", "high", "p-value", "cohort")

length(which(!is.na(myraw$MetS_score_category)))
length(which(!is.na(myfemale$MetS_score_category)))
length(which(!is.na(mymale$MetS_score_category)))

#ObsMetS~20traits male/female
K=1
for(i in c(colidx)){
myprs <- colnames(myraw)[i]
mymodel1<- paste0('MetS_score_category',  '~', myprs)
tmp1 <- myr2_95CI(myfemale, mymodel1)
r2_result[K, ]<-c("ObsMetS", "ObsMetS", 246, "female", myprs, tmp1, "ATTRaCT");  K=K+1; rm(tmp1)

tmp2 <- myr2_95CI(mymale, mymodel1)
r2_result[K, ]<-c("ObsMetS", "ObsMetS", 254, "male", myprs, tmp2, "ATTRaCT");  K=K+1; rm(tmp2)
}

#ObsMetS~GenMetS+age, GenMetS+ethnicity
mymodel2<- paste0('MetS_score_category',  '~', 'PRS_GenMetS', '+ AGE')
mymodel3<- paste0('MetS_score_category',  '~', 'PRS_GenMetS', '+ AGE + Ethnicity')

tmp2 <- myr2_95CI(myfemale, mymodel2)
r2_result[K, ]<-c("ObsMetS", "ObsMetS", 246, "female", "GenMetS+age", tmp2, "ATTRaCT");  K=K+1; rm(tmp2)
 
tmp3 <- myr2_95CI(myfemale, mymodel3)
r2_result[K, ]<-c("ObsMetS", "ObsMetS", 246, "female", "GenMetS+age+ethnicity", tmp3, "ATTRaCT");  K=K+1; rm(tmp3)

tmp2 <- myr2_95CI(mymale, mymodel2)
r2_result[K, ]<-c("ObsMetS", "ObsMetS", 254, "male", "GenMetS+age", tmp2, "ATTRaCT");  K=K+1; rm(tmp2)

tmp3 <- myr2_95CI(mymale, mymodel3)
r2_result[K, ]<-c("ObsMetS", "ObsMetS", 254, "male", "GenMetS+age+ethnicity ", tmp3, "ATTRaCT"); rm(tmp3)

##
r2_result[1:3, ]
dim(r2_result)
r2_result <- r2_result %>% data.frame()

idx<-which(!is.na(r2_result$target)) 
if(length(idx)> 0) r2_result <- r2_result[idx, ]

idx<- which(r2_result$predictor == "MetS_PRS")
r2_result$predictor[idx]<- "GenMetS"
write.table(r2_result, "R2_result/result_R2_20PRS_to_ObsMetS_ATTRaCT.txt", row.names=F, sep="\t", quote=F)

########################################################
other_trait_r2_result <- matrix(rep(NA, 2*Nrow*22*11 ), ncol=11)
colnames(other_trait_r2_result)<- c("target", "name", "N", "sex", "predictor", "R2", "95%CI", "low", "high", "p-value", "cohort")

mytraits<-c("BMI", "HEIGHT", "HR", "ECG_HR", "ECHO_HR",
            "SBP", "DBP","V11", "CRET", "HEMO", "WHITE_COUNT","LYMPHOCYTE", 
            "FastingGlu", "CHOLESTROL", "LDL", "HDL", "THDL", "TG",
            "FHL2","HSTNT", "NTPROUCNII")
ntraits <- length(mytraits)

K=1 
for(I in c(1:ntraits) )
{
  mytrait1=mytraits[I]; myname1 =mytraits[I];mystage1=mytraits[I];
  
  mymodel1<- paste0(mytrait1,  '~',  'PRS_GenMetS')

  tmp1 <- myr2_95CI(myfemale, mymodel1)
  idx1<-which(colnames(myfemale)== mytrait1);  
  N= length(which(!is.na(myfemale[, c(idx1)])))
  other_trait_r2_result[K, ]<-c(mytrait1,myname1,N, "female", "PRS_GenMetS", tmp1, "ATTRaCT"); K<- K+1; rm(tmp1)
  
  tmp2 <- myr2_95CI(mymale, mymodel1 )
  idx1<-which(colnames(mymale)== mytrait1);  
  N= length(which(!is.na(mymale[, c(idx1)])))
  other_trait_r2_result[K, ]<-c(mytrait1,myname1, N, "male", "PRS_GenMetS", tmp2, "ATTRaCT"); K<- K+1; rm(tmp2)
  
}
 
other_trait_r2_result <- other_trait_r2_result %>% data.frame()
idx<-which(!is.na(other_trait_r2_result$target))
if(length(idx) > 0) other_trait_r2_result <- other_trait_r2_result[idx, ]
dim(other_trait_r2_result)

idx<-which(other_trait_r2_result$predictor=="MetS_PRS")
if(length(idx) > 0) other_trait_r2_result$predictor[idx]<- "GenMetS"
write.table(other_trait_r2_result, "R2_result/result_R2_20PRS_to_otherTraits_ATTRaCT.txt", row.names=F, sep="\t", quote=F)


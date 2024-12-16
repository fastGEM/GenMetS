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
study_data<- fread("../../0_data/attract_mets_data.txt") %>% data.frame()
PGS20 <- fread("../../0_data/studysample_20PGS.txt")
study_data <- merge(study_data, PGS20, by.x="IID", by.y="sampleID")

study_data_women <- study_data %>% filter(GENDER=="Female"); dim(study_data_women)
study_data_men <- study_data %>% filter(GENDER=="Male"); dim(study_data_men)

N_all  <- dim(study_data)[1]
N_women <- dim(study_data_women)[1]
N_men  <- dim(study_data_men)[1]


###1. R2 for GenMetS to MetS traits
### standard deviation is from 1000 permutations.
library(boot)
##########################################myr2_95CI
myr2_95CI <- function(mydata, mymodel){
  set.seed(1234)
  r2 <- function(mydata, idx){
    fit<- lm( mymodel, data=mydata[idx, ])
    summary(fit)$r.square
  }
  
  myfun <- boot(mydata,r2, R=1000) #permutation, 1000
  #x<- boot.ci(myfun, type="norm")
  #tmplower <- x$normal[2]
  #tmpupper <- x$normal[3]
  
  x<- boot.ci(myfun, type="perc")
  tmplower <- x$percent[4]
  tmpupper <- x$percent[5]
  
  
  fit<- lm( mymodel, data=mydata)
  tmpr2 <- summary(fit)$r.square
  tmppv<- scientific(summary(fit)$coef[2,4], digits=4)
  
  #tmpsd <- sd(myfun$t)
  #tmplow <-  tmpr2 - 1.96 * tmpsd
  #tmphigh <- tmpr2 + 1.96 * tmpsd
  
  #if(tmpr2 > 0.0001) { tmpr2_formated<- round(tmpr2,4)}
  #else {tmpr2_formated <- scientific(tmpr2, digits=4) } 
  
  #if(tmplow > 0.0001) { tmplow <- round(tmplow,4)}
  #else {tmplow <- scientific(tmplow, digits=4) } 
  
  #if(tmphigh > 0.0001) { tmphigh <- round(tmphigh,4)}
  #else {tmphigh <- scientific(tmphigh, digits=4) } 
  
  
  tmpr2_formated <- formatC(tmpr2, format = "g", digits = 3)
  tmplower <- formatC(tmplower, format = "g", digits = 3)
  tmpupper <- formatC(tmpupper, format = "g", digits = 3)
  
  tmpci<-paste0("(", tmplower, ",", tmpupper, ")")
  
  tmpresult <-c(tmpr2_formated, tmpci, tmplower, tmpupper, tmppv)
  return(tmpresult)
}


###############################################
#ObsMetS~otherPRS(include GenMets)
#ObsMetS~GenMetS+age
#ObsMetS~GenMetS+age+ethnicity
#alltraits~GenMetS
#################################################
#all, women, men
pgslist<-fread("../../0_data/PGSlist.txt")
colidx <- which(colnames(study_data) %in% pgslist$PGS)
Nrow <- length(colidx)
cat("we are studying ", Nrow, " PGSs\n")

tmpsize=(3*Nrow*20+10)*11
r2_result <- matrix(rep(NA, tmpsize ), ncol=11)
colnames(r2_result)<- c("target", "name", "N", "sex", "predictor", "R2", "95%CI", "low", "high", "p-value", "cohort")

length(which(!is.na(study_data$MetS_score_category)))
length(which(!is.na(study_data_women$MetS_score_category)))
length(which(!is.na(study_data_men$MetS_score_category)))

#ObsMetS~20traits men/women
K=1
for(i in c(colidx)){
myprs <- colnames(study_data)[i]
mymodel1<- paste0('MetS_score_category',  '~', myprs)
tmp1 <- myr2_95CI(study_data_women, mymodel1)
r2_result[K, ]<-c("ObsMetS", "ObsMetS", N_women, "women", myprs, tmp1, "ATTRaCT");  K=K+1; rm(tmp1)

tmp2 <- myr2_95CI(study_data_men, mymodel1)
r2_result[K, ]<-c("ObsMetS", "ObsMetS", N_men, "men", myprs, tmp2, "ATTRaCT");  K=K+1; rm(tmp2)

tmp3 <- myr2_95CI(study_data, mymodel1)
r2_result[K, ]<-c("ObsMetS", "ObsMetS", N_all, "all", myprs, tmp3, "ATTRaCT");  K=K+1; rm(tmp3)

}

#ObsMetS~GenMetS+age, GenMetS+ethnicity
mymodel2<- paste0('MetS_score_category',  '~', 'PRS_GenMetS', '+ AGE')
mymodel3<- paste0('MetS_score_category',  '~', 'PRS_GenMetS', '+ AGE + Ethnicity')

tmp2 <- myr2_95CI(study_data_women, mymodel2)
r2_result[K, ]<-c("ObsMetS", "ObsMetS", N_women, "women", "GenMetS+age", tmp2, "ATTRaCT");  K=K+1; rm(tmp2)
 
tmp3 <- myr2_95CI(study_data_women, mymodel3)
r2_result[K, ]<-c("ObsMetS", "ObsMetS", N_women, "women", "GenMetS+age+ethnicity", tmp3, "ATTRaCT");  K=K+1; rm(tmp3)

tmp2 <- myr2_95CI(study_data_men, mymodel2)
r2_result[K, ]<-c("ObsMetS", "ObsMetS", N_men, "men", "GenMetS+age", tmp2, "ATTRaCT");  K=K+1; rm(tmp2)

tmp3 <- myr2_95CI(study_data_men, mymodel3)
r2_result[K, ]<-c("ObsMetS", "ObsMetS", N_men, "men", "GenMetS+age+ethnicity ", tmp3, "ATTRaCT"); rm(tmp3)

tmp2 <- myr2_95CI(study_data, mymodel2)
r2_result[K, ]<-c("ObsMetS", "ObsMetS", N_all, "all", "GenMetS+age", tmp2, "ATTRaCT");  K=K+1; rm(tmp2)

tmp3 <- myr2_95CI(study_data, mymodel3)
r2_result[K, ]<-c("ObsMetS", "ObsMetS", N_all, "all", "GenMetS+age+ethnicity ", tmp3, "ATTRaCT"); rm(tmp3)


##
r2_result[1:3, ]
dim(r2_result)
r2_result <- r2_result %>% data.frame()

idx<-which(!is.na(r2_result$target)) 
if(length(idx)> 0) r2_result <- r2_result[idx, ]

write.table(r2_result, "../../0_data/R2_result/result_R2_20PRS_to_ObsMetS_ATTRaCT.txt", row.names=F, sep="\t", quote=F)

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

  tmp1 <- myr2_95CI(study_data_women, mymodel1)
  idx1<-which(colnames(study_data_women)== mytrait1);  
  N= length(which(!is.na(study_data_women[, c(idx1)])))
  other_trait_r2_result[K, ]<-c(mytrait1,myname1,N, "women", "PRS_GenMetS", tmp1, "ATTRaCT"); K<- K+1; rm(tmp1)
  
  tmp2 <- myr2_95CI(study_data_men, mymodel1 )
  idx1<-which(colnames(study_data_men)== mytrait1);  
  N= length(which(!is.na(study_data_men[, c(idx1)])))
  other_trait_r2_result[K, ]<-c(mytrait1,myname1, N, "men", "PRS_GenMetS", tmp2, "ATTRaCT"); K<- K+1; rm(tmp2)
  
  tmp3 <- myr2_95CI(study_data , mymodel1 )
  idx1<-which(colnames(study_data )== mytrait1);  
  N= length(which(!is.na(study_data[, c(idx1)])))
  other_trait_r2_result[K, ]<-c(mytrait1,myname1, N, "all", "PRS_GenMetS", tmp3, "ATTRaCT"); K<- K+1; rm(tmp3)
  
}
 
other_trait_r2_result <- other_trait_r2_result %>% data.frame()
idx<-which(!is.na(other_trait_r2_result$target))
if(length(idx) > 0) other_trait_r2_result <- other_trait_r2_result[idx, ]
dim(other_trait_r2_result)

write.table(other_trait_r2_result, "../../0_data/R2_result/result_R2_20PRS_to_otherTraits_ATTRaCT.txt", row.names=F, sep="\t", quote=F)


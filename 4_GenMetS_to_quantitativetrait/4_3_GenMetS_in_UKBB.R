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
#data
#####################
study_data <- fread("../../0_data/UKB_s365495_mets_diseases_others.txt")

study_data$ID <-as.character(study_data$ID)
PGS20 <- fread("../../0_data/studysample_20PGS.txt")
study_data<- merge(study_data, PGS20, by.x="ID", by.y="sampleID")
colnames(study_data)[which(colnames(study_data)=="Age")] <- "age"

study_data_ASN <- study_data %>% dplyr::filter(ethnicity %in% c("Chinese", "Indian", "OtherAsian"));
dim(study_data_ASN) #8792
num_study_data_ASN <- dim(study_data_ASN)[1]

study_data_ASN_women <- study_data %>% dplyr::filter(sex=="Female" & ethnicity %in% c("Chinese", "Indian", "OtherAsian")); 
dim(study_data_ASN_women) #4189
num_study_data_ASN_women <- dim(study_data_ASN_women)[1]

study_data_ASN_men <- study_data %>% dplyr::filter(sex=="Male" & ethnicity %in% c("Chinese", "Indian", "OtherAsian")); 
dim(study_data_ASN_men)   #4603
num_study_data_ASN_men <- dim(study_data_ASN_men)[1]

study_data_EUR <- study_data %>% dplyr::filter(ethnicity == "European"); 
dim(study_data_EUR) #356703
num_study_data_EUR <- dim(study_data_EUR)[1]

study_data_EUR_women <- study_data %>% dplyr::filter(sex=="Female"& ethnicity == "European"); 
dim(study_data_EUR_women) #191,476
num_study_data_EUR_women <- dim(study_data_EUR_women)[1]

study_data_EUR_men <- study_data %>% dplyr::filter(sex=="Male"& ethnicity == "European"); 
dim(study_data_EUR_men) #165,227
num_study_data_EUR_men <- dim(study_data_EUR_men)[1]

library(boot)
##########################################myr2_95CI
myr2_95CI <- function(mydata, mymodel){
  set.seed(1234)
  r2 <- function(mydata, idx){
    fit<- lm( mymodel, data=mydata[idx, ])
    summary(fit)$r.square
  }
  
  myfun <- boot(mydata,r2, R=1000 ) #permutation, 1000
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


#UKBAsian, female, male
pgslist<-fread("../../0_data/PGSlist.txt")
colidx <- which(colnames(study_data) %in% pgslist$PGS)
Nrow <- length(colidx)
cat("we are studying ", Nrow, " PGSs\n")
tmpsize=(3*Nrow*20 +10)*11
r2_result <- matrix(rep(NA, tmpsize ), ncol=11)
colnames(r2_result)<- c("target", "name", "N", "sex", "predictor", "R2", "95%CI", "low", "high", "p-value", "cohort")

###############################################
#ObsMetS~otherPRS(include GenMets)
#ObsMetS~GenMetS+age
#ObsMetS~GenMetS+age+ethnicity
#alltraits~GenMetS
#################################################
#ObsMetS~otherPRS(include GenMets)
K=1
for(i in c(colidx))
{
  myprs <- colnames(study_data)[i]
  mymodel1<- paste0('ObsMetS',  '~', myprs)
  
  #ASN
  tmp1 <- myr2_95CI(study_data_ASN_women, mymodel1)
  r2_result[K, ]<-c("ObsMetS", "ObsMetS", num_study_data_ASN_women, "women", myprs, tmp1, "UKBB Asian"); K=K+1; rm(tmp1)

  tmp2 <- myr2_95CI(study_data_ASN_men, mymodel1)
  r2_result[K, ]<-c("ObsMetS", "ObsMetS", num_study_data_ASN_men, "men", myprs, tmp2, "UKBB Asian"); K=K+1; rm(tmp2)
  
  tmp3 <- myr2_95CI(study_data_ASN, mymodel1)
  r2_result[K, ]<-c("ObsMetS", "ObsMetS", num_study_data_ASN, "all", myprs, tmp3, "UKBB Asian"); K=K+1; rm(tmp3)

  #EUR 
  tmp1 <- myr2_95CI(study_data_EUR_women, mymodel1)
  r2_result[K, ]<-c("ObsMetS", "ObsMetS", num_study_data_EUR_women, "women", myprs, tmp1, "UKBB European"); K=K+1; rm(tmp1)
  
  tmp2 <- myr2_95CI(study_data_EUR_men, mymodel1)
  r2_result[K, ]<-c("ObsMetS", "ObsMetS", num_study_data_EUR_men, "men", myprs, tmp2, "UKBB European"); K=K+1; rm(tmp2)
  
  tmp3 <- myr2_95CI(study_data_EUR , mymodel1)
  r2_result[K, ]<-c("ObsMetS", "ObsMetS", num_study_data_EUR , "all", myprs, tmp3, "UKBB European"); K=K+1; rm(tmp3)
  
}


#ObsMetS~GenMetS+age
#ObsMetS~GenMetS+age+ethnicity
mymodel2<- paste0('ObsMetS',  '~', 'PRS_GenMetS + age')
mymodel3<- paste0('ObsMetS',  '~', 'PRS_GenMetS + age + ethnicity')

  tmp1 <- myr2_95CI(study_data_ASN_women, mymodel2)
  r2_result[K, ]<-c("ObsMetS", "ObsMetS", num_study_data_ASN_women, "women", "GenMetS+age", tmp1, "UKBB Asian")
  print(r2_result[K, ]); K=K+1; rm (tmp1) 

  tmp2 <- myr2_95CI(study_data_ASN_men, mymodel2)
  r2_result[K, ]<-c("ObsMetS", "ObsMetS", num_study_data_ASN_men, "men", "GenMetS+age", tmp2, "UKBB Asian"); 
  print(r2_result[K, ]);  K=K+1; rm(tmp2)

  tmp3 <- myr2_95CI(study_data_ASN, mymodel2)
  r2_result[K, ]<-c("ObsMetS", "ObsMetS", num_study_data_ASN, "all", "GenMetS+age", tmp3, "UKBB Asian"); 
  print(r2_result[K, ]);  K=K+1; rm(tmp3)
  
  tmp1 <- myr2_95CI(study_data_ASN_women, mymodel3)
  r2_result[K, ]<-c("ObsMetS", "ObsMetS", num_study_data_ASN_women, "women", "GenMetS+age+ethnicity", tmp1, "UKBB Asian"); 
  print(r2_result[K, ]);  K=K+1; rm(tmp1)

  tmp2 <- myr2_95CI(study_data_ASN_men, mymodel3)
  r2_result[K, ]<-c("ObsMetS", "ObsMetS", num_study_data_ASN_men, "men", "GenMetS+age+ethnicity", tmp2, "UKBB Asian"); 
  print(r2_result[K, ]);  K=K+1; rm(tmp2)

  tmp3 <- myr2_95CI(study_data_ASN, mymodel3)
  r2_result[K, ]<-c("ObsMetS", "ObsMetS", num_study_data_ASN, "all", "GenMetS+age+ethnicity", tmp3, "UKBB Asian"); 
  print(r2_result[K, ]);  K=K+1; rm(tmp3)
  

r2_result <- r2_result %>% data.frame()
r2_result[1:3, ]
dim(r2_result)

idx<-which(!is.na(r2_result$target) )
if(length(idx)>0) r2_result <- r2_result[idx, ]
write.table(r2_result, "../../0_data/R2_result/result_R2_20PRS_to_ObsMetS_UKBB.txt", row.names=F, sep="\t", quote=F)


############################### 
mytraits<-c( "WC", "SBP", "DBP", "TG", "HDL","HbA1c", "BMI")

nTrait<- length(mytraits)
other_trait_r2_result <- matrix(rep(NA, 4*Nrow*nTrait*11 ), ncol=11)
colnames(other_trait_r2_result)<- c("target", "name", "N", "sex", "predictor", "R2", "95%CI", "low", "high", "p-value", "cohort")

K=1 
for(I in c(1:nTrait) )
{
  mytrait1=mytraits[I]; myname1 =mytraits[I];mystage1=mytraits[I];

  colidx <- which(colnames(study_data)=="PRS_GenMetS")
  myprs <- colnames(study_data)[colidx]
  mymodel1<- paste0(mytrait1,  '~', myprs)
  
  mydata <- study_data_ASN_women %>% data.frame()
  tmp1 <- myr2_95CI(mydata, mymodel1)
  idx1<-which(colnames(mydata)== mytrait1);  
  N= length(which(!is.na(mydata[, c(idx1)])))
  other_trait_r2_result[K, ]<-c(mytrait1,myname1, N, "women", myprs, tmp1, "UKBB Asian")
  cat("K=", K, "\t", other_trait_r2_result[K, ], "\n"); K<- K+1; rm(tmp1)
  
  mydata <- study_data_ASN_men %>% data.frame()
  tmp1 <- myr2_95CI(mydata, mymodel1)
  idx1<-which(colnames(mydata)== mytrait1);  
  N= length(which(!is.na(mydata[, c(idx1)])))
  other_trait_r2_result[K, ]<-c(mytrait1,myname1, N, "men", myprs, tmp1, "UKBB Asian")
  cat("K=", K, "\t", other_trait_r2_result[K, ], "\n"); K<- K+1; rm(tmp1)
  
  mydata <- study_data_ASN %>% data.frame()
  tmp1 <- myr2_95CI(mydata, mymodel1)
  idx1<-which(colnames(mydata)== mytrait1);  
  N= length(which(!is.na(mydata[, c(idx1)])))
  other_trait_r2_result[K, ]<-c(mytrait1,myname1, N, "all", myprs, tmp1, "UKBB Asian")
  cat("K=", K, "\t", other_trait_r2_result[K, ], "\n"); K<- K+1; rm(tmp1)
  
  mydata <- study_data_EUR_women %>% data.frame()
  tmp1 <- myr2_95CI(mydata, mymodel1)
  idx1<-which(colnames(mydata)== mytrait1);  
  N= length(which(!is.na(mydata[, c(idx1)])))
  other_trait_r2_result[K, ]<-c(mytrait1,myname1, N, "women", myprs, tmp1, "UKBB European")
  cat("K=", K, "\t", other_trait_r2_result[K, ], "\n"); K<- K+1; rm(tmp1)
  
  mydata <- study_data_EUR_men %>% data.frame()
  tmp1 <- myr2_95CI(mydata, mymodel1)
  idx1<-which(colnames(mydata)== mytrait1);  
  N= length(which(!is.na(mydata[, c(idx1)])))
  other_trait_r2_result[K, ]<-c(mytrait1,myname1, N, "men", myprs, tmp1, "UKBB European")
  cat("K=", K, "\t", other_trait_r2_result[K, ], "\n"); K<- K+1; rm(tmp1)
  
  mydata <- study_data_EUR %>% data.frame()
  tmp1 <- myr2_95CI(mydata, mymodel1)
  idx1<-which(colnames(mydata)== mytrait1);  
  N= length(which(!is.na(mydata[, c(idx1)])))
  other_trait_r2_result[K, ]<-c(mytrait1,myname1, N, "all", myprs, tmp1, "UKBB European")
  cat("K=", K, "\t", other_trait_r2_result[K, ], "\n"); K<- K+1; rm(tmp1)
  
  }

other_trait_r2_result <- data.frame(other_trait_r2_result)
idx<-which(!is.na(other_trait_r2_result$target) )
if(length(idx)>0) other_trait_r2_result <- other_trait_r2_result[idx, ]
dim(other_trait_r2_result)
write.table(other_trait_r2_result, "../../0_data/R2_result/result_R2_20PRS_to_otherTraits_UKBB.txt", row.names=F, sep="\t", quote=F)



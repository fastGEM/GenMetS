######################
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
myraw <- fread("UKB_s365495_mets_diseases_others.txt")
table(myraw$ethnicity)

myraw$ID <-as.character(myraw$ID)
myPRS19 <- fread("../1_make_20PGS/studysample_20PRS.txt")
mydata2 <- merge(myraw, myPRS19, by.x="ID", by.y="sampleID")
myraw <- mydata2; rm(mydata2)
colnames(myraw)[which(colnames(myraw)=="Age")] <- "age"

mydata <- myraw
myAsian <- mydata %>% dplyr::filter(ethnicity %in% c("Chinese", "Indian", "OtherAsian"));
dim(myAsian) #8792
num_myAsian <- dim(myAsian)[1]

myAsianfemale <- mydata %>% dplyr::filter(sex=="Female" & ethnicity %in% c("Chinese", "Indian", "OtherAsian")); 
dim(myAsianfemale) #4189
num_myAsianfemale <- dim(myAsianfemale)[1]

myAsianmale <- mydata %>% dplyr::filter(sex=="Male" & ethnicity %in% c("Chinese", "Indian", "OtherAsian")); 
dim(myAsianmale)   #4603
num_myAsianmale <- dim(myAsianmale)[1]

myEuropean <- mydata %>% dplyr::filter(ethnicity == "European"); 
dim(myEuropean) #356703
num_myEuropean <- dim(myEuropean)[1]

myEuropeanfemale <- mydata %>% dplyr::filter(sex=="Female"& ethnicity == "European"); 
dim(myEuropeanfemale) #191,476
num_myEuropeanfemale <- dim(myEuropeanfemale)[1]

myEuropeanmale <- mydata %>% dplyr::filter(sex=="Male"& ethnicity == "European"); 
dim(myEuropeanmale) #165,227
num_myEuropeanmale <- dim(myEuropeanmale)[1]

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

#UKBAsian, female, male
length(which(!is.na(myAsian$GenMetS)))
length(which(!is.na(myAsianfemale$GenMetS)))
length(which(!is.na(myAsianmale$GenMetS)))


pgs20<-fread("PGS20.txt")
colidx <- which(colnames(myraw) %in% pgs20$PGS20)
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
  myprs <- colnames(myraw)[i]
  mymodel1<- paste0('ObsMetS',  '~', myprs)
  
  #ASN
  tmp1 <- myr2_95CI(myAsianfemale, mymodel1)
  r2_result[K, ]<-c("ObsMetS", "ObsMetS", num_myAsianfemale, "female", myprs, tmp1, "UKBB Asian"); K=K+1; rm(tmp1)

  tmp2 <- myr2_95CI(myAsianmale, mymodel1)
  r2_result[K, ]<-c("ObsMetS", "ObsMetS", num_myAsianmale, "male", myprs, tmp2, "UKBB Asian"); K=K+1; rm(tmp2)

  #EUR 
  tmp1 <- myr2_95CI(myEuropeanfemale, mymodel1)
  r2_result[K, ]<-c("ObsMetS", "ObsMetS", num_myEuropeanfemale, "female", myprs, tmp1, "UKBB European"); K=K+1
  
  tmp2 <- myr2_95CI(myEuropeanmale, mymodel1)
  r2_result[K, ]<-c("ObsMetS", "ObsMetS", num_myEuropeanmale, "male", myprs, tmp2, "UKBB European"); K=K+1
  
}


#ObsMetS~GenMetS+age
#ObsMetS~GenMetS+age+ethnicity
mymodel2<- paste0('ObsMetS',  '~', 'PRS_GenMetS + age')
mymodel3<- paste0('ObsMetS',  '~', 'PRS_GenMetS+ age + ethnicity')

  tmp1 <- myr2_95CI(myAsianfemale, mymodel2)
  r2_result[K, ]<-c("ObsMetS", "ObsMetS", num_myAsianfemale, "female", "GenMetS+age", tmp1, "UKBB Asian")
  print(r2_result[K, ]); K=K+1; rm (tmp1) 

  tmp2 <- myr2_95CI(myAsianmale, mymodel2)
  r2_result[K, ]<-c("ObsMetS", "ObsMetS", num_myAsianmale, "male", "GenMetS+age", tmp2, "UKBB Asian"); 
  print(r2_result[K, ]);  K=K+1; rm(tmp2)

  tmp1 <- myr2_95CI(myEuropeanfemale, mymodel2)
  r2_result[K, ]<-c("ObsMetS", "ObsMetS", num_myEuropeanfemale, "female", "GenMetS+age", tmp1, "UKBB European"); 
  print(r2_result[K, ]);  K=K+1; rm(tmp1)
  
  tmp2 <- myr2_95CI(myEuropeanmale, mymodel2)
  r2_result[K, ]<-c("ObsMetS", "ObsMetS", num_myEuropeanmale, "male", "GenMetS+age", tmp2, "UKBB European"); 
  print(r2_result[K, ]);  K=K+1; rm(tmp2)


  tmp1 <- myr2_95CI(myAsianfemale, mymodel3)
  r2_result[K, ]<-c("ObsMetS", "ObsMetS", num_myAsianfemale, "female", "GenMetS+age+ethnicity", tmp1, "UKBB Asian"); 
  print(r2_result[K, ]);  K=K+1; rm(tmp1)

  tmp2 <- myr2_95CI(myAsianmale, mymodel3)
  r2_result[K, ]<-c("ObsMetS", "ObsMetS", num_myAsianmale, "male", "GenMetS+age+ethnicity", tmp2, "UKBB Asian"); 
  print(r2_result[K, ]);  K=K+1; rm(tmp2)

  #UKB EUR no ethnicity

r2_result <- r2_result %>% data.frame()
r2_result[1:3, ]
idx<- which(r2_result$predictor == "PRS_GenMetS")
if(length(idx)>0) r2_result$predictor[idx]<- "GenMetS"
dim(r2_result)

idx<-which(!is.na(r2_result$target) )
if(length(idx)>0) r2_result <- r2_result[idx, ]
write.table(r2_result, "R2_result/result_R2_20PRS_to_ObsMetS_UKBB.txt", row.names=F, sep="\t", quote=F)


############################### 
mytraits<-c( "WC", "SBP", "DBP", "TG", "HDL","HbA1c", "BMI")

nTrait<- length(mytraits)
other_trait_r2_result <- matrix(rep(NA, 4*Nrow*nTrait*11 ), ncol=11)
colnames(other_trait_r2_result)<- c("target", "name", "N", "sex", "predictor", "R2", "95%CI", "low", "high", "p-value", "cohort")

K=1 
for(I in c(1:nTrait) )
{
  mytrait1=mytraits[I]; myname1 =mytraits[I];mystage1=mytraits[I];

  colidx <- which(colnames(myraw)=="PRS_GenMetS")
  myprs <- colnames(myraw)[colidx]
  mymodel1<- paste0(mytrait1,  '~', myprs)
  
  mydata <- myAsianfemale %>% data.frame()
  tmp1 <- myr2_95CI(mydata, mymodel1)
  idx1<-which(colnames(mydata)== mytrait1);  
  N= length(which(!is.na(mydata[, c(idx1)])))
  other_trait_r2_result[K, ]<-c(mytrait1,myname1, N, "female", myprs, tmp1, "UKBB Asian")
  cat("K=", K, "\t", other_trait_r2_result[K, ], "\n"); K<- K+1; rm(tmp1)
  
  mydata <- myAsianmale %>% data.frame()
  tmp2 <- myr2_95CI(mydata, mymodel1)
  idx1<-which(colnames(mydata)== mytrait1);  
  N= length(which(!is.na(mydata[, c(idx1)])))
  other_trait_r2_result[K, ]<-c(mytrait1,myname1, N, "male", myprs, tmp2, "UKBB Asian")
  cat("K=", K, "\t", other_trait_r2_result[K, ], "\n"); K<- K+1; rm(tmp2)

  mydata <- myEuropeanfemale %>% data.frame()
  tmp3<- myr2_95CI(mydata, mymodel1)
  idx1<-which(colnames(mydata)== mytrait1);  
  N= length(which(!is.na(mydata[, c(idx1)])))
  other_trait_r2_result[K, ]<-c(mytrait1,myname1,  N, "female", myprs, tmp3, "UKBB European")
  cat("K=", K, "\t", other_trait_r2_result[K, ], "\n"); K<- K+1; rm(tmp3)
  
  mydata <- myEuropeanmale %>% data.frame()
  tmp4 <- myr2_95CI(mydata, mymodel1)
  idx1<-which(colnames(mydata)== mytrait1);  
  N= length(which(!is.na(mydata[, c(idx1)])))
  other_trait_r2_result[K, ]<-c(mytrait1,myname1, N, "male", myprs, tmp4, "UKBB European")
  cat("K=", K, "\t", other_trait_r2_result[K, ], "\n"); K<- K+1; rm(tmp4)
}

other_trait_r2_result <- data.frame(other_trait_r2_result)
idx<-which(!is.na(other_trait_r2_result$target) )
if(length(idx)>0) other_trait_r2_result <- other_trait_r2_result[idx, ]
dim(other_trait_r2_result)
write.table(other_trait_r2_result, "R2_result/result_R2_20PRS_to_otherTraits_UKBB.txt", row.names=F, sep="\t", quote=F)



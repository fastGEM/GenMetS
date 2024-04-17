############################################
############################################

rm(list=ls())
library(dplyr)
library(data.table)
library(scales) 

mytmp <- fread("GUSTOchildren1073_mets_pc.txt")
colnames(mytmp)[which(colnames(mytmp)=="Age_yr6")]<- "age"
colnames(mytmp)[which(colnames(mytmp)=="GenMetS")]<- "predicted_mets"

myPRS20 <- fread("../1_make_20PGS/studysample_20PRS.txt")
myPRS20$sampleID<- gsub("B", "", myPRS20$sampleID)
mydata2 <- merge(mytmp, myPRS20, by.x="ID", by.y="sampleID")
mytmp <- mydata2; rm(mydata2)

mytmp_girl <- mytmp %>% dplyr::filter(sex=="Female") %>% data.frame()
mytmp_boy <- mytmp %>% dplyr::filter(sex=="Male") %>% data.frame()

tmp1<-mytmp_girl %>% dplyr::filter(!is.na(MetS_child_yr6))
summary(tmp1$age)
table(tmp1$ethnicity)

tmp2<-mytmp_boy %>% dplyr::filter(!is.na(MetS_child_yr6))
summary(tmp2$age)
table(tmp2$ethnicity)


###1. R2 for PredMetS to all related traits.
##########################################myr2_95CI
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
#mytmp, MetS_child_yr6, N=321, both boy and girl
pgs20<-fread("PGS20.txt")
colidx <- which(colnames(mytmp) %in% pgs20$PGS20)
Nrow <- length(colidx)
cat("we are studying ", Nrow, " PGSs\n")

tmpsize=(3*Nrow*20+10) *11
r2_result <- matrix(rep(NA, tmpsize), ncol=11)
colnames(r2_result)<- c("target", "name", "N", "sex", "predictor", "R2", "95%CI", "low", "high", "p-value", "cohort")


###############################################
#ObsMetS~otherPRS(include GenMets)
#ObsMetS~GenMetS+age
#ObsMetS~GenMetS+age+ethnicity
#alltraits~GenMetS
#################################################
#ObsMetS~20PRS
K=1
for(i in colidx)
{
myprs <- colnames(mytmp)[i]
mymodel1<- paste0('MetS_child_yr6',  '~', myprs)

tmp1 <- myr2_95CI(mytmp_girl, mymodel1)
r2_result[K, ]<-c("ObsMetS", "ObsMetS", "163", "girl", myprs, tmp1, "GUSTOchildren6yr"); K=K+1; rm(tmp1)

tmp2 <- myr2_95CI(mytmp_boy, mymodel1)
r2_result[K, ]<-c("ObsMetS", "ObsMetS", "158", "boy", myprs, tmp2, "GUSTOchildren6yr"); K=K+1; rm(tmp2)
}

#ObsMetS~GenMetS+age
#ObsMetS~GenMetS+age+ethnicity
mymodel2<- paste0('MetS_child_yr6',  '~', myprs, '+ age')
mymodel3<- paste0('MetS_child_yr6',  '~', myprs, '+ age + ethnicity')

tmp2 <- myr2_95CI(mytmp_girl, mymodel2)
r2_result[K, ]<-c("ObsMetS", "ObsMetS", "163", "girl", "GenMetS+age", tmp2, "GUSTOchildren6yr"); K=K+1; rm(tmp2)

tmp3 <- myr2_95CI(mytmp_girl, mymodel3)
r2_result[K, ]<-c("ObsMetS", "ObsMetS", "163", "girl", "GenMetS+age+ethnicity", tmp3, "GUSTOchildren6yr"); K=K+1; rm(tmp3)

tmp2 <- myr2_95CI(mytmp_boy, mymodel2)
r2_result[K, ]<-c("ObsMetS", "ObsMetS", "158", "boy", "GenMetS+age", tmp2, "GUSTOchildren6yr"); K=K+1; rm(tmp2)

tmp3 <- myr2_95CI(mytmp_boy, mymodel3)
r2_result[K, ]<-c("ObsMetS", "ObsMetS", "158", "boy", "GenMetS+age+ethnicity ", tmp3, "GUSTOchildren6yr"); rm(tmp3)

r2_result <- r2_result %>% data.frame()
r2_result[1:3, ]
idx<- which(r2_result$predictor == "predicted_mets"); r2_result$predictor[idx]<- "GenMetS"

idx<-which(!is.na(r2_result$target)) 
if(length(idx) > 0) r2_result<- r2_result[idx, ]
dim(r2_result)
write.table(r2_result, "R2_result/result_R2_20PRS_to_ObsMetS_GUSTOchildren.txt", row.names=F, sep="\t", quote=F)

###############################################################################
mytraits<-c("c_abdominal_yr6", "FLI_child_yr6", "TG_mmol_L_child_yr6", 
            "HDL_mmol_L_child_yr6","LDL_mmol_L_child_yr6", "INS_mU_L_child_yr6", 
            "HOMA_IR_child_yr6", "c_weight_birth","neonate_SAT_18Sept", 
            "neonate_dSAT_18Sept", "neonate_dSAT_18Sept",
            "Cord_Insulin_pgml", "Cord_Cpeptide_pgml", 
            "BMI_child_yr6",  "logHSCRP", "SAT_yr6", "DSAT_yr6", "SSAT_yr6", "VAT_yr6", 
	    "Liver_fat_yr6", "IMCL_yr6", "MetS_child_yr6")

mytrait_names<-c("Abdominal Circumference", "Fatty Liver Index", "TG(mmol/L)", "HDL(mmol/L)",
            "LDL(mmol/L)", "Insulin(mU/L)", "HOMA_IR", "BirthWeight", "neonate_SAT", "neonate_dSAT", "neonate_IAT",
            "Cord_Insulin(pgml)", "Cord_Cpeptide(pgml)", "BMI",  "logHSCRP(mg/L)", 
            "SAT_yr6", "DSAT_yr6", "SSAT_yr6", "VAT_yr6", "Liver_fat_yr6", "IMCL_yr6", "ObSMetS")


mytrait_stage<-c( "yr6", "yr6", "yr6", "yr6","yr6", "yr6", "yr6", "birth",
             "birth", "birth", "birth","birth", "birth", "yr6", "yr6", "yr6", "yr6", "yr6", "yr6", "yr6", "yr6")

Ntrait=length(mytraits)
other_trait_r2_result <- matrix(rep(NA, 2*Nrow*Ntrait*12 ), ncol=12)
colnames(other_trait_r2_result)<- c("target", "name", "stage", "N", "sex", "predictor", "R2", "95%CI", "low", "high", "p-value", "cohort")


K=1 
for(I in c(1:Ntrait) )
{
  mytrait1=mytraits[I]; myname1 =mytrait_names[I];mystage1=mytrait_stage[I];
  
  colidx <- which(colnames(mytmp) == "PRS_GenMetS")  

    myprs <- colnames(mytmp)[colidx]
    mymodel1<- paste0(mytrait1,  '~', myprs)
  
    tmp1 <- myr2_95CI(mytmp_girl, mymodel1)
    idx1<-which(colnames(mytmp_girl)== mytrait1);  
    N= length(which(!is.na(mytmp_girl[, c(idx1)])))
    other_trait_r2_result[K, ]<-c(mytrait1,myname1, mystage1, N, "girl", myprs, tmp1, "GUSTOchildren"); K<- K+1; rm(tmp1)

    tmp2 <- myr2_95CI(mytmp_boy, mymodel1)
    idx1<-which(colnames(mytmp_boy)== mytrait1);  
    N= length(which(!is.na(mytmp_boy[, c(idx1)])))
    other_trait_r2_result[K, ]<-c(mytrait1,myname1, mystage1, N, "boy", myprs, tmp2, "GUSTOchildren"); K<- K+1; rm(tmp2)
    K<- K+1
}

other_trait_r2_result<- other_trait_r2_result %>% data.frame()
idx<-which(!is.na(other_trait_r2_result$target)) 
if(length(idx) > 0) other_trait_r2_result<- other_trait_r2_result[idx, ]
dim(other_trait_r2_result)

idx<-which(other_trait_r2_result$predictor == "predicted_mets")
other_trait_r2_result$predictor[idx]<- "GenMetS"
write.table(other_trait_r2_result, "R2_result/result_R2_20PRS_to_otherTraits_GUSTOchildren.txt",row.names=F, sep="\t", quote=F)




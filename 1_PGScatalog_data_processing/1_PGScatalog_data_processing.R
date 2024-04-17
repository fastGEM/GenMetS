#######################
#use system() to run
#bash command
######################

library(data.table)
library(dplyr)

myinput="PGSlist_for_SNP_aggregation.txt"
myoutputdir="1_processed_data/"


mydata<- fread(myinput) %>% data.frame()
myPGSid <- mydata$PRSid
num_PGSid <- length(myPGSid)
cat("there are ", num_PGSid, "PGSids\n")

##################################################################download PGScatalog data
for (i in c(1:num_PGSid))
{
  file=paste0(myPGSid[i], ".txt")
  tmp=paste0(myoutputdir, file)

  if(! file.exists( tmp ) )
{
  cat ("try ", file, "\n")

  myurl= paste0("https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/", myPGSid[i], "/ScoringFiles/", file, ".gz")
  cmd1 = paste0("wget ", myurl)
  cmd2 = paste0("gunzip ", file, ".gz")
  cmd3 = paste0("mv ", file, " ", myoutputdir)
  system(cmd1)
  system(cmd2)
  system(cmd3)
}
}


#######################################################################PGS data QC
myinput="PGSlist_for_SNP_aggregation.txt"
myinputdir="1_processed_data/"

mydata<- read.csv(myinput, header=T, sep="\t")
myPGSid <- mydata$PRSid
num_PGSid <- length(myPGSid)
cat("there are ", num_PGSid, "PGSids\n")

for (i in c(1:num_PGSid))
{
myinputfile =paste0(myinputdir, myPGSid[i], ".txt")
myoutputfile=gsub(".txt", ".QC.txt", myinputfile);

mydata <- fread(myinputfile)
mydata[1:3, ]

flag=0						#at least three basic columns. 
if(c("rsID")          %in% colnames(mydata) && 
   c("effect_allele") %in% colnames(mydata) &&
   c("effect_weight") %in% colnames(mydata) ) {flag=1; }
if(flag==0){ cat("failed in 3 column check = ", myinputfile, "\n")
next;}

#For redundant SNPs, select the biggest effect
x<- mydata %>% group_by(rsID) %>% filter(effect_weight == max(effect_weight)) 
dim(mydata); dim(x)
mydata<-x

colnames(mydata)[which(colnames(mydata)=="rsID")]<- "SNP"
colnames(mydata)[which(colnames(mydata)=="effect_allele")]<- "A1"
colnames(mydata)[which(colnames(mydata)=="effect_weight")]<- "BETA"
mydata$A1 <- toupper(mydata$A1)

if(c("reference_allele") %in% colnames(mydata)) {
colnames(mydata)[which(colnames(mydata)=="reference_allele")]<- "A2"
mydata$A2 <- toupper(mydata$A2)
myout<- mydata[, c("SNP", "A1", "A2",  "BETA")] 
} else{ 
myout<- mydata[, c("SNP", "A1",  "BETA")] 
}

#2. QC
#remove SNP="."
mydata1 <- myout %>% filter(!(SNP %in% c(".", "") ) )

#remove ambigious if A2
if(c("A2") %in% colnames(myout)) {
mydata2 <- mydata1 %>% filter(
		(A1=="A" & A2=="T") |
		(A1=="T" & A2=="A") |
		(A1=="G" & A2=="C") |
		(A1=="C" & A2=="G") )

mydata3<- mydata1 %>% filter(!(SNP %in% mydata2$SNP))
fwrite(mydata3, myoutputfile, sep="\t")
} else{
fwrite(mydata1, myoutputfile, sep="\t")
}

}

#######################################################PGS data aggregation
inputfile="PGSlist_for_SNP_aggregation.txt"
PRSQCdir="1_processed_data/"
outputfile="2_result/SNPlist.txt"	 

myPGS<-read.csv(inputfile, header=T, sep="\t", comment.char="#")
myPGS[1:3, ]

myfile1=paste0(PRSQCdir,myPGS$PRSid[1], ".QC.txt")
cat("myfile = ", myfile1, "\n")

mydata <-data.frame(fread(myfile1,header=T,sep="\t"))
mydata<- mydata[, c("SNP", "A1", "BETA")]
mydata$SNP_ALT<- paste0(mydata$SNP,"_", mydata$A1)
mydata$tag <- myPGS$Phenotype1[1]
mydata$BETA<- scale(mydata$BETA)
mydata$PGSID <- myPGS$PRSid[1]

for(I in c(2:dim(myPGS)[1]))
{
  myfile1=paste0(PRSQCdir, myPGS$PRSid[I], ".QC.txt")
  cat("myfile = ", myfile1, "\n")
  
  tmp <-data.frame(fread(myfile1,header=T,sep="\t"))
  dim(tmp) 

  tmp<- tmp[, c("SNP", "A1", "BETA")]
  tmp$SNP_ALT<- paste0(tmp$SNP,"_", tmp$A1)
  tmp$tag <- myPGS$Phenotype1[I]
  tmp$BETA<- scale(tmp$BETA)
  tmp$PGSID <- myPGS$PRSid[I]
  
  mydata <- rbind(mydata,tmp); rm(tmp)
}

dim(mydata)
mydata[1:3, ]

table(mydata$tag)
y<- x$SNP_ALT
y<- gsub("[_AGCT]", "", y)
write.table(y, outputfile, row.names=F, col.names=F, sep="\t", quote=F)
cat("output = ", outputfile, "\n")


#######################################################PGS meta information
rm(list=ls())
library(data.table)
library(dplyr)
library(readxl)
library(tidyverse)
library(stringr)


metadatafile= "pgs_all_metadata_6Jan2022.xlsx"
myPGSfile="PGSlist_for_SNP_aggregation.txt"

mydata1<- read_excel(metadatafile, sheet=1)
mydata2<- read_excel(metadatafile, sheet=2)
mydata3<- read_excel(metadatafile, sheet=3)
mydata4<- read_excel(metadatafile, sheet=4)
mydata5<- read_excel(metadatafile, sheet=5)
mydata6<- read_excel(metadatafile, sheet=6)
mydata7<- read_excel(metadatafile, sheet=7)
mydata8<- read_excel(metadatafile, sheet=8)

colnames(mydata2)

mydata2<- mydata2[, c(1,4,5,8)]
colnames(mydata2)<- c("PGPID", "Journal", "Date", "PMID") #PGP

colnames(mydata4)
mydata4<- mydata4[, c(1,3,9, 11,15,17)]
colnames(mydata4)<- c("PGSID", "Reported Trait", "Number of Variants", "PGPID", 
                      "Ancestry_GWAS", "Ancestry_Evaluation")

colnames(mydata5)
mydata5<- mydata5[, c(1,3)]
colnames(mydata5)<- c("PGSID", "GWAS_Samples")

colnames(mydata7)
mydata7<- mydata7[, c(2,3)]
colnames(mydata7)<- c("PGSID", "Evaluation_Samples")

tmp1<- merge(mydata5, mydata7, by.x="PGSID", by.y="PGSID", all.x=T)
tmp2<- merge(mydata4, tmp1, by.x="PGSID", by.y="PGSID")
tmp3<- merge(tmp2, mydata2, by.x="PGPID", by.y="PGPID", all.x=T)
mydata<- tmp3 
rm(tmp1, tmp2, tmp3)
dim(mydata)
colnames(mydata)

myPGS<- read.csv(myPGSfile, header=T, sep="\t")
dim(myPGS)

idx<- match(myPGS$PRSid, mydata$PGSID)
myout <- mydata[idx, ]
dim(myout)
myout <- myout %>% data.frame() %>% dplyr::select(-PGPID)

myout[1:3, ]
write.table(myout, "2_result/Mets_PRS_in_metadata.txt", row.names=F, sep="\t", quote=F)


##############
#Pan Hong
#Feb16,2024
#last visit Mar15,2024
#upon BJJ returned 
#the result.
##############
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

mydisease<-c("T2D", "CAD", "HTN","HF", "stroke", "MI")

attract <- fread("sup_table_ATTRaCT_GenMetS.txt")
ukb_asn<- fread("sup_table_UKB_ASN_GenMetS.txt")


#################BBJ

bbj0<-fread("../../0_data/BBJ/BBJ.binary_trait.demographic.txt") 
###female
bbj1<- fread("../../0_data/BBJ/GenMetS4627.binary_trait.female.txt")
bbj1$Cohort <- "BBJ"
bbj1$Sex <- "Women"
bbj0female <- bbj0[, c(1,2)]; colnames(bbj0female)<-c("TRAIT", "case")
bbj1 <- merge(bbj1, bbj0female); rm(bbj0female)

bbj2<- fread("../../0_data/BBJ/GenMetS4627.binary_trait.male.txt")
bbj2$Cohort<- "BBJ"
bbj2$Sex <- "Men"
bbj0male <- bbj0[, c(1,3)]; colnames(bbj0male)<-c("TRAIT", "case")
bbj2 <- merge(bbj2, bbj0male); rm(bbj0male)

bbj <- rbind(bbj1, bbj2)
bbj$Predictor<- "GenMetS"
bbj$"Major Ancestry" <- "EAS"
bbj$Control <- bbj$case
bbj$Samples <- "matched"
bbj <- bbj %>% rename(Disease = TRAIT, OddsRatio = OR, "p-value" = P,
                      "95%CI lower"= CI_LOW, "95%CI upper"=CI_UP, 
                      Case=case) %>%
  dplyr::select(Cohort, "Major Ancestry", Sex, Disease, Predictor, Samples,
                Case, Control, OddsRatio, "95%CI lower", "95%CI upper",
                "p-value")

bbj <- bbj %>% dplyr::filter(Disease %in% mydisease )


myresult <- rbind(ukb_asn, attract, bbj)
idx<-which(myresult$Disease=="T2D") ; myresult$Disease[idx]<- "Type 2 diabetes"
idx<-which(myresult$Disease=="MI") ; myresult$Disease[idx]<- "Myocardial infarction"
idx<-which(myresult$Disease=="HF") ; myresult$Disease[idx]<- "Heart failure"
idx<-which(myresult$Disease=="stroke") ; myresult$Disease[idx]<- "Stroke"
idx<-which(myresult$Disease=="CAD") ; myresult$Disease[idx]<- "Coronary artery disease"

myresult <- myresult %>% dplyr::filter(!is.na(Disease))
unique(myresult$Disease)
unique(myresult$Cohort)
  

myresult$Disease <- factor(myresult$Disease,
                         levels= c("Type 2 diabetes","Stroke",
                                   "Heart failure", "Coronary artery disease",
                                   "Myocardial infarction", "Hypertension")) 
myresult$Cohort <- factor(myresult$Cohort, levels = c("UKB ASN", "ATTRaCT", "BBJ"))
myresult$Sex <- factor(myresult$Sex, levels = c("Women", "Men"))


myresult1 <- myresult %>%  arrange(Cohort, Sex)   

write.table(myresult, "../../0_data/OddsRatio_result/sup_table_GenMetS_3cohorts.txt", row.names=F, sep="\t", quote=F)


#F8766D
rm(list=ls())
mycolors<- c("#6F99ADFF", "#FFC000FF") #yellow for female

result_title="GenMetS in prediction of cardiometabolic diseases"
result_img1="SubFig_GenMetS_3cohort.tiff"
result_img2="SubFig_GenMetS_3cohort.pdf"

myresult<- fread("../../0_data/OddsRatio_result/sup_table_GenMetS_3cohorts.txt")
myresult$disease <- factor(myresult$disease,
                           levels=rev(c("Type 2 diabetes", "Coronary artery disease", 
                                        "Hypertension","Heart failure",
                                        "Stroke","Myocardial infarction")))

myresult$Cohort <- factor(myresult$Cohort,
                          levels= rev(c("BBJ", "ATTRaCT", "UKB ASN")) )


dplot<- myresult
dplot$OddsRatio <- as.numeric(dplot$OddsRatio)
dplot$cohort <- as.factor(dplot$Cohort)
dplot$high <- as.numeric(dplot$"95%CI upper")
dplot$low <- as.numeric(dplot$"95%CI lower")
dplot$pp <- -log10(dplot$"p-value")
pd <- position_dodge(width=0.5)

myplot <- ggplot(dplot, aes(x = OddsRatio, y = Disease, color = Sex, shape = Sex)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", size = 0.5) +
  geom_point(size = 2, position = pd) +
  geom_errorbarh(aes(xmax = high, xmin = low), height = 0.2, size = 0.5, position = pd) +
  scale_color_manual(values = mycolors) +
  scale_shape_manual(values = c(15, 19)) +
  facet_wrap(~Cohort, ncol = 3) +
  theme_minimal(base_size = 12) + theme_bw() + 
  theme(legend.position = "bottom",
        legend.title = element_text(colour="black", size = 10),
        legend.text = element_text(colour="black",size = 9),
        axis.title = element_text(colour="black",size = 12),
        strip.text.x = element_text(colour="black",size = 10),
        panel.grid.minor = element_blank(),
        panel.border = element_blank()) +
  labs(y = "", x = "OR (95% CI)", title = "") +
  guides(color = guide_legend(title = "Sex"), shape = guide_legend(title = "Sex"))

# Print the plot
print(myplot)

ggsave(result_img1, width =6, height =4, units = "in",dpi = 600, compression = "lzw")
ggsave(result_img2, width =6, height =4, units = "in",dpi = 600)

 
 
 

  

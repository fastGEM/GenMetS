##############
#Pan Hong
#Feb10,2024
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

 
############################UKB WOMEN
myresult_raw1<- fread("../../0_data/OddsRatio_result/OddsRatio_diseases_UKB.txt") %>% data.frame()
myresult_raw2<- fread("../../0_data/OddsRatio_result/OddsRatio_diseases_ATTRaCT.txt") %>% data.frame()
myresult_raw <- rbind(myresult_raw1, myresult_raw2)
mydiseases=c("T2DM", "CAD", "HTN","HF", "Stroke", "MI")
idx<- which(myresult_raw$disease == "t2dm"); myresult_raw$disease[idx]="T2DM"
idx<- which(myresult_raw$disease == "hypertension"); myresult_raw$disease[idx]="HTN"
idx<- which(myresult_raw$disease %in% c("stroke", "STROKE")); myresult_raw$disease[idx]="Stroke"

women1_title="UKB ASN women"
men1_title="UKB ASN men"

women2_title="ATTRaCT women"
men2_title="ATTRaCT men"

idx<- which(myresult_raw$predictor=="PRS_GenMetS")
myresult_raw$predictor[idx]<- "GenMetS"
myresult_raw <- myresult_raw %>% dplyr::filter(!predictor %in% c("PRS_Walree2022", "PRS_Lind2019",
                                                                 "PRS_PGS001108_BasalMetabolicRate", "PRS_PGS000807_T2D"))
unique(myresult_raw$predictor)

myresult_raw$predictor <- factor(myresult_raw$predictor,
                                 levels=rev(c(  "GenMetS",  "PRS_PGS000828_WC1", "PRS_PGS001227_WC2",
                                                "PRS_PGS000661_LDL1", "PRS_PGS000891_LDL2",  "PRS_PGS000660_HDL1",
                                                "PRS_PGS000686_HDL2", "PRS_PGS000659_TG1" ,  "PRS_PGS000699_TG2" ,
                                                "PRS_PGS000306_GLU1", "PRS_PGS000684_GLU2",  "PRS_PGS000685_HbA1c",
                                                "PRS_PGS001133_DBP1" , "PRS_PGS000302_DBP2", "PRS_PGS000301_SBP1", 						   "PRS_PGS001134_SBP2"  )), 
                                 labels=rev(c("GenMetS", "PGS000828_WC#", "PGS001227_WC#",
                                              "PGS000661_LDL", "PGS000891_LDL", "PGS000660_HDL",
                                              "PGS000686_HDL#", "PGS000659_TG" , "PGS000699_TG#" ,
                                              "PGS000306_GLU#",  "PGS000684_GLU#", "PGS000685_HbA1c#",
                                              "PGS001133_DBP#" , "PGS000302_DBP#",  "PGS000301_SBP#", 
					     "PGS001134_SBP#" )))

unique(myresult_raw$disease)
myresult_raw<- myresult_raw %>% dplyr::filter(disease %in% mydiseases)
myresult_raw$disease <- factor(myresult_raw$disease,
                               levels=c("T2DM", "CAD", "HTN","HF", "Stroke", "MI"),
                               labels=c("Type 2 diabetes", "Coronary artery disease", "Hypertension",
                                        "Heart failure","Stroke", "Myocardial infarction"))
write.table(myresult_raw, "SubFig4_data.txt", row.names=F, sep="\t", quote=F)
OddsRatio_plot <- function(myinputdata, mytitle)
{
dplot <- myinputdata
mycolors<- c(rep("black", 15), "#F8766D")
dplot$OddsRatio <- as.numeric(dplot$OddsRatio)
dplot$high <- as.numeric(dplot$high)
dplot$low <- as.numeric(dplot$low)
pd <- position_dodge(width=0.5)
myplot<- ggplot(dplot, aes(x =OddsRatio, y = predictor , color=predictor )) +
  geom_vline(aes(xintercept = 1), linewidth = .25, linetype = "dashed") +
  geom_point(shape = 18, size = 2, position=pd) +  
  geom_errorbarh(aes(xmax = high, xmin = low), linewidth =0.5, position=pd, height = .4) +
  #xlim(c(0, 2)) +
  scale_color_manual(values = mycolors) +
  facet_wrap(~disease, ncol = 3) +
  theme_bw()+theme(axis.text=element_text(colour="black",size=8),
                   strip.text.x = element_text(colour="black", size = 8), 
                   legend.text=element_text(colour="black", size=8)) +
  theme(panel.grid.minor = element_blank()) +
  theme(legend.position = "")+
  ylab("") + xlab("OR (95% CI)")  +
  ggtitle(mytitle)

#print(myplot)
return(myplot)
}

my_UKB_ASN_women <- myresult_raw %>% dplyr::filter(cohort=="UKB ASN women")
myplot1 <- OddsRatio_plot(my_UKB_ASN_women, women1_title); rm(my_UKB_ASN_women)
my_UKB_ASN_men <- myresult_raw %>% dplyr::filter(cohort=="UKB ASN men")
myplot2 <- OddsRatio_plot(my_UKB_ASN_men, men1_title); rm(my_UKB_ASN_men)

my_ATTRACT_women <- myresult_raw %>% dplyr::filter(cohort=="ATTRaCT women")
myplot3 <- OddsRatio_plot(my_ATTRACT_women, women2_title); rm(my_ATTRACT_women)
my_ATTRACT_men <- myresult_raw %>% dplyr::filter(cohort=="ATTRaCT men")
myplot4 <- OddsRatio_plot(my_ATTRACT_men, men2_title); rm(my_ATTRACT_men)


library(ggpubr)
combined_plot1 <- ggarrange(myplot1,myplot2, nrow=2, ncol=1,
                           labels = c("a", "b"), 
                           font.label=list(color="black",size=24)) 
print(combined_plot1)
ggsave("SubFig_benchmark_15PGS_UKB_ASN_men_women.tiff", combined_plot1, width =8, height =12, units = "in",dpi = 300, compression = "lzw")
ggsave("SubFig_benchmark_15PGS_UKB_ASN_men_women.pdf", combined_plot1, width =8, height =12, units = "in" )


combined_plot2 <- ggarrange(myplot3,myplot4, nrow=2, ncol=1,
                            labels = c("a", "b"), 
                            font.label=list(color="black",size=24)) 
print(combined_plot2)
ggsave("SubFig_benchmark_15PGS_ATTRACT_men_women.tiff", combined_plot2, width =8, height =12, units = "in",dpi = 300, compression = "lzw")
ggsave("SubFig_benchmark_15PGS_ATTRACT_men_women.pdf", combined_plot2, width =8, height =12, units = "in" )


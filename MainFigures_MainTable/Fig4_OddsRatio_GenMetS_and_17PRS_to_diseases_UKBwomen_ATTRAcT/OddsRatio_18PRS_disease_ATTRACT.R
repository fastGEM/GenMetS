######################
#Pan Hong

#last update Nov13,2023
#test 17PRS + GenMetS = 18PRS
#in association cardiometabolic diseases
#in ATTRaCT
#PRS_PGS000306_GLU1 is null in ATTRACT
#Nov22,2023
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


myresult_raw<- fread("OddsRatio_18PRS_to_diseases_ATTRACT.txt") %>% data.frame()

idx<- which(myresult_raw$predictor=="PPRS_PGS000301_SBP1")
myresult_raw$predictor[idx]<- "PRS_PGS000301_SBP1"

idx<- which(myresult_raw$predictor=="PredMetS")
myresult_raw$predictor[idx]<- "GenMetS"

myresult_raw <- myresult_raw %>% dplyr::filter(!is.na(predictor))
myresult_raw <- myresult_raw %>% dplyr::filter(disease != "DEPRESSION")
idx<-which(myresult_raw$disease == "STROKE")
myresult_raw$disease[idx]<- "Stroke"

unique(myresult_raw$disease)

myresult_raw$predictor <- factor(myresult_raw$predictor,
                                 levels=rev(c(  "GenMetS",                         
                                                "PRS_PGS000828_WC1",
                                                "PRS_PGS001227_WC2",
                                                "PRS_PGS000661_LDL1",              
                                                "PRS_PGS000891_LDL2",
                                                "PRS_PGS000660_HDL1",
                                                "PRS_PGS000686_HDL2",
                                                "PRS_PGS000659_TG1" ,              
                                                "PRS_PGS000699_TG2" ,
                                                "PRS_PGS000684_GLU2",
                                                "PRS_PGS001133_DBP" ,              
                                                "PRS_PGS000301_SBP1",
                                                "PRS_PGS001134_SBP2",              
                                                "PRS_PGS000807_T2D" ,
                                                "PRS_PGS001108_BasalMetabolicRate",
                                                "PRS_PGS000685_HbA1c", 
                                                "PRS_PGS000958_BP"         
                                 )),
                                 labels=rev(c("GenMetS",
                                              "PGS000828_WC#",
                                              "PGS001227_WC#",
                                              "PGS000661_LDL",              
                                              "PGS000891_LDL",
                                              "PGS000660_HDL",
                                              "PGS000686_HDL#",
                                              "PGS000659_TG" ,              
                                              "PGS000699_TG#" ,
                                              "PGS000684_GLU#",
                                              "PGS001133_DBP#" ,              
                                              "PGS000301_SBP#",
                                              "PGS001134_SBP#",              
                                              "PGS000807_T2D" ,
                                              "PGS001108_BasalMetabolicRate#",
                                              "PGS000685_HbA1c#", 
                                              "PGS000958_BP#"    
                                 )))


myresult_raw$disease <- factor(myresult_raw$disease,
levels=c("T2DM", "CAD", "HF","HTN", "Stroke", "MI", "PVD", "COPD"))

#F8766D
mycolors<- c(rep("black", 16), "#F8766D")

dplot<- myresult_raw
dplot$boxLabels <- paste0(dplot$disease, ":", dplot$predictor)
dplot$OddsRatio <- as.numeric(dplot$OddsRatio)
dplot$high <- as.numeric(dplot$high)
dplot$low <- as.numeric(dplot$low)

pd <- position_dodge(width=0.5)


ggplot(dplot, aes(x =OddsRatio, y = predictor , color=predictor )) +
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") +
  geom_point(shape = 18, size = 2, position=pd) +  
  geom_errorbarh(aes(xmax = high, xmin = low), linewidth =1, position=pd, height = .2) +
  scale_color_manual(values = mycolors) +
  facet_wrap(~disease, ncol = 4) +
  theme_bw()+theme(axis.text=element_text(colour="black",size=8),
                   legend.text=element_text(colour="black", size=8)) +
  theme(panel.grid.minor = element_blank()) +
  theme(legend.position = "")+
  ylab("") +
  xlab("OR (95% CI)")  # ggtitle("PredMetS vs ObsMetS in case/control")

ggsave("OddsRatio_18PRS_in_ATTRACT_disease_case_control.tiff", 
       width =8, height =6, units = "in",dpi = 300, compression = "lzw")

ggsave("OddsRatio_18PRS_in_ATTRACT_disease_case_control.pdf", 
       width =8, height =6, units = "in")



 



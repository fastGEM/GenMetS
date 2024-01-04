#############################################
#Pan Hong
#
#18PRS=GenMetS+17PRSfromliterature
#benchmarking GenMetS and 17 PRSs in explaining
#Observed MetS by R2. 
#Nov 22, 2023
##############################################
rm(list=ls())
library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
theme_set(theme_pubr())

blues <- brewer.pal(9, "Blues")
blue_range <- colorRampPalette(blues)

dplot <- fread("data_GenMetS_to_ObsMetS_for_plot.txt")

#change the display order
dplot$cohort <- factor(dplot$cohort,
                     levels=c( "S-PRESTO,GUSTO", "UKB Asian", "UKB European",
                               "ATTRaCT", "GUSTOchildren"))


dplot$predictor <- factor(dplot$predictor,
                                 levels=c(  "GenMetS",                         
                                                "PRS_PGS000828_WC1",
                                                "PRS_PGS001227_WC2",
                                                "PRS_PGS000661_LDL1",              
                                                "PRS_PGS000891_LDL2",
                                                "PRS_PGS000660_HDL1",
                                                "PRS_PGS000686_HDL2",
                                                "PRS_PGS000659_TG1" ,              
                                                "PRS_PGS000699_TG2" ,
                                                "PRS_PGS000306_GLU1",              
                                                "PRS_PGS000684_GLU2",
                                                "PRS_PGS001133_DBP" ,              
                                                "PRS_PGS000301_SBP1",
                                                "PRS_PGS001134_SBP2",              
                                                "PRS_PGS000807_T2D" ,
                                                "PRS_PGS001108_BasalMetabolicRate",
                                                "PRS_PGS000685_HbA1c", 
                                                "PRS_PGS000958_BP"         
                                 ),
                                 labels=c("GenMetS",
                                              "PGS000828_WC#",
                                              "PGS001227_WC#",
                                              "PGS000661_LDL",              
                                              "PGS000891_LDL",
                                              "PGS000660_HDL",
                                              "PGS000686_HDL#",
                                              "PGS000659_TG" ,              
                                              "PGS000699_TG#" ,
                                              "PGS000306_GLU#",              
                                              "PGS000684_GLU#",
                                              "PGS001133_DBP#" ,              
                                              "PGS000301_SBP#",
                                              "PGS001134_SBP#",              
                                              "PGS000807_T2D" ,
                                              "PGS001108_Basal#",
                                              "PGS000685_HbA1c#", 
                                              "PGS000958_BP#"    
                                 ))



ggplot(data=dplot, aes(x=predictor, y=R2, fill=predictor)) +
  geom_bar(stat="identity", color="black", position=position_dodge(), width=0.8)+
  geom_errorbar(aes(ymin=R2, ymax=high), width=.2,position=position_dodge(.9)) +
  #scale_fill_manual(values = blue_range(3)) +
  facet_wrap(vars(sex,cohort), ncol=5 )+
  theme_bw() + 
  labs(y=expression ("Coefficient determination of Observed MetS "~(R^2)))  +
  theme(
    legend.position = "bottom",
    legend.key=element_rect(colour="black"),
    legend.key.size = unit(0.4, 'cm'),
    strip.text = element_text(size = 8, face = "bold",hjust = 0.5, vjust = 0.9),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(0,4), "cm") 
  )
 


pdfv1="barplot_18PRS_to_ObsMetS.pdf"
imgv1="barplot_18PRS_to_ObsMetS.tiff"
ggsave(imgv1, width = 8.4, height =5.34,units = "in", dpi = 300, compression = "lzw")
ggsave(pdfv1, width = 8.4, height =5.34,units = "in")

  

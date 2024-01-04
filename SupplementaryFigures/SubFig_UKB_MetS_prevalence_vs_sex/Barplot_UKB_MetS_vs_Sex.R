##################
#Pan Hong
#
#latest update: July 17, 2023
##################

rm(list=ls())

library(data.table)
library(dplyr) 
library(RColorBrewer)
library(ggpubr)
library(ggplot2)
 

myresult1 <- fread("UKB_MetS_prevalence_vs_sex.txt") %>% data.frame()

  ggplot(myresult1, aes(x=Ethnicity, y=freq, fill=Sex))+  
  geom_bar(position = position_dodge(), stat = "identity",color='white', width=0.8)+ coord_flip()+
  #scale_y_continuous(labels = scales::percent))
  geom_text(aes(label=paste0(round(freq*100,2),"%")), 
            position = position_dodge( 0.6), color="black", size=3) +  
  xlab("") + ylab("MetS Severity Rate") + theme_bw() +ylim(0, 0.8 ) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        text=element_text(family="serif", colour="black", size=12),
        axis.title=element_text(family="serif", colour="black", size=12),
        axis.text =element_text(family="serif", color="black", size=12),
        legend.title = element_blank(),
        legend.text=element_text(family="serif", color="black", size=12))
  
  
  ggsave("UKB_MetS_vs_Sex.tiff", width = 7.635,height =4.7175,units = "in")

  ggsave("UKB_MetS_vs_Sex.pdf", width = 7.635,height =4.7175,units = "in")
  
  
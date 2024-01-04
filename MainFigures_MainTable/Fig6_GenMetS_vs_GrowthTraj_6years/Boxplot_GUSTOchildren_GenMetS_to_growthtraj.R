############################################
#Pan Hong
#boxplot for GenMetS vs GrowthTraj for
#GUSTO children at 6 years

#revisited on Oct 20,2023
#revisited on Nov 29, 2023
############################################

rm(list=ls())
library(dplyr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(ggpubr)

my_wd <- "/Users/TEHAL/OneDrive - A STAR/Misc/MetS-PH/" 
#my_wd <- "C:/Users/tehal/OneDrive - A STAR/Misc/MetS-PH/"
setwd(my_wd)

##########################
myfile="GenMetS_Traj_children6years.txt"
mydata<- fread(myfile) %>% data.frame()
dim(mydata)
mydata[1:6, 1:6]

mydata$Traj <- factor(mydata$Traj_class,
                      levels=c("1", "2", "3", "5", "4"),
                      labels=c("Normal Low", "Normal", "Normal High", 
                               "Late Accelerated", "Early Accelerated"))

pData1<- data.frame(PredMetS=scale(mydata$predicted_mets), Traj=mydata$Traj) %>% 
	dplyr::filter(!is.na(Traj))

table(pData1$Traj)
table(mydata$Traj)

mycolorcode <-  c("#C0C0C0", "#C0C0C0","#C0C0C0", "#808080", "#808080") 
mycolorcode <- adjustcolor(mycolorcode,alpha.f=0.5)

ggboxplot(pData1, x = "Traj", y = "PredMetS",
               color = "black", palette =mycolorcode, fill = "Traj",
               add = "jitter",
               xlab="",
               ylab="zscore of GenMetS",
               width = 0.8,
               size = 0.3,
               #title="Child Growth Trajectory Patterns",
               legend="none",
               font.label = list(size = 10, face = "bold", family="serif"),
               orientation = "horizontal",
          add.params = list(size = 0.8, alpha = 1)
) +
  theme(
    plot.title = element_text(size = 14, face="bold", hjust = 0.5),
    axis.title.x = element_text(size = 14, family = "serif"),
    axis.title.y = element_text(size = 11, family = "serif"),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 13, family = "serif"),
    axis.ticks.y=element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_line(color = 'black', size = 0.8)
  )+
  stat_compare_means( comparisons =list(c("Normal Low", "Normal"), c("Normal", "Normal High"),
                                        c("Normal","Late Accelerated"), 
                                        c("Normal", "Early Accelerated") )) +
  stat_compare_means( method="kruskal.test", label.y = 2.8, label.x = 0.65, size=4) 


ggsave("Growth_Trajectory_class_vs_MetS.tiff",
       width =6.3, height =5, units = "in",dpi = 600, compression = "lzw")

ggsave("Growth_Trajectory_class_vs_MetS.pdf",width=6.3, height=5)
 


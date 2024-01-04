################################
#Pan Hong
#May 21, 2023
################################


library(dplyr)
library(data.table)
library(ggplot2)

trait="multimorbidity"
feature_importance_fig=paste0(trait, "_M2model_SHAP_feature_importance.tiff")

myresult <- fread("Feature_importance_multimorbidity_m2model.txt") 
idx<-which(myresult$feature== "SocialEconomicalStatus" )
myresult$feature[idx]="SocioeconomicStatus"
idx<-which(myresult$feature== "PredMetS" )
myresult$feature[idx]="GenMetS"


dplot<- myresult
dplot$boxLabels = paste0( dplot$feature)
dplot$mean <- as.numeric(dplot$mean)
dplot$sd <- as.numeric(dplot$sd)
dplot$high <- as.numeric(dplot$mean) 
 

mytitle=paste0("Feature importance in prediction of ", trait)
dplot$mean <- round(dplot$mean, digits=3)


library(ggplot2)
#library(flextable)
library(grid)
library(cowplot)
library(tidyverse)

myresult <- myresult %>% arrange(desc(mean))
mydf <- tibble(features = myresult$feature,
               importance = myresult$mean)


p1<- ggplot(dplot, aes(x=reorder(feature, mean), y=mean, color=feature )) + 
  geom_bar(stat="identity", width=0.5,fill="steelblue",
           color="steelblue", position=position_dodge()) +
  #geom_text(aes(label=mean) hjust=-(dplot$mean+dplot$sd), color="black", size=4)+
  #geom_text(aes(label=mean), hjust=-dplot$high, color="black", size=3.5)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean + sd), 
                width=.2,position=position_dodge(.9), color="blue") +  
  # Set origin of axes to zero                            
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip() +
  labs(title=mytitle, x="", y = "mean of abs(SHAP)")+
  theme_classic() +
  theme(axis.text = element_text(size=12,color="black")) +
  theme(legend.position="none")

print(p1)

 

ggsave(feature_importance_fig,  
       width =6.2, height = 4.14, units = "in",dpi = 300, compression = "lzw")

 
 

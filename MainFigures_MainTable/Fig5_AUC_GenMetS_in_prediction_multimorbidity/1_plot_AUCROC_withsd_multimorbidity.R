################################
#Pan Hong
#May 2, 2023
#last update May 13, 2023
#exclude ObsMetS (m1, m2)
################################

rm(list=ls())
library(dplyr)
library(data.table)
library(ggplot2)


trait="multimorbidity"
model_ROCAUC_fig=paste0("GenMetS_prediction_AUC_for_", trait, ".tiff")
mytitle=paste0("ROC for prediction of ", trait, " case and control")

dplot<- fread("GenMetS_prediction_multimorbidity.txt")
dplot <-dplot %>% data.frame() %>%
  mutate(model = case_when(
                         model == 'm2' ~ 'GenMetS+covariates, AUC(sd)=0.69(0.03)',
                         model == 'm4' ~ 'covariates, AUC(sd)=0.66(0.02)'))

 

dplot$model <- factor(dplot$model, levels = c(
                                              "GenMetS+covariates, AUC(sd)=0.69(0.03)",
                                              "covariates, AUC(sd)=0.66(0.02)"))

ggplot(data = dplot, aes(x = fpr, model=model)) + 
  geom_line(aes(y = mean, color = model), size = 1) + 
  geom_ribbon(aes(y =mean , ymin = mean - sd, ymax = mean + sd, fill = model), alpha = .2) +
  xlab("FPR") + ylab("TPR")+ 
  ggtitle(mytitle) +
  theme_bw() +  
  #theme(legend.key = element_blank()) + 
  #theme(plot.margin=unit(c(1,3,1,1),"cm"))+
  #theme(legend.position = c(0.6,0.2), legend.direction = "vertical") +
  theme(legend.title = element_blank()) + 
  #theme_classic() +
  theme(axis.text = element_text(size=12,color="black")) +
  theme(legend.text=element_text(size=12,color="black")) +
  theme(legend.position = c(0.6,0.2), legend.direction = "vertical") 

 
ggsave(model_ROCAUC_fig,  
       width =6.2, height = 4.14, units = "in",dpi = 600, compression = "lzw")

 
 

##############
#Pan Hong
#May 20, 2024
#############

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

#forest plot for the odds ratio

colors=c("#0072B2", "#D55E00")


odds_plot <- function(dplot) {
  pd <- position_dodge(width=0.5)
  mycolors<- c(rep("black",3),  "#F8766D") #yellow for female
  
  pp<- ggplot(dplot, aes(x =OddsRatio, y = predictor , color=predictor )) +
    geom_vline(aes(xintercept = 1), linewidth = .25, linetype="dashed") +  
    geom_point(size = 3, shape=19, position=pd) +  
    geom_errorbarh(aes(xmax = high, xmin = low), linewidth =0.5, position=pd, height = .2) +
    scale_color_manual(values = mycolors) +
    facet_wrap(~disease, ncol=3) +
    theme_bw()+
    theme(axis.text=element_text(colour="black",size=8),
                   strip.text.x = element_text(colour="black", size = 8), 
                   legend.text=element_text(colour="black", size=8)) +
    theme(panel.grid.major = element_blank(),  # Remove major gridlines
          panel.grid.minor = element_blank(),  # Remove minor gridlines
          panel.background = element_rect(fill = "white", color = NA),  # Set background to white
          plot.background = element_rect(fill = "white", color = NA),
          panel.spacing.y = unit(0.1, "lines")) +    # Set plot background to white
  
    theme(legend.position = "bottom")+
    ylab("") + xlab("OR (95% CI)")  + ggtitle("") +
    guides(color = guide_none())  # This line removes the color guide
 return(pp)
}

###############################this_data
study_data <- fread("result_META_benchmark.csv") %>% data.frame() %>%
  rename(OddsRatio = Combined_OR,
         high = Upper_CI,
         low  = Lower_CI,
         disease = Disease) %>%
  dplyr::filter(!is.na(disease)) %>%
  dplyr::filter(Cohort=="META") %>%
  dplyr::select(-Weight, -Combined_N)  

 
#t2dm
t2dm<- study_data %>% dplyr::filter(disease=="Type 2 diabetes")  
t2dm$predictor <- factor(t2dm$predictor,
levels=rev(c("GenMetS","Walree2022_MetS","Lind2019_MetS","Suzuki2024_T2DM")))
p1<- odds_plot(t2dm)
print(p1)
#rm(t2dm, dplot1)

################HF
HF<- study_data %>% dplyr::filter(disease=="Heart failure") 
HF$predictor <- factor(HF$predictor,
  levels=rev(c("GenMetS","Walree2022_MetS","Lind2019_MetS","Levin2022_HF")))
p2<-odds_plot(HF)
print(p2)
#rm(HF, dplot1)

##########Hyertension
HTN<- study_data %>% dplyr::filter(disease=="Hypertension")  
HTN$predictor <- factor(HTN$predictor,
                       levels=rev(c("GenMetS","Walree2022_MetS","Lind2019_MetS","Sinnott2021_HTN")))
p3<-odds_plot(HTN)
print(p3)
 
##########4MI
MI<- study_data %>% dplyr::filter(disease=="Myocardial infarction")  
MI$predictor <- factor(MI$predictor,
                        levels=rev(c("GenMetS","Walree2022_MetS","Lind2019_MetS","Hartiala2021_MI")))

p4<- odds_plot(MI)
print(p4)
 


##########5Stroke
stroke<- study_data %>% dplyr::filter(disease=="Stroke") 

stroke$predictor <- factor(stroke$predictor,
                         levels=rev(c("GenMetS",
                                      "Walree2022_MetS", 
                                      "Lind2019_MetS", 
                                      "Mishra2022_stroke")))

p5<- odds_plot(stroke)
print(p5) 

##########6CAD
CAD<- study_data %>% dplyr::filter(disease=="Coronary artery disease") %>%
  dplyr::filter(predictor %in% c("GenMetS", "Walree2022_MetS", "Lind2019_MetS", "Aragam2022_CAD"))
CAD$predictor <- factor(CAD$predictor,
                             levels=rev(c("GenMetS",
                                          "Walree2022_MetS", 
                                          "Lind2019_MetS", 
                                          "Aragam2022_CAD")))


p6<- odds_plot(CAD)
print(p6)

 
combined_plot <- ggarrange(p1, p2, p3, p4, p5, p6,
                           nrow = 2, ncol = 3,
                           common.legend = TRUE, legend = "bottom")
print(combined_plot)

 

ggsave("Fig4b.tiff",  
       width =8.2, height =4.2, units = "in",dpi = 300, compression = "lzw")

ggsave("Fig4b.pdf",  
       width =8.2, height =4.2, units = "in",dpi = 300)  


 
 
 
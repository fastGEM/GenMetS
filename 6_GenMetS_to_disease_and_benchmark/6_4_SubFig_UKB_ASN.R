##############
#Pan Hong
#May 20, 2024
#############
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
rm(list=ls())

#forest plot for the odds ratio

odds_plot <- function(dplot) {
  pd <- position_dodge(width=0.5)
  mycolors<- c(rep("black",6), "#6F99ADFF", "#FFC000FF") #yellow for female
  
  pp<- ggplot(dplot, aes(x =OddsRatio, y = predictor , color=pred_cohort, shape=cohort  )) +
  geom_vline(aes(xintercept = 1), linewidth = .25, linetype = "dashed") +
  geom_point(size = 2, position=pd) +  
  geom_errorbarh(aes(xmax = high, xmin = low), linewidth =0.5, position=pd, height = .4) +
  scale_color_manual(values = mycolors) +
  scale_shape_manual(values = c(15, 19)) +
  facet_wrap(~disease, ncol = 3) +
  theme_bw()+theme(axis.text=element_text(colour="black",size=8),
                   strip.text.x = element_text(colour="black", size = 8), 
                   legend.text=element_text(colour="black", size=8)) +
  theme(panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom")+
  ylab("") +
  xlab("OR (95% CI)")  +
  ggtitle("") +
  guides(color = guide_none())  # This line removes the color guide
 return(pp)
}

###############################UKB_ASN
UKB_ASN<- fread("../../0_data/OddsRatio_result/OR_benchmark_diseasePGS_in_UKB_ASN.txt") %>% data.frame() %>%
  dplyr::filter(!is.na(disease))
unique(UKB_ASN$predictor)
unique(UKB_ASN$cohort)

#rename predictor as AuthorYear_disease
idx<-which(UKB_ASN$predictor=="HF_meta"); UKB_ASN$predictor[idx]<- "Levin2022_HF"
idx<-which(UKB_ASN$predictor=="MI_EUR"); UKB_ASN$predictor[idx]<- "Hartiala2021_MI"
idx<-which(UKB_ASN$predictor=="T2DM_EAS"); UKB_ASN$predictor[idx]<- "Suzuki2024_T2DM"
idx<-which(UKB_ASN$predictor=="CAD_EUR"); UKB_ASN$predictor[idx]<- "Aragam2022_CAD"
idx<-which(UKB_ASN$predictor=="Hypertension_EUR"); UKB_ASN$predictor[idx]<- "Sinnott2021_HTN"
idx<-which(UKB_ASN$predictor=="Stroke_EAS"); UKB_ASN$predictor[idx]<- "Mishra2022_stroke"
unique(UKB_ASN$predictor)

UKB_ASN$pred_cohort<- paste0(UKB_ASN$predictor, ":", UKB_ASN$cohort)

#t2dm
unique(UKB_ASN$disease)
t2dm<- UKB_ASN %>% dplyr::filter(disease=="t2dm") %>%
  dplyr::filter(predictor %in% c("GenMetS", "Walree2022_MetS", "Lind2019_MetS", "Suzuki2024_T2DM")) %>%
  data.frame()


t2dm$disease="Type 2 diabetes"
t2dm$predictor <- factor(t2dm$predictor,
levels=rev(c("GenMetS","Walree2022_MetS","Lind2019_MetS","Suzuki2024_T2DM")))

t2dm$pred_cohort <- factor(t2dm$pred_cohort,
levels=rev(c("GenMetS:UKB ASN women","GenMetS:UKB ASN men",
             "Walree2022_MetS:UKB ASN women", "Walree2022_MetS:UKB ASN men",
             "Lind2019_MetS:UKB ASN women", "Lind2019_MetS:UKB ASN men",
             "Suzuki2024_T2DM:UKB ASN women", "Suzuki2024_T2DM:UKB ASN men")))

dplot1<- t2dm
dplot1$OddsRatio <- as.numeric(dplot1$OddsRatio)
dplot1$high <- as.numeric(dplot1$high)
dplot1$low <- as.numeric(dplot1$low)
dplot1$cohort <- as.factor(dplot1$cohort)
dplot1$pred_cohort<- as.factor(dplot1$pred_cohort)

p1<- odds_plot(dplot1)
print(p1)
rm(dplot1)

################HF
HF<- UKB_ASN %>% dplyr::filter(disease=="HF") %>%
  dplyr::filter(predictor %in% c("GenMetS", "Walree2022_MetS", "Lind2019_MetS", "Levin2022_HF")) %>%
  data.frame()

HF$disease="Heart Failure"
HF$predictor <- factor(HF$predictor,
  levels=rev(c("GenMetS","Walree2022_MetS","Lind2019_MetS","Levin2022_HF")))

HF$pred_cohort <- factor(HF$pred_cohort,
  levels=rev(c("GenMetS:UKB ASN women","GenMetS:UKB ASN men",
               "Walree2022_MetS:UKB ASN women", "Walree2022_MetS:UKB ASN men",
               "Lind2019_MetS:UKB ASN women", "Lind2019_MetS:UKB ASN men",
               "Levin2022_HF:UKB ASN women", "Levin2022_HF:UKB ASN men")))



dplot1<- HF
dplot1$OddsRatio <- as.numeric(dplot1$OddsRatio)
dplot1$high <- as.numeric(dplot1$high)
dplot1$low <- as.numeric(dplot1$low)
dplot1$cohort <- as.factor(dplot1$cohort)
dplot1$pred_cohort <- as.factor(dplot1$pred_cohort)

p2<-odds_plot(dplot1)
print(p2)
rm(dplot1)

##########Hyertension
HTN<- UKB_ASN %>% dplyr::filter(disease=="hypertension") %>%
  dplyr::filter(predictor %in% c("GenMetS", "Walree2022_MetS", "Lind2019_MetS", "Sinnott2021_HTN"))

HTN$disease="Hypertension"
HTN$predictor <- factor(HTN$predictor,
                       levels=rev(c("GenMetS","Walree2022_MetS","Lind2019_MetS","Sinnott2021_HTN")))

HTN$pred_cohort <- factor(HTN$pred_cohort,
                         levels=rev(c("GenMetS:UKB ASN women","GenMetS:UKB ASN men",
                                      "Walree2022_MetS:UKB ASN women", "Walree2022_MetS:UKB ASN men",
                                      "Lind2019_MetS:UKB ASN women", "Lind2019_MetS:UKB ASN men",
                                      "Sinnott2021_HTN:UKB ASN women", "Sinnott2021_HTN:UKB ASN men")))

dplot1<- HTN
dplot1$OddsRatio <- as.numeric(dplot1$OddsRatio)
dplot1$high <- as.numeric(dplot1$high)
dplot1$low <- as.numeric(dplot1$low)
dplot1$cohort <- as.factor(dplot1$cohort)
dplot1$pred_cohort <- as.factor(dplot1$pred_cohort)

p3<-odds_plot(dplot1)
print(p3)
rm(dplot1)

##########4MI
MI<- UKB_ASN %>% dplyr::filter(disease=="MI") %>%
  dplyr::filter(predictor %in% c("GenMetS", "Walree2022_MetS", "Lind2019_MetS", "Hartiala2021_MI"))

MI$disease="Myocardial infarction"
MI$predictor <- factor(MI$predictor,
                        levels=rev(c("GenMetS","Walree2022_MetS","Lind2019_MetS","Hartiala2021_MI")))

MI$pred_cohort <- factor(MI$pred_cohort,
                          levels=rev(c("GenMetS:UKB ASN women","GenMetS:UKB ASN men",
                                       "Walree2022_MetS:UKB ASN women", "Walree2022_MetS:UKB ASN men",
                                       "Lind2019_MetS:UKB ASN women", "Lind2019_MetS:UKB ASN men",
                                       "Hartiala2021_MI:UKB ASN women", "Hartiala2021_MI:UKB ASN men")))

                       
dplot1<- MI
dplot1$OddsRatio <- as.numeric(dplot1$OddsRatio)
dplot1$high <- as.numeric(dplot1$high)
dplot1$low <- as.numeric(dplot1$low)
dplot1$cohort <- as.factor(dplot1$cohort)
dplot1$pred_cohort <- as.factor(dplot1$pred_cohort)


p4<- odds_plot(dplot1)
print(p4)
rm(dplot1)


##########5Stroke
stroke<- UKB_ASN %>% dplyr::filter(disease=="stroke") %>%
  dplyr::filter(predictor %in% c("GenMetS", "Walree2022_MetS", "Lind2019_MetS", "Mishra2022_stroke"))

stroke$disease="Stroke"
stroke$predictor <- factor(stroke$predictor,
                       levels=rev(c("GenMetS","Walree2022_MetS","Lind2019_MetS","Mishra2022_stroke")))

stroke$pred_cohort <- factor(stroke$pred_cohort,
                         levels=rev(c("GenMetS:UKB ASN women","GenMetS:UKB ASN men",
                                      "Walree2022_MetS:UKB ASN women", "Walree2022_MetS:UKB ASN men",
                                      "Lind2019_MetS:UKB ASN women", "Lind2019_MetS:UKB ASN men",
                                      "Mishra2022_stroke:UKB ASN women", "Mishra2022_stroke:UKB ASN men")))


dplot1<- stroke
dplot1$OddsRatio <- as.numeric(dplot1$OddsRatio)
dplot1$high <- as.numeric(dplot1$high)
dplot1$low <- as.numeric(dplot1$low)
dplot1$cohort <- as.factor(dplot1$cohort)
dplot1$pred_cohort <- as.factor(dplot1$pred_cohort)

p5<- odds_plot(dplot1)
print(p5)
rm(dplot1)

##########6CAD
CAD<- UKB_ASN %>% dplyr::filter(disease=="CAD") %>%
  dplyr::filter(predictor %in% c("GenMetS", "Walree2022_MetS", "Lind2019_MetS", "Aragam2022_CAD"))

CAD$disease="Coronary artery disease"
CAD$predictor <- factor(CAD$predictor,
                           levels=rev(c("GenMetS","Walree2022_MetS","Lind2019_MetS","Aragam2022_CAD")))

CAD$pred_cohort <- factor(CAD$pred_cohort,
                             levels=rev(c("GenMetS:UKB ASN women","GenMetS:UKB ASN men",
                                          "Walree2022_MetS:UKB ASN women", "Walree2022_MetS:UKB ASN men",
                                          "Lind2019_MetS:UKB ASN women", "Lind2019_MetS:UKB ASN men",
                                          "Aragam2022_CAD:UKB ASN women", "Aragam2022_CAD:UKB ASN men")))


dplot1<- CAD
dplot1$OddsRatio <- as.numeric(dplot1$OddsRatio)
dplot1$high <- as.numeric(dplot1$high)
dplot1$low <- as.numeric(dplot1$low)
dplot1$cohort <- as.factor(dplot1$cohort)
dplot1$pred_cohort <- as.factor(dplot1$pred_cohort)

p6<- odds_plot(dplot1)
print(p6)
rm(dplot1)

 
combined_plot <- ggarrange(p1, p2, p3, p4, p5, p6,
                           nrow = 2, ncol = 3,
                           common.legend = TRUE, legend = "bottom")
print(combined_plot)

final_plot <- annotate_figure(combined_plot,
                              top = text_grob("UKB ASN", size = 14, face = "bold"))

# Print the final plot
print(final_plot)

ggsave("SubFig_Benchmark_UKB_ASN.tiff",  
       width =8, height =5, units = "in",dpi = 600, compression = "lzw")

ggsave("SubFig_Benchmark_UKB_ASN.pdf",  
       width =8, height =5, units = "in",dpi = 600)  




#
mydata <- rbind(t2dm,HF, HTN, MI, stroke, CAD)
idx<-which(mydata$disease == "Heart Failure")
mydata$disease[idx]<- "Heart failure"
mydata$Sex <- "Women"
idx <- which(mydata$cohort == "UKB ASN men")
mydata$Sex[idx] <- "Men"
mydata$Cohort <- "UKB ASN"
mydata$"Major Ancestry" <- "Asian"
colnames(mydata)
table(mydata$disease)

 
mydata$disease <- factor(mydata$disease,
                         levels= c("Type 2 diabetes","Stroke",
                                   "Heart failure", "Coronary artery disease",
                                   "Myocardial infarction", "Hypertension")) 


mydata1 <- mydata %>% 
  rename(Disease = disease, Samples = samples, Predictor = predictor,
         Case = case, Control = control, 
         "95%CI lower" = low, "95%CI upper" = high,
         "p-value" = "p.value",) %>%
  dplyr::select(Cohort,	"Major Ancestry", Sex,	Disease, Predictor, Samples,	Case,	Control,
                OddsRatio,	"95%CI lower"	,"95%CI upper",	"p-value") %>% 
  arrange(factor(Sex, levels = c("Women", "Men")))  %>%
  arrange(Disease)

mydata1

mydata2<- mydata1 %>% dplyr::filter(Predictor=="GenMetS")  


mydata2

write.table(mydata1, "sup_table_UKB_ASN_benchmark.txt", row.names=F, sep="\t", quote=F) 
write.table(mydata2, "sup_table_UKB_ASN_GenMetS.txt", row.names=F, sep="\t", quote=F) 




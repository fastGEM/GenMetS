############################################
#Pan Hong
#Mar21,2023
#estimated MetS vs child growth trajectory
#
#Five distinct trajectory classes were identified from previous study. 
#The three stable trajectories, 
#stable normal, low, high, stayed within the normal range of BMI
#The early-acceleration (EA) trajectory showed BMI acceleration 
#immediately after birth and crossed obesity threshold after 2years of age. 
#The late-acceleration [LA] trajectory was close the Normal 
#trajectory in the first year, started accelerating after age 1year 
#and was close to the obesity threshold by age 6years. 

#revisited on Oct 20,2023
#revisited on Nov 29, 2023
############################################

rm(list=ls())
library(dplyr); library(ggplot2)
library(data.table); library(readxl)
library(scales); library(ggpubr)
library(gridExtra)

setwd("C:/Users/tehal/OneDrive - A STAR/Misc/")
#setwd("/Users/TEHAL/OneDrive - A STAR/Misc/")

mytmp  <- fread("Growth_trajectory_and_GenMetS_data.txt") %>% data.frame()
head(mytmp); dim(mytmp)
hist(scale(mytmp$GenMetS), xlab="zscore of GenMetS score",
     main="Histogram of GenMetS \n in GUSTO chilren (N=1073)")

#zBMI_birth, zbmi3,zbmi4,..zbmi16 representing z-score sex adjusted BMI for
#age= c("at birth", "3w", "3m", "6m", "9m", "12m", "15m", 
#"18m", "24m", "36m", "48m", "54m","60m", "66m", "72m")


#100 times boostrap imputation to calculated the 95% CI. 

myr2_95CI <- function(mydata, mymodel){
  set.seed(1234)
  library(boot)
  r2 <- function(mydata, idx){
    fit<- lm( mymodel, data=mydata[idx, ])
    summary(fit)$r.square
  }
  myfun <- boot(mydata,r2, R=100) #permutation, 100
  
  fit<- lm( mymodel, data=mydata)
  tmpr2 <- summary(fit)$r.square
  
  if(tmpr2 > 0.0001) { tmpr2_formated<- round(tmpr2,4)}
  else {tmpr2_formated <- scientific(tmpr2, digits=4) } 
  
  tmppv<- scientific(summary(fit)$coef[2,4], digits=4)
  
  r2_permute <- myfun$t
  tmpsd <- sd(myfun$t)
  tmplow <- tmpr2 - 1.96 * tmpsd
  tmphigh <- tmpr2 + 1.96 * tmpsd
  if(tmplow > 0.0001) { tmplow <- round(tmplow,4)}
  else {tmplow <- scientific(tmplow, digits=4) } 
  
  if(tmphigh > 0.0001) { tmphigh <- round(tmphigh,4)}
  else {tmphigh <- scientific(tmphigh, digits=4) } 
  
  tmpci<-paste0("(", tmplow, ",", tmphigh, ")")
  
  tmpresult <-c(tmpr2_formated, tmpci, tmplow, tmphigh, tmppv)
  
  results_list <- list(r2_perm = r2_permute, sum_result = tmpresult)
  return(results_list)
}


################################all
r2_result <- matrix(rep(NA, 2*15*8 ), ncol=8)
r2_perm_GenMetS <- matrix(rep(NA, 100*15), ncol = 15)
r2_perm_PGS <- matrix(rep(NA, 100*15), ncol = 15)
colnames(r2_result)<- c("target", "predictor", "N", "R2", "95%CI", "low", "high", "p-value")
K=1

for(i in c(91:105)){
  mybmi1 <- colnames(mytmp)[i]
  mymodel1 <- paste0(mybmi1, "~ GenMetS")
  
  tmp1 <- myr2_95CI(mytmp, mymodel1)
  N <- length(!is.na(mytmp[, i]))
  r2_result[K, ]<-c(mybmi1, "GenMetS",   N, tmp1$sum_result)
  r2_perm_GenMetS[,K] <- tmp1$r2_perm 
  K=K+1
}

for(i in c(91:105)){
  mybmi1 <- colnames(mytmp)[i]
  mymodel1 <- paste0(mybmi1, "~ BMI_PRS")
  
  tmp1 <- myr2_95CI(mytmp, mymodel1)
  N <- length(!is.na(mytmp[, i]))
  r2_result[K, ]<-c(mybmi1, "BMI_PGS",  N, tmp1$sum_result)
  r2_perm_PGS[,K-15] <- tmp1$r2_perm
  K=K+1
}


dplot2 <- r2_result[1:30, ] %>% data.frame()
tmp1 <- data.frame(
  target=unique(dplot2$target),
  age=c("at birth", "3w", "3m", "6m", "9m", "12m", "15m", 
        "18m", "24m", "36m", "48m", "54m","60m", "66m", "72m")
)
dplot2 <- merge(dplot2, tmp1, by="target")
dplot2$age <- factor(dplot2$age, levels = tmp1$age, labels = tmp1$age)                               
dplot2$R2 <- as.numeric(dplot2$R2)
dplot2$low <- as.numeric(dplot2$low)
dplot2$high <- as.numeric(dplot2$high)

#############
## perform 1-sided t-test between GenMetS and PGS
pval_plot <- matrix(rep(NA, 3*15), ncol = 3)
pval_plot[,1] <- t(c("at birth", "3w", "3m", "6m", "9m", "12m", "15m", 
                   "18m", "24m", "36m", "48m", "54m","60m", "66m", "72m"))
for (i in 1:15) {
    pval_test <- signif(t.test(r2_perm_GenMetS[,i], r2_perm_PGS[,i], 
                             alternative = "greater")$p.value, 4) 
    pval_sym <- ifelse(pval_test < 0.0005, "***", ifelse(pval_test < 0.005, "**", 
                                                  ifelse(pval_test < 0.05, "*", "")))
    pval_plot[i,2] <- pval_test
    pval_plot[i,3] <- pval_sym
      
}

colnames(pval_plot) <- c("age", "pval", "signif_level")
pval_plot

dplot3 <- merge(dplot2, pval_plot, by = "age")
dplot3$age <- factor(dplot3$age, levels = tmp1$age, labels = tmp1$age)
dplot3$signif_level[dplot3$predictor == "BMI_PGS"] <- ""
head(dplot3)
dim(dplot3)
 
myplot <- ggplot(data=dplot3, aes(x=age, y=R2, fill=predictor, 
                        label = signif_level)) +
  geom_bar(stat="identity", color="black", position=position_dodge(), width=0.8)+
  geom_errorbar(aes(ymin=R2, ymax=high), width=.4,position=position_dodge(.8)) +
  geom_text(aes(label=signif_level, y = 0.061), 
 vjust = -1.5, hjust = 1.1, position = position_dodge(width = 1), size = 7)+
  scale_fill_manual(values = c("grey90", "grey20")) +
  ylim(0, 0.07) +
  theme_bw() + 
  labs(y=expression ("Coefficient determination of Observed MetS "~(R^2)))  +
  theme(axis.text=element_text(colour="black",size=12),
        legend.text=element_text(colour="black", size=12),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank()
  ) +
  theme(panel.grid.minor = element_blank()) +
  theme(legend.position = "bottom") +
  #ylab("R2") +
  xlab("") +
  ggtitle("Explanation on Longitudinal Measures of Child BMI")

myplot

write.table(dplot3, "Information-to-plot-R2-for-comparison.txt", row.names = F,
            sep = "\t")
ggsave("Comparison_BMIPRS_GenMetS_to_longitudinalBMI_child.tiff", myplot,
       width =7.71, height =5.5, units = "in",dpi = 300, compression = "lzw")
cat("plot Comparison_BMIPRS_GenMetS_to_longitudinalBMI_child.tiff \n")

ggsave("Comparison_BMIPRS_GenMetS_to_longitudinalBMI_boys_and_girls.pdf", myplot,
       width =7.71, height = 5.5, units = "in",dpi = 300 )


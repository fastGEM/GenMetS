######################
#Pan Hong

#GenMetS vs disease age-of-onset by c-index
#Dec 12, 2023
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

library(survcomp) #for c-index
library(MatchIt)  #for matching
library(survival) #for survival analysis
sessionInfo()

#############################plot
library(dplyr)
library(forestploter)
library(grid)

#UKB women
tm <- forest_theme(base_size = 10,
                   # Confidence interval point shape, line type/color/width
                   ci_pch = 15,
                   ci_col = "#999999",
                   ci_fill = "#FFC000FF",
                   ci_alpha = 1,
                   ci_lty = 1,
                   ci_lwd = 1.5,
                   ci_Theight = 0.2, # Set an T end at the end of CI 
                   # Reference line width/type/color
                   refline_lwd = 1,
                   refline_lty = "dashed",
                   refline_col = "grey20",
                   # Vertical line width/type/color
                   vertline_lwd = 1,
                   vertline_lty = "dashed",
                   vertline_col = "grey20",
                   # Change summary color for filling and borders
                   summary_fill = "#4575b4",
                   summary_col = "#4575b4",
                   # Footnote font size/face/color
                   footnote_cex = 0.8,
                   footnote_fontface = "italic",
                   footnote_col = "black") 

dt <- fread("myresult_cindex_ukb_women.txt") %>% data.frame() 
dt$NonEvents <- dt$SampleSize-dt$Events
dt$`C-index (95% CI)` <- ifelse(is.na(dt$cindex), "",
                                sprintf("%.2f (%.2f to %.2f)",
                                        dt$cindex, dt$lower, dt$upper))
dt$plot <- paste(rep("   ", 6), collapse = " ")

dt <- dt  %>% arrange(desc(cindex))

colnames(dt)[which(colnames(dt)=="disease")]<- "Disease"
colnames(dt)[which(colnames(dt)=="EventTime")]<- "Event Age"
colnames(dt)[which(colnames(dt)=="pvalue")]<- "p-value"
dt1 <- dt %>% dplyr::select(Disease, NonEvents, Events, 'Event Age', plot, `C-index (95% CI)`, 'p-value')


colnames(dt1)[5]<- ""
p<- forest(dt1[,c(1,4, 3, 2, 5:7)],
           est = dt$cindex,
           lower = dt$lower, 
           upper = dt$upper,
           sizes = 0.5,
           ci_column = 5,
           ref_line = 0.5,
           arrow_lab = c("GenMetS lower", "GenMetS higher"),
           xlim = c(0.1, 0.9),
           ticks_at = c(0.1, 0.5, 0.9),
           footnote = "Event Age is the lowest 10% age of 
                        disease onset (years)",
           title="UKB ASN women",
           theme = tm)

print(p)
ggplot2::ggsave("GenMetS_diseaseAOO_ukb_women.tiff",
                plot=p, 
                width = 7, units = "in",
                dpi = 300, compression = "lzw")

pdfv1="GenMetS_diseaseAOO_ukb_women.pdf"
ggsave(pdfv1, p, width = 8, height =5.34,units = "in")

#############################men
dt <- fread("myresult_cindex_ukb_men.txt") %>% data.frame() 
dt$NonEvents <- dt$SampleSize-dt$Events
dt$`C-index (95% CI)` <- ifelse(is.na(dt$cindex), "",
                                sprintf("%.2f (%.2f to %.2f)",
                                        dt$cindex, dt$lower, dt$upper))
dt$plot <- paste(rep("   ", 6), collapse = " ")

dt <- dt  %>% arrange(desc(cindex))

tm <- forest_theme(base_size = 10,
                   # Confidence interval point shape, line type/color/width
                   ci_pch = 15,
                   ci_col = "#999999",
                   ci_fill = "#6F99ADFF",
                   ci_alpha = 1,
                   ci_lty = 1,
                   ci_lwd = 1.5,
                   ci_Theight = 0.2, # Set an T end at the end of CI 
                   # Reference line width/type/color
                   refline_lwd = 1,
                   refline_lty = "dashed",
                   refline_col = "grey20",
                   # Vertical line width/type/color
                   vertline_lwd = 1,
                   vertline_lty = "dashed",
                   vertline_col = "grey20",
                   # Change summary color for filling and borders
                   summary_fill = "#4575b4",
                   summary_col = "#4575b4",
                   # Footnote font size/face/color
                   footnote_cex = 0.8,
                   footnote_fontface = "italic",
                   footnote_col = "black") 


colnames(dt)[which(colnames(dt)=="disease")]<- "Disease"
colnames(dt)[which(colnames(dt)=="EventTime")]<- "Event Age"
colnames(dt)[which(colnames(dt)=="pvalue")]<- "p-value"
dt1 <- dt %>% dplyr::select(Disease, NonEvents, Events, 'Event Age', plot, `C-index (95% CI)`, 'p-value')


colnames(dt1)[5]<- ""
p<- forest(dt1[,c(1,4, 3, 2, 5:7)],
           est = dt$cindex,
           lower = dt$lower, 
           upper = dt$upper,
           sizes = 0.5,
           ci_column = 5,
           ref_line = 0.5,
           arrow_lab = c("GenMetS lower", "GenMetS higher"),
           xlim = c(0.1, 0.9),
           ticks_at = c(0.1, 0.5, 0.9),
           footnote = "Event Age is the lowest 10% age of 
                        disease onset (years)",
           title = "UKB ASN men", 
           theme = tm)

print(p)
ggplot2::ggsave("GenMetS_diseaseAOO_ukb_men.tiff",
                plot=p, 
                width = 7, units = "in",
                dpi = 600, compression = "lzw")

pdfv1="GenMetS_diseaseAOO_ukb_men.pdf"
ggsave(pdfv1, p, width = 8, height =5.34,units = "in")

############ATTRACT
#women
tm <- forest_theme(base_size = 10,
                   # Confidence interval point shape, line type/color/width
                   ci_pch = 15,
                   ci_col = "#999999",
                   ci_fill = "#FFC000FF",
                   ci_alpha = 1,
                   ci_lty = 1,
                   ci_lwd = 1.5,
                   ci_Theight = 0.2, # Set an T end at the end of CI 
                   # Reference line width/type/color
                   refline_lwd = 1,
                   refline_lty = "dashed",
                   refline_col = "grey20",
                   # Vertical line width/type/color
                   vertline_lwd = 1,
                   vertline_lty = "dashed",
                   vertline_col = "grey20",
                   # Change summary color for filling and borders
                   summary_fill = "#4575b4",
                   summary_col = "#4575b4",
                   # Footnote font size/face/color
                   footnote_cex = 0.8,
                   footnote_fontface = "italic",
                   footnote_col = "black") 




dt <- fread("myresult_cindex_ukb_women.txt") %>% data.frame() 
dt$NonEvents <- dt$SampleSize-dt$Events
dt$`C-index (95% CI)` <- ifelse(is.na(dt$cindex), "",
                                sprintf("%.2f (%.2f to %.2f)",
                                        dt$cindex, dt$lower, dt$upper))
dt$plot <- paste(rep("   ", 6), collapse = " ")

dt <- dt  %>% arrange(desc(cindex))

colnames(dt)[which(colnames(dt)=="disease")]<- "Disease"
colnames(dt)[which(colnames(dt)=="EventTime")]<- "Event Age"
colnames(dt)[which(colnames(dt)=="pvalue")]<- "p-value"
dt1 <- dt %>% dplyr::select(Disease, NonEvents, Events, 'Event Age', plot, `C-index (95% CI)`, 'p-value')


colnames(dt1)[5]<- ""
p<- forest(dt1[,c(1,4, 3, 2, 5:7)],
           est = dt$cindex,
           lower = dt$lower, 
           upper = dt$upper,
           sizes = 0.5,
           ci_column = 5,
           ref_line = 0.5,
           arrow_lab = c("GenMetS lower", "GenMetS higher"),
           xlim = c(0.1, 0.9),
           ticks_at = c(0.1, 0.5, 0.9),
           footnote = "Event Age is the lowest 10% age of 
                        disease onset (years)",
           title="UKB ASN women",
           theme = tm)

print(p)
ggplot2::ggsave("GenMetS_diseaseAOO_ukb_women.tiff",
                plot=p, 
                width = 7, units = "in",
                dpi = 300, compression = "lzw")

pdfv1="GenMetS_diseaseAOO_ukb_women.pdf"
ggsave(pdfv1, p, width = 8, height =5.34,units = "in")

#############################men
dt <- fread("myresult_cindex_ukb_men.txt") %>% data.frame() 
dt$NonEvents <- dt$SampleSize-dt$Events
dt$`C-index (95% CI)` <- ifelse(is.na(dt$cindex), "",
                                sprintf("%.2f (%.2f to %.2f)",
                                        dt$cindex, dt$lower, dt$upper))
dt$plot <- paste(rep("   ", 6), collapse = " ")

dt <- dt  %>% arrange(desc(cindex))

tm <- forest_theme(base_size = 10,
                   # Confidence interval point shape, line type/color/width
                   ci_pch = 15,
                   ci_col = "#999999",
                   ci_fill = "#6F99ADFF",
                   ci_alpha = 1,
                   ci_lty = 1,
                   ci_lwd = 1.5,
                   ci_Theight = 0.2, # Set an T end at the end of CI 
                   # Reference line width/type/color
                   refline_lwd = 1,
                   refline_lty = "dashed",
                   refline_col = "grey20",
                   # Vertical line width/type/color
                   vertline_lwd = 1,
                   vertline_lty = "dashed",
                   vertline_col = "grey20",
                   # Change summary color for filling and borders
                   summary_fill = "#4575b4",
                   summary_col = "#4575b4",
                   # Footnote font size/face/color
                   footnote_cex = 0.8,
                   footnote_fontface = "italic",
                   footnote_col = "black") 


colnames(dt)[which(colnames(dt)=="disease")]<- "Disease"
colnames(dt)[which(colnames(dt)=="EventTime")]<- "Event Age"
colnames(dt)[which(colnames(dt)=="pvalue")]<- "p-value"
dt1 <- dt %>% dplyr::select(Disease, NonEvents, Events, 'Event Age', plot, `C-index (95% CI)`, 'p-value')


colnames(dt1)[5]<- ""
p<- forest(dt1[,c(1,4, 3, 2, 5:7)],
           est = dt$cindex,
           lower = dt$lower, 
           upper = dt$upper,
           sizes = 0.5,
           ci_column = 5,
           ref_line = 0.5,
           arrow_lab = c("GenMetS lower", "GenMetS higher"),
           xlim = c(0.1, 0.9),
           ticks_at = c(0.1, 0.5, 0.9),
           footnote = "Event Age is the lowest 10% age of 
                        disease onset (years)",
           title = "UKB ASN men", 
           theme = tm)

print(p)
ggplot2::ggsave("GenMetS_diseaseAOO_ukb_men.tiff",
                plot=p, 
                width = 7, units = "in",
                dpi = 600, compression = "lzw")

pdfv1="GenMetS_diseaseAOO_ukb_men.pdf"
ggsave(pdfv1, p, width = 8, height =5.34,units = "in")



####################
#Pan Hong
#Oct 5, 2023
#international criteria
#for early disease onset
#t2dm <40
#hypertension <=55
#CAD <55
#HF < 55
#NALFD <20
#stroke <55
##################################
library(dplyr)
library(forestploter)
library(grid)

 
tm <- forest_theme(base_size = 10,
                   # Confidence interval point shape, line type/color/width
                   ci_pch = 15,
                   ci_col = "#762a83",
                   ci_fill = "black",
                   ci_alpha = 0.8,
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
                   footnote_cex = 0.7,
                   footnote_fontface = "italic",
                   footnote_col = "black") 




dt <- fread("result_cindex.txt") %>% data.frame() 
dt$NonEvents <- dt$SampleSize-dt$Events
dt$`C-index (95% CI)` <- ifelse(is.na(dt$cindex), "",
                           sprintf("%.2f (%.2f to %.2f)",
                                   dt$cindex, dt$lower, dt$upper))
dt$plot <- paste(rep("   ", 6), collapse = " ")

dt <- dt  %>% arrange(pvalue)

colnames(dt)[which(colnames(dt)=="disease")]<- "Disease"
colnames(dt)[which(colnames(dt)=="EventTime")]<- "Event Age"
colnames(dt)[which(colnames(dt)=="pvalue")]<- "p-value"
dt1 <- dt %>% dplyr::select(Disease, NonEvents, Events, 'Event Age', plot, `C-index (95% CI)`, 'p-value')


colnames(dt1)[5]<- ""
p<- forest(dt1[,c(1:4, 5:7)],
            est = dt$cindex,
            lower = dt$lower, 
            upper = dt$upper,
            sizes = 0.5,
            ci_column = 5,
            ref_line = 0.5,
            arrow_lab = c("GenMetS lower", "GenMetS higher"),
            xlim = c(0.2, 0.8),
            ticks_at = c(0.3, 0.5, 0.7),
            footnote = "Event Age is the lowest 10% age of 
                        disease onset (years)",
            theme = tm)

print(p)
ggplot2::ggsave("GenMetS_disease_onset_10percentile_with_table.tiff",
                plot=p, 
                width = 7, units = "in",
                dpi = 600, compression = "lzw")
 
pdfv1="GenMetS_disease_onset_10percentile_with_table.pdf"
ggsave(pdfv1, p, width = 8, height =5.34,units = "in")

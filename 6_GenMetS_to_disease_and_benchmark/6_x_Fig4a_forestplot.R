rm(list=ls())
library(dplyr)
library(data.table)
library(forestplot)
library(ggplotify)
library(patchwork)

base_data<- fread("result_META_disease_prediction.csv")
base_data$mean <- round(base_data$Combined_OR,2)
base_data$lower <- round(base_data$Lower_CI,2)
base_data$upper <- round(base_data$Upper_CI,2)
base_data$OR1 <- paste0(base_data$mean,"(", base_data$lower, ",",
                        base_data$upper, ")")
base_data$N <- base_data$Combined_N

disease_name <-data.frame(
  disease= c("T2DM", "CAD", "HTN","HF", "Stroke", "MI"),
  fullname=c("Type 2 diabetes", "Coronary artery disease", "Hypertension","Heart failure",
             "Stroke","Myocardial infarction")) 


 
#################################### 
ForestPlot_for_one_disease_v1 <- function(this_disease){
  summary_data <- base_data %>% 
    dplyr::filter(Disease == this_disease) %>%
    dplyr::filter(Cohort == "META") %>% 
    dplyr::select(mean, lower, upper, OR1, N, Weight, Cohort) 
  
  data <- base_data %>% 
    dplyr::filter(Disease == this_disease) %>%
    dplyr::filter(Cohort != "META") %>% 
    dplyr::select(mean, lower, upper, OR1, N, Weight, Cohort) 
  
    idx <- which(disease_name$disease==this_disease)
  this_plot <- data |>
    forestplot(mean = mean,
               upper = upper,
               lower = lower,
               labeltext = c(Cohort, OR1, N),
               #clip = c(0.5, 3.5),
               title = disease_name$fullname[idx],
               align = "ccccl",
               xlog = FALSE,
               xlab = "OR (95%CI)",
               xticks = seq(0.5, 2.5, 0.5),
               fn.ci_norm = fpDrawCircleCI,
               graphwidth = unit(2, "cm"),
               boxsize = 0.2,
               zero = 1,
               colgap = unit(0.5, "mm"),  # space between columns
               lineheight = unit(0.2, "cm"),  # space between rows
               vertices = TRUE,
               lwd.zero = 1.5,
               lwd.xaxis = 1.5,
               lwd.ci = 1.0,
               # Reduce text sizes to squeeze space
               txt_gp = fpTxtGp(label = gpar(cex = 0.7),  # Reduced label text size
                                ticks = gpar(cex = 0.7),  # Reduced tick mark size
                                xlab = gpar(cex = 0.7),   # Reduced x-axis label size
                                title = gpar(cex = 0.8))) |>
    fp_add_lines() |> 
    fp_set_style(box = "#F8766D",
                 line = "#F8766D",
                 summary = "black")   |> 
    fp_add_header(Cohort = c("Cohort"),
                  OR1 = c("OR(95%CI)"),
                  N = c("N")
                 ) |> 
    fp_append_row(mean = summary_data$mean,
                  lower = summary_data$lower,
                  upper = summary_data$upper,
                  N = summary_data$N,
                 
                  Cohort = "Summary",
                  OR1 = summary_data$OR1,
                  #boxsize = 0.3,  # Reduce the box size to make it smaller
                  is.summary = TRUE) #|> 
   # fp_set_zebra_style("#EFEFEF")  
  
  print(this_plot)
  return(this_plot)
}

ForestPlot_for_one_disease_v1("T2DM")

ForestPlot_for_one_disease_v2 <- function(this_disease){
  
  summary_data <- base_data %>% 
    dplyr::filter(Disease == this_disease) %>%
    dplyr::filter(Cohort == "META") %>% 
    dplyr::select(mean, lower, upper, OR1, N ) 
  
  data <- base_data %>% 
    dplyr::filter(Disease == this_disease) %>%
    dplyr::filter(Cohort != "META") %>% 
    dplyr::select(mean, lower, upper, OR1, N) 
  
  idx<-which(disease_name$disease == this_disease)
  this_plot <- data |>
    forestplot(mean = mean,
               upper = upper,
               lower = lower,
               labeltext = c(OR1, N ),
               #clip = c(0.5, 3.0),
               title = disease_name$fullname[idx],
               align = "cccl",
               xlog = FALSE,
               xlab = "OR (95%CI)",
               xticks = seq(0.5, 2.5, 0.5),
               fn.ci_norm = fpDrawCircleCI,
               graphwidth = unit(2, "cm"),
               boxsize = 0.2,
               zero = 1,
               # Reduce the gap between columns
               colgap = unit(0.5, "mm"),  # Decreased from 1mm to 0.5mm
               # Increase line height to reduce space between rows
               lineheight = unit(0.2, "cm"),  # Make rows closer together
               vertices = TRUE,
               lwd.zero = 1.5,
               lwd.xaxis = 1.5,
               lwd.ci = 1.0,
               # Reduce text sizes to squeeze space
               txt_gp = fpTxtGp(label = gpar(cex = 0.7),  # Reduced label text size
                                ticks = gpar(cex = 0.7),  # Reduced tick mark size
                                xlab = gpar(cex = 0.7),   # Reduced x-axis label size
                                title = gpar(cex = 0.8))) |>
    fp_add_lines() |> 
    fp_set_style(box = "#F8766D",
                 line = "#F8766D",
                 summary = "black")   |> 
    fp_add_header(OR1 = c("OR(95%CI)"),
                  N = c("N")) |> 
    fp_append_row(mean = summary_data$mean,
                  lower = summary_data$lower,
                  upper = summary_data$upper,
                  N = summary_data$N,
                  OR1 = summary_data$OR1,
                  #boxsize = 0.3,  # Reduce the box size to make it smaller
                  is.summary = TRUE) #|> 
    #fp_set_zebra_style("#EFEFEF")  
  
  print(this_plot)
  return(this_plot)
}
 
ForestPlot_for_one_disease_v2("HF")
 

T2DM_plot <- ForestPlot_for_one_disease_v1("T2DM")
HF_plot <- ForestPlot_for_one_disease_v2("HF")
HTN_plot <- ForestPlot_for_one_disease_v2("HTN")
MI_plot <- ForestPlot_for_one_disease_v1("MI")
Stroke_plot <- ForestPlot_for_one_disease_v2("Stroke")
CAD_plot <- ForestPlot_for_one_disease_v2("CAD")


p1 <- grid2grob(print(T2DM_plot))
p2 <- grid2grob(print(HF_plot))
p3 <- grid2grob(print(HTN_plot))
p4 <- grid2grob(print(MI_plot))
p5 <- grid2grob(print(Stroke_plot))
p6 <- grid2grob(print(CAD_plot))

library(gridExtra)
library(grid)


pdf("Fig4a.pdf", width = 8.2, height = 4.2 )
grid.arrange(p1, p2, p3,
             p4, p5, p6, 
             ncol = 3, 
             widths = c(1.34, 1.1, 1.1)) 
dev.off()

tiff("Fig4a.tiff", width = 8.2, height = 4.2, units = "in",
     res = 300, compression = "lzw")
grid.arrange(p1, p2, p3,
             p4, p5, p6, 
             ncol = 3, 
             widths = c(1.34, 1.1, 1.1)) 
dev.off()




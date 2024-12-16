rm(list=ls())
library(data.table)
library(dplyr)
library(forestplot) 
#
 
dat1 <- fread("data_table.txt")
#three significnat digits. 
format_three_digits <- function(x) {
  if (x >=10) {
    sprintf("%.1f", x)
  } else if (x >= 1) {
    # Format with three significant figures for values >= 1
    sprintf("%.2f", x)
  } else if (x >= 0.1) {
    # Format with three decimal places for values between 0.01 and 1
    sprintf("%.3f", x)
  } else if (x >= 0.01) {
    sprintf("%.4f", x)
  } else {
    # Scientific notation with three significant digits for values < 0.01
    sprintf("%.2e", x)
  }
}

formatted_r2 <- sapply(dat1$R2*100, format_three_digits) ; formatted_r2
formatted_pv <- sapply(dat1$p.value, format_three_digits) ; formatted_pv

myplot1 <- dat1 |> 
  mutate(R2_s = formatted_r2) |>
  mutate(pv_s = formatted_pv) |> 
  filter(sex == "all") |>
  mutate(cohort2 = c("UKB", 
                     strsplit(unique(dat1$cohort)[1], " ")[[1]][2],
                     rep("", 6),
                     "UKB",
                     strsplit(unique(dat1$cohort)[2], " ")[[1]][2],
                     rep("", 6),
                      unique(cohort)[3], rep("", 8),
                      unique(cohort)[4], "Children", rep("", 7))) |>
  forestplot( mean = R2,
              lower = low,
              upper = high,
              labeltext = c(cohort2, target, N, R2_s, pv_s),
              title = "All",
              align = "clccc",
              #clip = c(-0.025,0.2),
              #line.margin = .15,
              xlog = FALSE,
              xlab = "prediction R2 with 95% CI",
              xticks = seq(-0.025, 0.2, 0.05),
              fn.ci_norm = fpDrawCircleCI,
              graph.pos = 3,
              boxsize = 0.5,
              lineheight = unit(2, "cm"),
              colgap = unit(3, "mm"),
              #legend = c("Female","Male"),
              #legend_args = fpLegend(pos=list(x=0.9,y=0.03),
              #                       gp = gpar(fill = "lightgrey", 
              #                                 col = "lightgrey")),
              #graphwidth = unit(7, "cm"),
              vertices = TRUE,
              lwd.zero = 1.25,
              lwd.xaxis = 1.25,
              lwd.ci = 1.5,
              #title = "GenMetS Association with Continuous Traits",
              txt_gp = fpTxtGp(label = gpar(cex = 1.1), 
                               ticks = gpar(cex = 1), 
                               xlab = gpar(cex = 1.1),
                               title = gpar(cex = 1.3, align="center",
                                            fontface = "bold"))) |>
  fp_add_header(cohort2 = "Cohort\n",
                target = "Target\n",
                N = "N\n",
                R2_s = "R2(%)\n",
                #CI2 = "CI (95%)\n",
                pv_s = "p-value\n") |>
  fp_set_style(box = "#F8766D",
               #line = "#FFC000FF",
               hrz_lines = "#999999") |>
  #fp_decorate_graph(grid = structure(c(-0.027), gp = gpar(lty = 1, col = "lightgrey"))) |>
  #fp_set_zebra_style("#EFEFEF") |> 
  fp_add_lines(h_2 = gpar(lty = 1, col = "#070504", lwd = 1.2),
               h_10 = gpar(lty = 1),
               h_18 = gpar(lty = 1),
               h_27 = gpar(lty = 1))
               #h_36 = gpar(lwd = 1, columns = 4:6, col = "#070504"))

 myplot1

 

library(ggplotify)
library(patchwork)

p1 <- grid2grob(print(myplot1))


p_all <- wrap_elements(p1)  
  plot_layout(widths = c(5.5, -0.75, 4),guides = "collect") 
p_all


pdf("Fig_R2_GenMetS_to_traits_by_ALL.pdf", width = 8.2, height = 8.5)
par(mar=c(5,4,2,2))
p_all
dev.off()

tiff("Fig_R2_GenMetS_to_traits_by_ALL.tiff", width = 8.2, height = 8.5, units = "in",
       res = 300, compression = "lzw")
p_all
dev.off()




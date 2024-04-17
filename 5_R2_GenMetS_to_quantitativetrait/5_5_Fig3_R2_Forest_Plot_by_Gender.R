rm(list=ls())
library(data.table)
library(dplyr)
library(forestplot)

#setwd(mywd)
data_all <- fread("R2_result.txt")
dim(data_all); head(data_all)

colnames(data_all)[10]<- "p.value"
# check if any duplicate row
duplicated(data_all)
data_all[duplicated(data_all),]

# remove duplicate row
data <- data_all[duplicated(data_all) == "FALSE",]
dim(data); head(data)

# filter for only predictor == "GenMetS"
data <- data[data$predictor == "GenMetS",]
dim(data); head(data)

# change Abdominal circumference to AC*
# and waist circumference to WC*
data$target[data$target == "Abdominal Circumference"] <- "AC"
data$target[data$target == "Waist Circumference"] <- "WC"
 
# for UKBB European
dat1 <- data[cohort=="UKBB European"]
dat1 <- dat1[order(target, sex),]

# Find the row with "ObsMetS"
obsMetS_row <- dat1[target == "ObsMetS", ]

# Remove "ObsMetS" from dat1 to avoid duplication
dat1 <- dat1[target != "ObsMetS", ]

# sort by R2 of girls
sorted_girls <- dat1[sex == "female"][order(R2, decreasing = TRUE)]

# Extract the sorted target names
sorted_targets <- sorted_girls$target

# Reorder the original data frame based on the sorted target names
dat1 <- dat1[target %in% sorted_targets]
dat1 <- dat1[order(match(target, sorted_targets), sex),]

# Add "ObsMetS" row on top
dat1 <- rbind(obsMetS_row, dat1)


# for UKBB Asian
dat2 <- data[cohort=="UKBB Asian"]
dat2 <- dat2[order(target,sex),]

# Find the row with "ObsMetS"
obsMetS_row <- dat2[target == "ObsMetS", ]

# Remove "ObsMetS" from dat1 to avoid duplication
dat2 <- dat2[target != "ObsMetS", ]

# sort by R2 of girls
sorted_girls <- dat2[sex == "female"][order(R2, decreasing = TRUE)]

# Extract the sorted target names
sorted_targets <- sorted_girls$target

# Reorder the original data frame based on the sorted target names
dat2 <- dat2[target %in% sorted_targets]
dat2 <- dat2[order(match(target, sorted_targets), sex),]

# Add "ObsMetS" row on top
dat2 <- rbind(obsMetS_row, dat2)
dat1 <- rbind(dat1, dat2)

# ATTRaCT
dat3 <- data[cohort=="ATTRaCT"]
dat3 <- dat3[order(target,sex),]

# find and replace some words
dat3$target <- gsub("CHOLESTROL", "Cholestrol", dat3$target)
dat3$target <- gsub("FastingGlucose", "Fasting Glucose", dat3$target)

# Find the row with "ObsMetS"
obsMetS_row <- dat3[target == "ObsMetS", ]

# Remove "ObsMetS" from dat1 to avoid duplication
dat3 <- dat3[target != "ObsMetS", ]

# sort by R2 of girls
sorted_girls <- dat3[sex == "female"][order(R2, decreasing = TRUE)]

# Extract the sorted target names
sorted_targets <- sorted_girls$target

# Reorder the original data frame based on the sorted target names
dat3 <- dat3[target %in% sorted_targets]
dat3 <- dat3[order(match(target, sorted_targets), sex),]

# Add "ObsMetS" row on top
dat3 <- rbind(obsMetS_row, dat3)

dat1 <- rbind(dat1, dat3)
dim(dat1); head(dat1); tail(dat1)

# GUSTO
dat4 <- data[(cohort=="GUSTOchildren" | cohort == "GUSTOchildren6yr")]
dat4 <- dat4[order(target,sex),]

dat4$cohort <- "GUSTO"

# remove fatty liver index
dat4 <- dat4[-which(dat4$target == "Fatty Liver Index"),]

# Find the row with "ObsMetS"
obsMetS_row <- dat4[target == "ObsMetS", ]

# Remove "ObsMetS" from dat1 to avoid duplication
dat4 <- dat4[target != "ObsMetS", ]

# sort by R2 of girls
sorted_girls <- dat4[sex == "girl"][order(R2, decreasing = TRUE)]

# Extract the sorted target names
sorted_targets <- sorted_girls$target

# Reorder the original data frame based on the sorted target names
dat4 <- dat4[target %in% sorted_targets]
dat4 <- dat4[order(match(target, sorted_targets), sex),]

# Add "ObsMetS" row on top
dat4 <- rbind(obsMetS_row, dat4)

dat1 <- rbind(dat1, dat4)
dat1$sex[dat1$sex == "boy"] <- "male"
dat1$sex[dat1$sex == "girl"] <- "female"
dim(dat1); head(dat1); tail(dat1)
table(dat1$cohort)

#dat1$R2 <- format(signif(dat1$R2,4), scientific = T, digits = 3)
#dat1$low <- format(signif(dat1$low,4), scientific = T, digits = 3)
#dat1$high <- format(signif(dat1$high,4), scientific = T, digits = 3)
dat1$p.value <- format(signif(dat1$p.value,4), scientific = T, digits = 3)
#dat1$CI2 <- paste0("(", dat1$low, ",", dat1$high, ")")
  
# Plot
dat1$target[dat1$target == "ObsMetS"] <- "MetS"

myplot1 <- dat1 |> 
  mutate(R2_s = format(ifelse(signif(dat1$R2*100,4) < 0.01, 
                              0.00, signif(dat1$R2*100,4)), scientific = F, digit = 3),
         CI2 = paste0("(", format(signif(dat1$low,4), 
                                  scientific = T, digits = 3),
                      ",", format(signif(dat1$high,4), 
                                  scientific = T, digits = 3), ")")) |>
  filter(sex == "female") |>
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
              labeltext = c(cohort2, target, N, R2_s, p.value),
              title = "Female",
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
              colgap = unit(2, "mm"),
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
                p.value = "p-value\n") |>
  fp_set_style(box = "#FFC000FF",
               #line = "#FFC000FF",
               hrz_lines = "#999999") |>
  #fp_decorate_graph(grid = structure(c(-0.027), gp = gpar(lty = 1, col = "lightgrey"))) |>
  #fp_set_zebra_style("#EFEFEF") |> 
  fp_add_lines(h_2 = gpar(lty = 1, col = "#070504", lwd = 1.2),
               h_10 = gpar(lty = 1),
               h_18 = gpar(lty = 1),
               h_27 = gpar(lty = 1))
               #h_36 = gpar(lwd = 1, columns = 4:6, col = "#070504"))

#myplot1

#R2_s = format(signif(dat1$R2,4), scientific = T, digits = 3)

myplot2 <- dat1 |> 
  mutate(R2_s = format(ifelse(signif(dat1$R2*100,4) < 0.01, 0.00, 
                       signif(dat1$R2*100,4)), scientific = F, digit = 3),
         CI2 = paste0("(", format(signif(dat1$low,4), 
                                  scientific = T, digits = 3),
                      ",", format(signif(dat1$high,4), 
                                  scientific = T, digits = 3), ")")) |>
  filter(sex == "male") |>
  forestplot( mean = R2,
              lower = low,
              upper = high,
              labeltext = c(N, R2_s, p.value),
              title = "Male                            ",
              #clip = c(-0.025,0.2),
              #line.margin = .15,
              align = "ccc",
              lineheight = unit(2, "cm"),
              xticks = seq(-0.025, 0.125, 0.05),
              xlog = FALSE,
              xlab = "prediction R2 with 95% CI",
              #col = fpColors(box = "#6F99ADFF", lines = "#6F99ADFF"),
              fn.ci_norm = fpDrawCircleCI,
              graphwidth = unit(5, "cm"),
              graph.pos = 1,
              boxsize = 0.5,
              colgap = unit(2, "mm"),
              #legend = c("Female","Male"),
              #legend_args = fpLegend(pos=list(x=0.9,y=0.03),
              #                       gp = gpar(fill = "lightgrey", 
              #                                 col = "lightgrey")),
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
  fp_add_header(N = "N\n",
                R2_s = "R2(%)\n",
                #CI2 = "CI (95%)",
                p.value = "p-value\n") |>
  fp_set_style(box = "#6F99ADFF",
               #line = "#6F99ADFF",
               hrz_lines = "#999999") |>
  #fp_decorate_graph(grid = structure(c(-0.03), gp = gpar(lty = 2, col = "lightgrey"))) |>
  #fp_set_zebra_style("#EFEFEF") |> 
  fp_add_lines(h_2 = gpar(lty = 1, col = "#070504", lwd = 1.2),
               h_10 = gpar(lty = 1),
               h_18 = gpar(lty = 1),
               h_27 = gpar(lty = 1))
               #h_36 = gpar(lwd = 1, columns = 2:4, col = "#070504"))


#tg = gridExtra::tableGrob(dat1[dat1$sex == "female", c("cohort", "target")],
#                          theme = ttheme_minimal())  
#h = grid::convertHeight(sum(tg$heights), "in", TRUE)
#w = grid::convertWidth(sum(tg$widths), "in", TRUE)

library(ggplotify)
library(patchwork)

p1 <- grid2grob(print(myplot1))
p2 <- grid2grob(print(myplot2))


p_both <- wrap_elements(p1) + plot_spacer() + 
  wrap_elements(p2) + 
  plot_layout(widths = c(5.5, -0.75, 4),guides = "collect") 
p_both

#setwd("C:/Users/tehal/OneDrive - A STAR/Misc/")
pdf("Fig3_Forest_plot_by_gender.pdf", width = 12.5, height = 8.5)
par(mar=c(5,4,2,2))
p_both
dev.off()

tiff("Fig3_Forest_plot_by_gender.tiff", width = 12.5, height = 8.5, units = "in",
       res = 300, compression = "lzw")
p_both
dev.off()

write.table(dat1, "data_for_table.txt", row.names=F, sep="\t", quote=F)



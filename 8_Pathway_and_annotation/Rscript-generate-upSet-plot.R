rm(list=ls())
library(UpSetR)
library(data.table)

data <- fread("/Users/TEHAL/Downloads/GenMetS_SNP_source.txt", 
              header = T)
dim(data)
head(data)

table(data$trait)

#### keep only unique trait in each SNP
tmp <- NULL

for (i in 1:nrow(data)) {
  tmp[i] <- paste(unique(unlist(strsplit(data$trait[i], ";"))), 
                  collapse = ";")
  
}

head(tmp); length(tmp)
length(unique(tmp))
table(tmp)

####################################
#### unique traits
trait.uniq <- unique(unlist(strsplit(unique(tmp), ";")))
length(trait.uniq)
trait.uniq

##############################################
### create a matrix to indicate traits
trait_matrix <- matrix(0, ncol = 7, nrow = nrow(data))
dim(trait_matrix); head(trait_matrix)
colnames(trait_matrix) <- c("BP", "GLU", "HDL", "TG", "WC", "LDL", "Basal")

for (i in 1:length(tmp)) {
  trait_matrix[i,colnames(trait_matrix) %in% unlist(strsplit(tmp[i], ";"))] <- 1
}

dim(trait_matrix); head(trait_matrix)
tmp[1:20]
trait_matrix[1:20,]

data_merge <- data.frame(data, tmp, trait_matrix)
colnames(data_merge)[13] <- "trait_uniq"
dim(data_merge); head(data_merge)

######### upset plot
# Dataset
mydata <- data_merge[,c("snpid", "BP", "GLU", "HDL", "TG",
                        "WC", "LDL", "Basal")]
rownames(mydata) <- mydata$snpid
mydata <- mydata[,-1]
head(mydata); dim(mydata)

## make combination matrix
library(ComplexHeatmap)
m_comb <- make_comb_mat(mydata)
ss = set_size(m_comb)
cs = comb_size(m_comb)

# Plot
UpSet(m_comb, pt_size = unit(3, "mm"), lwd = 2,
      comb_col = c("red", "blue", "black", "orange",
                   "pink", "violet", "grey")[comb_degree(m_comb)],
      left_annotation = upset_left_annotation(m_comb),
      row_names_side = "left",
      comb_order = order(comb_size(m_comb), decreasing = T),
      set_order = order(set_size(m_comb), decreasing = TRUE))

setwd("/Users/TEHAL/OneDrive - A STAR/Misc/")
#pdf("GetMetS-SNPs-UpSet-plot.pdf", width = 11.5, height = 6.5)
tiff("GenMetS-SNPs-UpSet-Plot.tiff", width = 7.5, height = 5, units = "in", res = 300 )
UpSet(m_comb, pt_size = unit(2.5, "mm"), lwd = 2, 
      #left_annotation = upset_left_annotation(m_comb),
      #row_names_side = "left",
      top_annotation = HeatmapAnnotation(
        "Intersection size" = anno_barplot(cs, 
                                           ylim = c(0, max(cs)*1.1),
                                           border = FALSE, 
                                           gp = gpar(fill = "black"), 
                                           height = unit(4.5, "cm")),
        annotation_name_side = "left", 
        annotation_name_rot = 90),
      left_annotation = rowAnnotation(
        "Set size" = anno_barplot(-ss, baseline = 0,
                                  axis_param = list(at = c(0, -500, -1000, -1500, -2000),
                                                    labels = c(0, 500, 1000, 1500, 2000),
                                                    labels_rot = 0), border = FALSE, 
                                  gp = gpar(fill = "black"), 
                                  width = unit(3.5, "cm")),
        set_name = anno_text(set_name(m_comb), 
                             location = 0.5, 
                             just = "center",
                             width = max_text_width(set_name(m_comb)) + unit(2.5, "mm"))
      ),   
      #row_names_side = "left",
      right_annotation = NULL,
      show_row_names = F,
      comb_order = order(comb_size(m_comb), decreasing = T),
      set_order = order(set_size(m_comb), decreasing = TRUE))
dev.off()

# filter by combination size
UpSet(m_comb[comb_size(m_comb) >= 25], pt_size = unit(3, "mm"), lwd = 2, 
      left_annotation = upset_left_annotation(m_comb[comb_size(m_comb) >= 25]),
      row_names_side = "left",
      comb_order = order(comb_size(m_comb[comb_size(m_comb) >= 25]), decreasing = T))
      #set_order = order(set_size(m_comb), decreasing = TRUE))



##########################################################
###### filter for top 500 SNPs by weight
subdata <- head(data_merge[order(abs(data_merge$weight), decreasing = T),],500)
dim(subdata); head(subdata)

# Dataset
mysubdata <- subdata[,c("snpid", "BP", "GLU", "HDL", "TG",
                        "WC", "LDL", "Basal")]
rownames(mysubdata) <- mysubdata$snpid
mysubdata <- mysubdata[,-1]
head(mysubdata); dim(mysubdata)

## make combination matrix
m_comb2 <- make_comb_mat(mysubdata)

ss = set_size(m_comb2)
cs = comb_size(m_comb2)

#pdf("Top500-SNPs-by-absolute-weight.pdf", width = 11.5, height = 6)
tiff("Top500-SNPs-by-absolute-weight.tiff", width = 7, height = 5, units = "in", res = 300)
UpSet(m_comb2, pt_size = unit(3.5, "mm"), lwd = 2, 
      top_annotation = HeatmapAnnotation(
        "Intersection size" = anno_barplot(cs, 
                       ylim = c(0, max(cs)*1.1),
                       border = FALSE, 
                       gp = gpar(fill = "black"), 
                       height = unit(4.5, "cm")),
        annotation_name_side = "left", 
        annotation_name_rot = 90),
      left_annotation = rowAnnotation(
        "Set size" = anno_barplot(-ss, baseline = 0,
                     axis_param = list(at = c(0, -50, -100, -150),
                     labels = c(0, 50, 100, 150),
                     labels_rot = 0), border = FALSE, 
                     gp = gpar(fill = "black"), 
                     width = unit(3.5, "cm")),
        set_name = anno_text(set_name(m_comb2), 
                             location = 0.5, 
                             just = "center",
                             width = max_text_width(set_name(m_comb2)) + unit(2.5, "mm"))
      ),   
      #row_names_side = "left",
      right_annotation = NULL,
      show_row_names = F,
      comb_order = order(comb_size(m_comb2), decreasing = T),
      set_order = order(set_size(m_comb2), decreasing = TRUE))  
dev.off()



        

  
  
  
  
  
  
  
  
  
  
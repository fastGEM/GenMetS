rm(list=ls())
library(data.table)
library(ggplot2)
library(dplyr); library(grid)

mywd <- "/Users/TEHAL/OneDrive - A STAR/Misc/MetS-PH/FUMA_gene2func128645/"

setwd(mywd)
data <- fread("GS.txt")
data$Prop <- data$N_overlap/data$N_genes

#filter for KEGG gene set
Kegg <- data[data$Category == "KEGG",]
Kegg <- Kegg[order(Kegg$adjP),]
dim(Kegg); head(Kegg[,c(1:6,9)])

# sort Kegg geneset by p-value
Kegg_plot <- Kegg[,] |> arrange(desc(adjP))
Kegg_plot$GeneSet2 <- factor(gsub("KEGG_", "", Kegg_plot$GeneSet), 
                             levels = gsub("KEGG_", "", Kegg_plot$GeneSet))
dim(Kegg_plot)

p1 <- Kegg_plot |> 
  ggplot(aes(x=GeneSet2, y= Prop)) + 
  geom_col(fill='red') + 
  theme_bw()+
  coord_flip() + 
  scale_y_reverse(name= "Proportion overlap", expand = expand_scale(mult= c(c(0.05,0)))) +
  scale_x_discrete() +
  theme(panel.spacing.x = unit(0.05, "mm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size=8, color = "black"),
        axis.line = element_line(colour = "black", 
        size = 0.7, linetype = "solid"),
        panel.border = element_blank()
        ) +
  theme(plot.margin = unit(c(5.5, 0, 5.5, 5.5), "pt"))

p2 <- Kegg_plot |> 
  ggplot(aes(x=GeneSet2, y= -log10(adjP))) + 
  geom_col(fill='blue') + 
  scale_y_continuous(name = "-log10 adjP", expand = expand_scale(mult= c(c(0,0.05)))) +
  scale_x_discrete() +
  coord_flip() +
  theme(panel.spacing.x = unit(0, "mm"))+ 
  theme_minimal() +
  theme(axis.title.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(colour = "black", 
        size = 0.7, linetype = "solid"),
        plot.margin = unit(c(5.5, 5.5, 5.5, -5.6), "pt"))

grid.newpage()
#grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
p <- arrangeGrob(p1, p2, nrow = 1, widths = c(1.5, 1))
grid.draw(p)
#Kegg.plot <-  gridExtra::grid.arrange(p1, p2, nrow = 1, widths = c(1.5, 1))

ggsave("FUMA-Pathway-analysis-results-Kegg.pdf", p, dpi = 300) 
#width =10.5, height = 5.5, units = "in",
#dpi = 300)

tiff("FUMA-Pathway-analysis-results-Kegg.tiff", width = 12, height = 7, units = "in",
     res = 300, compression = "lzw")
grid.draw(cbind(ggplotGrob(p1), ggplotGrob(p2), size = "last"))
dev.off()

# filter for reactome geneset
reactome <- data[data$Category == "Reactome",]
reactome <- reactome[order(reactome$adjP),]
dim(reactome); head(reactome[,c(1:6,9)])

# sort canonical geneset by p-value
reactome_plot <- reactome[reactome$adjP < 0.00001,] |> arrange(desc(adjP))
reactome_plot$GeneSet2 <- factor(gsub("REACTOME_", "", reactome_plot$GeneSet), 
                                levels = gsub("REACTOME_", "", reactome_plot$GeneSet))
dim(reactome_plot)


p1c <- reactome_plot |> 
  ggplot(aes(x=GeneSet2, y= Prop)) + 
  geom_col(fill='red') + 
  theme_bw()+
  coord_flip() + 
  scale_y_reverse(name= "Proportion overlap", expand = expand_scale(mult= c(c(0.05,0)))) +
  scale_x_discrete() +
  theme(panel.spacing.x = unit(0, "mm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size=8, color = "black"),
        axis.line = element_line(colour = "black", 
                                 size = 0.7, linetype = "solid"),
        panel.border = element_blank()
  ) +
  theme(plot.margin = unit(c(5.5, 0, 5.6, 5.5), "pt"))

p2c <- reactome_plot |> 
  ggplot(aes(x=GeneSet2, y= -log10(adjP))) + 
  geom_col(fill='blue') + 
  scale_y_continuous(name = "-log10 adjP", expand = expand_scale(mult= c(c(0,0.05)))) +
  scale_x_discrete() +
  coord_flip() +
  theme(panel.spacing.x = unit(0, "mm"))+ 
  theme_minimal() +
  theme(axis.title.y=element_blank(), 
        axis.text.y=element_blank(),
        axis.line.y = element_blank(), 
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line.x = element_line(colour = "black", 
                                   size = 0.7, linetype = "solid"),
        #plot.margin = unit(c(5.2, 5.3, 5.5, -5.75), "pt"))
        plot.margin = unit(c(5.5,5.85,5.5,-2), "pt"))

grid.newpage()
#grid.draw(cbind(ggplotGrob(p1c), ggplotGrob(p2c), size = "last"))
reactome.plot <-  gridExtra::grid.arrange(p1c, p2c, nrow = 1, widths = c(1.85, 1))
c <- arrangeGrob(p1c, p2c, nrow = 1, widths = c(1.9, 1))
grid.draw(c)

#pdf("FUMA-Pathway-analysis-results.pdf", width = 12, height = 7)
#gridExtra::grid.arrange(p1, p2, nrow = 1, widths = c(1.5, 1))
#dev.off()

ggsave("FUMA-Pathway-analysis-results-Reactome.pdf", c, dpi = 300) 

tiff("FUMA-Pathway-analysis-results-Reactome.tiff", width = 12, height = 7, units = "in",
     res = 300, compression = "lzw")
grid.draw(cbind(ggplotGrob(p1c), ggplotGrob(p2c), size = "last"))
dev.off()

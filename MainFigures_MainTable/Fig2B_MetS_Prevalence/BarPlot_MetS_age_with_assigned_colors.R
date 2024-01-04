##################################################
#Pan Hong
#R circular stacked barplot
#updated on Dec 5, 2023

#updated on Dec8, 2023
#flat barplot.
###################################################
rm(list=ls())
library(dplyr)
library(data.table)
library(ggplot2)
 
mydata <- fread("cohort_ethnicity_MetS_statistics.txt")
#mydata <- mydata %>% dplyr::filter(ethnicity != "OtherAsian")
table(mydata$eth_cohort, mydata$MetS_score_category)



#need to make up the missing rows
#all Chinese==> Northeast Asian
#all Malay  ==> Southeast Asian
#all Indian ==> Southwest Asian
#all European ==> Northern European
#all OtherAsian in sgp ==> Southeast Asian
#all OtherAsian in UKB ==> Southwest Asian
mydata<- mydata %>% add_row(eth_cohort="Southeast Asian:ATTRaCT:Male", 
                   MetS_score_category=0, 
                   n=0, 
                   freq=0.0, 
                   cohort="ATTRaCT", 
                   ethnicity="Southeast Asian", 
                   sex="Male") %>%
            add_row(eth_cohort="Southeast Asian:UKB:Male", 
                    MetS_score_category=0, 
                    n=0, 
                    freq=0.0, 
                    cohort="UKB", 
                    ethnicity="Southeast Asian", 
                    sex="Female") %>%
            add_row(eth_cohort="Southwest Asian:ATTRaCT:Male", 
                    MetS_score_category=0, 
                    n=0, 
                    freq=0.0, 
                    cohort="ATTRaCT", 
                    ethnicity="Southwest Asian", 
                    sex="Male")  
table(mydata$eth_cohort, mydata$MetS_score_category)

mydata$eth_cohort <- gsub("Northeast Asian", "NEA", mydata$eth_cohort)
mydata$eth_cohort <- gsub("Southeast Asian", "SEA", mydata$eth_cohort)
mydata$eth_cohort <- gsub("Southwest Asian", "SWA", mydata$eth_cohort)
mydata$eth_cohort <- gsub("Northern European", "EUR", mydata$eth_cohort)
 

data <- mydata
data$individual<- mydata$eth_cohort
data$group <- mydata$cohort
data$observation <- paste0("s", mydata$MetS_score_category)
data$value <- mydata$freq
data$sex <- mydata$sex
data <- data %>% dplyr::select(individual, group, observation, value, sex)

data$group <- factor(data$group, levels=c("SGP", "UKB", "ATTRaCT"),
                     labels=c("Discovery", "UKB", "ATTRaCT"))

# Set a number of 'empty bar' to add at the end of each group
empty_bar <- 1
nObsType <- nlevels(as.factor(data$observation))
tmp <- nlevels(as.factor(data$group))
to_add <- data.frame( matrix(NA, empty_bar*tmp*nObsType, ncol(data)) )
colnames(to_add) <- colnames(data)
to_add$group <- rep(levels(as.factor(data$group)), each=empty_bar*nObsType )
data <- rbind(data, to_add)
data <- data %>% arrange(group, individual)
data$id <- rep( seq(1, nrow(data)/nObsType) , each=4)

# Get the name and the y position of each label
label_data <- data %>% group_by(id, individual) %>% summarize(tot=sum(value))
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     
label_data$hjust <- ifelse( angle < -120, 1.1, -0.1)
label_data$angle <- ifelse( angle < -120, angle+180, angle)

label_data$individual<- gsub(":ATTRaCT|:UKB|:SGP", "", label_data$individual)
label_data$individual<- gsub(":Female|:Male", "", label_data$individual)
x<- label_data %>% group_by(individual)
label_data[10:16, ] 

data$observation<- paste0(data$sex, ":", as.character(data$observation))
data$observation <- factor(data$observation ,
                           levels=c(paste0("Female:", c("s3", "s2", "s1", "s0")), 
                                    paste0("Male:", c("s3", "s2", "s1", "s0"))),
                           labels=c(paste0("Female:", c("3", "2", "1", "0")), 
                                    paste0("Male:", c("3", "2", "1", "0"))))


# prepare a data frame for base lines
base_data <- data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

#myversion
#mycolorcode <- rev(paletteer_c("ggthemes::Red", 4))
#mycolorcode1 <- c("#AE123AFF", "#F8826BFF",  "#FFBEB2FF",  "#FFFFFFFF")
#mycolorcode2 <- rev(paletteer_c("ggthemes::Blue", 4))
#mycolorcode2 <- c("#2A5783FF", "#7EAED3FF", "#B9DDF1FF", "#FFFFFFFF")
#mycolorcode <-c(mycolorcode1, mycolorcode2)

#Evelyn 1
mycolorcode1<- c("#DC0000FF", "#DC000080", "#DC000025", "#DC000000")
mycolorcode2<- c('#91D1C2FF', '#91D1C290', '#91D1C240', '#91D1C200')
mycolorcode <-c(mycolorcode1, mycolorcode2)

#Evelyn 2
mycolorcode1<-c("#EFC000FF", "#EFC00090", "#EFC00040", "#EFC00000")
mycolorcode2<-c("#6F99ADFF", "#6F99AD95", "#6F99AD40", "#6F99AD00")
mycolorcode <-c(mycolorcode1, mycolorcode2)

library("scales")
library("ggsci")
show_col(mycolorcode)
 

# Make the plot
legend_title <- "MetS"
ggplot(data) +      
  # Add the stacked bar
  geom_bar(aes(x=as.factor(id), y=value, fill=(observation)), colour="black" , 
           stat="identity",width=0.8) +
  scale_fill_manual(legend_title, values=mycolorcode)  +
  # Add text showing the value of each 100/75/50/25 lines
  ggplot2::annotate("text", x = rep(20,5), 
                    y = c(0, 0.25, 0.50, 0.75, 1.0), 
                    label = c("0", "0.25", "0.5", "0.75", "1.0") , 
                    color="black", size=3 , angle=0, fontface="bold", hjust=0.4) +
  ylim(-0.3,max(label_data$tot+0.4, na.rm=T)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.key=element_rect(colour="black"),
    legend.key.size = unit(0.5, 'cm'),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(0,5), "cm") 
  ) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE)) +
  # Add labels on top of each bar
  geom_text(data=label_data, aes(x=id, y=tot, label=individual, hjust=-0.1), 
            color="black", fontface="bold",alpha=0.9, size=4, 
            angle= 90, inherit.aes = FALSE ) +
  #bottom labels. 
  geom_text(data=base_data, aes(x = title, y =-0.1, label=group), 
            hjust=c(0.5,0.5,0.2), colour = "black", alpha=0.9, size=4, fontface="bold", 
            inherit.aes = FALSE)


pdfv1="Fig2B.pdf"
imgv1="Fig2B.tiff"
ggsave(imgv1, width = 6.71,height =4,units = "in", dpi = 300, compression = "lzw")
ggsave(pdfv1, width = 6.71, height =4,units = "in")

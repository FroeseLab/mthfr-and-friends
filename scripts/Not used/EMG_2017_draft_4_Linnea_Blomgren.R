##### INTRO #####
# File description ----

#Make some type of filtering think to filter out each individual Substrateple OK!
#Filter data based on well.Position  OK!
#Take out one Temperature furthest in column 1 and fluorescence filling up OK!
#Make a new table with fluorescence vs its own temperature -- later  OK!
# Define the inclination for each data OK!
#Combine dublicates 
#Show one raw flot and one fixed
#Remove background checks on fixes 

# Libraries ----
library(dplyr)
library(tidyr)

library(inflection)
library(data.table)
library(DataCombine)
library(scales)
library(sicegar) #Curve fitting 
require("ggrepel")
library(grid)



library(ggplot2)
library(cowplot)

library(trackViewer)

# Data files ----
getwd()
# Set directory 
setwd(getwd())
#Set the folder I will be working in and save my data in directly , will also se checked box in files. 
#Can also just go there and save file path. 


data_sample<-data.frame(Site=c(  9, 10, 18, 20, 21, 23, 25, 26, 29, 30, 103, 394),
           aa="Ser",
           color="#95aaab")

data_sample<-data_sample %>%
rbind(data.frame(Site=c(34,94,451),
      aa="Thr",
      color="#b795be"))

data_sample<-data_sample %>%
  rbind(data.frame(Site=c(90),
                   aa="Tyr",
                   color="#d1af97"))
#set size of circles with exception for T34A 
data_sample$cex<-0.5
data_sample$cex[13]<-1
#set size of circles with exception for T34A 
data_sample$score<-c(0.9,1.5)
data_sample$score[13]<-1.9

#set size of circles with exception for T34A 
data_sample$border<-"gray"
data_sample$border[13]<-"black"





sample.gr <- GRanges(
  seqnames = "chr1",
  IRanges(
    start=data_sample$Site,
    width = 1
    #names = data_sample$Site #names above the lollipops 
  ),
  color = data_sample$color,
  fill = data_sample$color,
  shape="circle",
  cex=data_sample$cex, # circle size 
  lwd= 1, #line width of it all 
  score=data_sample$score,#make T34A pop out higher above 
  height=1

)

sample.gr$border <- data_sample$border

features <- GRanges("chr1", ranges =
                      IRanges(
                        start = c(1, 48, 338, 413),
                        end =   c(47, 337, 412, 656),
                        names = paste0(c("Ser-rich region", "CD", "Linker", "RD"))
                      ),
                    fill=c("pink", "#51C6E6", "red", "#DFA32D"),
                   height=0.03
                
                    )


legend <- list(
  ## legend fill color
  labels = unique(data_sample$aa),
  col = "black",
  fill = unique(data_sample$color, )
)

# Open PNG device
png("data/output/lolliplot.png", width=10, height=7, units = "in",  res=600)
lolliplot(sample.gr, features,xaxis=c(1,48,338,413,656),  legend=legend, yaxis=FALSE, ylab=" ", cex=0.5) #cex affects how close circles are 
dev.off()




?GRanges
?lolliplot
?IRanges

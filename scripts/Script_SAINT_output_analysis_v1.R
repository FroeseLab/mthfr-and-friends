library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(DataCombine)
library(UniProt.ws)
library(gridExtra)
library("ggrepel")

#Read the input file 
# Data files ----
getwd()

#setwd(getwd()) #sets file location as directory 
setwd("~/Programs/Git/mthfr-and-friends")

#Import raw csv file. Manually converted from output Scaffold .txt file prior. 
data_raw_SAINT<- read.csv("data/raw/21241_baitPrey_data_FCA0.0, SP0.9, WD-1.0, IS-1.0.tab.csv",  header = TRUE, check.names = FALSE) #it is just the E1 SYPRO 5x from the original raw data from 220420


#Retrieve gene names using mapUniProt
#data_accession_names<-mapUniProt( from="UniProtKB_AC-ID", to='RefSeq_Protein', query=unique(data_raw_SAINT$Prey) )
data_accession_names<-mapUniProt( from="UniProtKB_AC-ID", to='Gene_Name', query=unique(data_raw_SAINT$Prey) )
#allToKeys() #check which types of accession names we can retrieve 
  
#Add column gene_name for converted Uniprot accession numbers to Gene names 
data_raw_SAINT<-data_raw_SAINT %>%
  mutate(gene_name=
    FindReplace(data_raw_SAINT, "Prey", data_accession_names, from = "From", to = "To", exact = TRUE, vector = TRUE))



#Plotting----
##Plot SAINT SP vs FC-A all data ----
plot_SP_FCA <- ggplot(data = data_raw_SAINT, aes(x = SP, y = FCA, label = gene_name)) +
  geom_point(aes(color = Bait), show.legend = FALSE) +
  
  # #hihglight sample names
  geom_text_repel(
    data = filter(
      data_raw_SAINT,
      gene_name == "MTHFR" | gene_name == "MTHFD1" | gene_name == "GCN1" | gene_name == "PRKDC" | gene_name == "FASN"
    ),
    aes(x = SP, y = FCA, label = gene_name),
    size = 2.5,
    nudge_y = 0.5,
    segment.size = 0.2,
    arrow = arrow(length = unit(0.01, "npc")),
    box.padding = 0.2,
    point.padding = 0.3,
    min.segment.length = 0,
    max.overlaps = 30,
    show.legend = FALSE
  ) +
  
   theme_minimal() +
  ggtitle("SAINT analysis")+
  theme(plot.title = element_text(hjust = 0.5))

plot_SP_FCA

ggsave(
  "plot_SP_FCA.png",
  plot = plot_SP_FCA,
  height = 2000,
  width = 2000,
  units = "px"
)


##Plot SAINT SP vs FC-A with SP >0.95 and FCA>2 ----
plot_SP_FCA <- ggplot(data = filter(data_raw_SAINT, SP > 0.95, FCA > 2), aes(x = SP, y = FCA)) +
  geom_point(aes(color = Bait), show.legend = FALSE) +
  
  # #hihglight sample names
  geom_text_repel(
    aes(label = gene_name),
    size = 2.5,
    nudge_y = 1,
    segment.size = 0.2,
    arrow = arrow(length = unit(0.01, "npc")),
    box.padding = 0.2,
    point.padding = 0.3,
    min.segment.length = 0,
    max.overlaps = 20,
    show.legend = FALSE
  ) +
  
  theme_minimal() +
  ggtitle("SAINT analysis")+
  theme(plot.title = element_text(hjust = 0.5))
plot_SP_FCA

ggsave(
  "plot_SP_0.95_FCA_2.png",
  plot = plot_SP_FCA,
  height = 2000,
  width = 2000,
  units = "px"
)



##Make heatmap ----
#Set up plot settings
plot_heatmap <- ggplot(data = filter(data_raw_SAINT, SP >= 0.95, FCA >= 1.99) ,
                       aes(x = Bait, y = gene_name)) +
  #Gradient scale
  scale_fill_gradient(
    low = "blue",
    high = "red",
    guide = guide_colorbar(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      label.position = "bottom",
      ticks = TRUE
    ),
  ) +
  #Theme
  scale_x_discrete(expand = c(0, 0)) +  # Remove space on the x-axis
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0),
    # Tilt x-axis title
    legend.position = "top",
    legend.title = element_text(angle = 0, vjust = 1),
    legend.text = element_text(angle = 0),
    panel.background = element_rect(fill = "gray"),
    # Gray background
    panel.grid.major = element_blank(),
    # No major grid lines
    panel.grid.minor = element_blank(),
    # No minor grid lines
    panel.border = element_rect(color = "black", fill = NA),
    # Box around the plot
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

#Generate different plots visualising different SAINT results 
plot_heatmap_FCA<-plot_heatmap +
  geom_tile(aes(fill=FCA))
plot_heatmap_FCA

plot_heatmap_SP<-plot_heatmap +
  geom_tile(aes(fill=SP))   
plot_heatmap_SP
   
plot_heatmap_Abundance<-plot_heatmap +
  geom_tile(aes(fill=Abundance))   
plot_heatmap_Abundance

#visualise all plots together
plot_heatmaps<-grid.arrange(plot_heatmap_FCA, plot_heatmap_SP, plot_heatmap_Abundance, nrow=1)
ggsave("plot_heatmaps.png",plot=plot_heatmaps,  height = 2000, width=3000, units = "px")

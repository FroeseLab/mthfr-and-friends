library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(DataCombine)
library(UniProt.ws)
library(gridExtra)
library("ggrepel")

library(STRINGdb)



#Read the input file 
# Data files ----
getwd()

#setwd(getwd()) #sets file location as directory 
setwd("~/Programs/Git/mthfr-and-friends")

#SAINT raw data---
#Import raw .tsv file of the "raw SAINT results". When downloaded a .output file isgiven. Just add .tsv file format 
data_raw_SAINT<- read.table("data/raw/21262.output.tsv",  header = TRUE, check.names = FALSE) 
#Import the unfiltered matrix from "Bait-prey heatmap" tab. Downloaded as a .tab file. Just add .tsv file format
data_baitPrey<-read.table("data/raw/21262_baitPrey_data.tab.tsv",  header = TRUE, check.names = FALSE) 

#Merge SAINT files 
#Note preys that are found in raw data but not found in bait-prey data will not be merged. E.g. A0A075B6S2
data_SAINT_merged<- merge(data_raw_SAINT,data_baitPrey)




# #Retrieve gene names using mapUniProt
# data_accession_names<-mapUniProt( from="UniProtKB_AC-ID", to='Gene_Name', query=unique(data_raw_SAINT$Prey) )
# #allToKeys() #check which types of accession names we can retrieve 
#   
# #Add column gene_name for converted Uniprot accession numbers to Gene names 
# data_SAINT_merged<-data_SAINT_merged %>%
#   mutate(gene_name=
#     FindReplace(data_SAINT_merged, "Prey", data_accession_names, from = "From", to = "To", exact = TRUE, vector = TRUE))


#Plotting----
##SP vs FC-A all data ----
plot_SP_FCA <- ggplot(data = data_SAINT_merged, aes(x = SP, y = FCA, label = PreyGene)) +
  geom_point(aes(color = Bait), show.legend = FALSE) +
  
  # #hihglight sample names
  geom_text_repel(
    data = filter(
      data_SAINT_merged,
      PreyGene == "MTHFR" | PreyGene == "MTHFD1" | PreyGene == "GCN1" | PreyGene == "PRKDC" | PreyGene == "FASN"
    ),
    aes(x = SP, y = FCA, label = PreyGene),
    size = 2.5,
    nudge_y = 0.1,
    nudge_x = -0.1,
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
  "plot_SP_FCA_check.png",
  plot = plot_SP_FCA,
  height = 2000,
  width = 2000,
  units = "px",
  path ="data/output"
)

##SP vs FC-A facet ----
plot_SP_FCA <- ggplot(data = data_SAINT_merged, aes(x = SP, y = FCA, label = PreyGene)) +
  geom_hline(yintercept = 2, color="pink")+
  annotate("text", x = 0.1, y = 2.2 , label = "FCA = 2", color="pink")+
   geom_point(color="gray", show.legend = FALSE) +
  
  # #hihglight sample names
  geom_text_repel(
    data = filter(
      data_SAINT_merged,
      PreyGene == "MTHFR" | PreyGene == "MTHFD1" | PreyGene == "GCN1" | PreyGene == "PRKDC" | PreyGene == "FASN"
    ),
    aes(x = SP, y = FCA, label = PreyGene),
    size = 3.5,
    nudge_y = 0.2,
    nudge_x = -0.1,
    segment.size = 0.2,
    arrow = arrow(length = unit(0.01, "npc")),
    box.padding = 1,
    point.padding = 0.3,
    min.segment.length = 0,
    max.overlaps = 30,
    show.legend = FALSE
  ) +
 
  #Theme 
  facet_grid(~Bait)+
  theme_minimal() +
  ggtitle("SAINT analysis")+
  scale_y_continuous(trans='log10')+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))

plot_SP_FCA

ggsave(
  "plot_SP_FCA_facet.png",
  plot = plot_SP_FCA,
  height = 2000,
  width = 4000,
  units = "px",
  path ="data/output"
  
)





## SP vs FC-A filtered ----
plot_SP_FCA <- ggplot(data = filter(data_SAINT_merged, SP > 0.95, FCA > 2), aes(x = SP, y = FCA)) +
  geom_point(aes(color = Bait), show.legend = FALSE) +
  
  # #hihglight sample names
  geom_text_repel(
    aes(label = PreyGene),
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
  units = "px",
  path ="data/output"
)

##Heatmap ----
#Set up plot settings
plot_heatmap <- ggplot(data = filter(data_SAINT_merged, FDR<=0.1, FCA>=2) ,
                       aes(x = Bait, y = PreyGene)) +
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
    panel.background = element_rect(fill = "white"),
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
ggsave("plot_heatmaps.png",plot=plot_heatmaps,  height = 2000, width=3000, units = "px", path ="data/output" )


#Modify data ----
data_SAINT_merged_modified<- data_SAINT_merged

##Kinases -------------------------------------------------------------------------------
data_kinases<-read.csv("data/raw/uniprot_kinases.csv",  header = TRUE, check.names = FALSE)
data_SAINT_merged_modified<-data_SAINT_merged_modified %>%
  mutate(kinase =data_SAINT_merged$PreyGene %in% data_kinases$gene_name )
#intersect(data_SAINT_merged$PreyGene,data_kinases$gene_name)


## STRING--------------------------------------------------------------------------------
#Setup, usally use thresholwd above 400  
string_db <- STRINGdb$new( version="12.0", species=9606,
                           score_threshold=100, network_type="full", input_directory="")

#Gene names to look up in STRING
#dat_example = data.frame(list(genes=c('MTHFR', 'GCN1')))
dat_example = data.frame(list(genes= unique(data_SAINT_merged_modified$PreyGene)))

#Retrieve STRING accession ID for each gene name
dat_examples_mapped = string_db$map(dat_example, 'genes')

#Checks score for interactions between defined gene names, and re-formatting the file
#Interactions contain combined gene neighborhood, fusion and co-occurrence
data_string<- string_db$get_interactions(dat_examples_mapped$STRING_id) %>% 
  merge(dat_examples_mapped %>%
          rename(from_gene=genes),
        by.x='from', by.y='STRING_id') %>%
  merge(dat_examples_mapped %>%
          rename(to_gene=genes),
        by.x='to', by.y='STRING_id')

#Get only interactions with MTHFR
data_STRING_MTHFR<-data_string %>%
    unique() %>%
    filter(from_gene=="MTHFR" | to_gene=="MTHFR") %>%
    dplyr::select(c(combined_score, from_gene, to_gene)) %>%
    # Reshape it to have the interacting gene
    tidyr::gather('direction', 'gene',-combined_score ) %>%
    filter(gene != 'MTHFR')

#Merge to input combined scores and direction of interaction 
data_SAINT_merged_modified<- merge(data_SAINT_merged_modified,data_STRING_MTHFR, by.x ="PreyGene", by.y="gene", all=TRUE)



##Crapome----------------------------------------------------------------------------------
#Generated from entering the Prey accession ID into "Workflow 1: Query proteins and retrieve profiles"at https://reprint-apms.org/?q=chooseworkflow
data_crapome<-read.csv("data/raw/Crapome_1716919084_gp.csv",  header = TRUE, check.names = FALSE) 

# Split the 'Num of Expt. (found/total)' column into two new columns
split_values <- strsplit(as.character(data_crapome$`Num of Expt. (found/total)`), " / ")

# Add the new columns to the data frame and convert to numeric
data_crapome$found <- as.numeric(sapply(split_values, `[`, 1))
data_crapome$total <- as.numeric(sapply(split_values, `[`, 2))

#calculate which of the Prey proteins that appear less than 50% of the crapome experiments 
data_crapome<- data_crapome %>%
  mutate(crapome_frequency=found/total)
#Merge the frequency data with the data_SAINT_nerged_nodified
#Remove the extra entries in the crapome file, these comes from presence of pulled down proteins in empty vector. 
data_SAINT_merged_modified<-merge(data_SAINT_merged_modified,data_crapome[,c('Mapped Gene Symbol',"crapome_frequency")], by.x="PreyGene", by.y ='Mapped Gene Symbol', all.x=TRUE, all.y=FALSE)

#If data is crapome or not 
data_SAINT_merged_modified = data_SAINT_merged_modified %>%
  mutate(is_not_crapome = crapome_frequency <= 0.5 | is.na(crapome_frequency))

#Plotting ----------------------------------------------------------------------------

##SP vs FC-A facet old ----
plot_SP_FCA <- ggplot(data = filter(data_SAINT_merged_modified), aes(x = SP, y = FCA, label = PreyGene)) +
  geom_hline(yintercept = 2, color="pink")+
  annotate("text", x = 0.1, y = 2.2 , label = "FCA = 2", color="pink")+
  geom_point(data=filter(data_SAINT_merged_modified, crapome_frequency<=0.5),aes(size=combined_score), color="black", show.legend = FALSE) +
  geom_point(data=filter(data_SAINT_merged_modified, crapome_frequency>0.5),aes(size=combined_score), color="grey", show.legend = FALSE) +
  geom_point(data=filter(data_SAINT_merged_modified, kinase==TRUE), color="red", shape=4, show.legend = FALSE) +
  
    
  # #hihglight sample names
  geom_text_repel(
    data = filter(
      data_SAINT_merged_modified,
      PreyGene == "MTHFR" | PreyGene == "MTHFD1" | PreyGene == "GCN1" | PreyGene == "PRKDC" | PreyGene == "FASN"
    ),
    aes(x = SP, y = FCA, label = PreyGene),
    size = 3.5,
    nudge_y = 0.2,
    nudge_x = -0.1,
    segment.size = 0.2,
    arrow = arrow(length = unit(0.01, "npc")),
    box.padding = 1,
    point.padding = 0.3,
    min.segment.length = 0,
    max.overlaps = 30,
    show.legend = FALSE
  ) +
  
  #Theme 
  facet_grid(~Bait)+
  theme_minimal() +
  ggtitle("SAINT analysis")+
  xlim(c(0,1))+
  scale_y_continuous(trans='log10')+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))

plot_SP_FCA

ggsave(
  "plot_SP_FCA_modified.png",
  plot = plot_SP_FCA,
  height = 2000,
  width = 4000,
  units = "px",
  path ="data/output"
  
)

(plot_SP_FCA
  + coord_cartesian(xlim = c(0.75,1))
  
  )

##SP vs FC-A facet ----


ggplot(data = data_SAINT_merged_modified,
       aes(
         x = (SP),
         y = FCA,
         label = PreyGene,
         shape = is_not_crapome,
         color = kinase
       )) +
  facet_grid(~ Bait) +
  geom_point(size = 1) +
  scale_color_manual(values = c("#fb8072", "#bebada")) +
  scale_shape_manual(values = c(4, 18)) +
  coord_cartesian(xlim = c(0.75, 1), ylim = c(1, 10)) +
  geom_text_repel(
    data = filter(pdat, FCA > 2 &
                    (is_not_crapome) & SP > 0.95 | PreyGene == 'MTHFR', ),
    aes(label = PreyGene),
    color = 'black',
    max.overlaps = 140,
    size = 3,

    box.padding = 0.5,

  ) +
  scale_y_log10() +
  theme_bw()

ggsave(
  "plot_SP_FCA_modified.png",
  height = 2000,
  width = 4000,
  units = "px",
  path ="data/output"
  
)

##SP vs FC-A facet Color  ----


# Reorder the data frame to highlight kinase, then no crapome and crapome in the back
data_SAINT_merged_modified <- data_SAINT_merged_modified %>%
  arrange(is_not_crapome) %>% # Change this to your desired ordering logic
  arrange(kinase)

data_SAINT_merged_modified <- data_SAINT_merged_modified %>%
  mutate(
    nudge_x = ifelse(PreyGene %in% c("MTHFR"), 0, -0.05),
    nudge_y = ifelse(PreyGene %in% c("MTHFR"), 0, 0.05)
  )
  
ggplot(data = data_SAINT_merged_modified,
       aes(
         x = SP,
         y = FCA,
         label = PreyGene,
         fill = is_not_crapome,
         color = kinase
       )) +
  facet_grid(~ Bait) +
  geom_point(size = 2,
             shape=21) +
  scale_color_manual(values = c("FALSE" = "white", "TRUE" = "red"),
                     labels = c("TRUE" = "Kinase", "FALSE"=" ")) +
  scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "gray"),
                    labels = c("TRUE" = "Not crapome", "FALSE" = "Crapome")) +
  coord_cartesian(xlim = c(0.75, 1), ylim = c(1, 10)) +
  geom_text_repel(
    data = filter(pdat, FCA > 2 &
                    (is_not_crapome) & SP > 0.95 | PreyGene == 'MTHFR' ),
    aes(label = PreyGene),
    color = 'black',
    max.overlaps = 140,
    size = 3,
    nudge_y = 0.05,
    nudge_x = -0.05,
    box.padding = 0.4,
    
  ) +
  scale_y_log10() +
  theme_bw()+
  labs(fill = NULL, color = NULL) +  # Remove legend titles
  guides(fill = guide_legend(title = NULL),  # Remove the fill legend title
         color = guide_legend(title = NULL))  # Remove the color legend title

ggsave(
  "plot_SP_FCA_modified_v2.png",
  height = 2000,
  width = 4000,
  units = "px",
  path ="data/output"
  
)


ggtitle("SAINT analysis")+
  xlim(c(0,1))+
  scale_y_continuous(trans='log10')+
  theme(plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))

plot_SP_FCA



ggsave(
  "plot_SP_FCA_modified.png",
  plot = plot_SP_FCA,
  height = 2000,
  width = 4000,
  units = "px",
  path ="data/output"
  
)







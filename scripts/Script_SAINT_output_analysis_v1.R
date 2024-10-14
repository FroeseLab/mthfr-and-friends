library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(DataCombine)
library(UniProt.ws)
library(grid)
library(gridExtra)
library("ggrepel")
library(UniprotR)
library(STRINGdb)
library(tidyverse)
library(scrutiny)
library(ggpubr)


#Read the input file
# Data files ----
getwd()

#setwd(getwd()) #sets file location as directory
setwd("~/Programs/Git/mthfr-and-friends")

#SAINT raw data---
#Data from SAINT analysis using "data/interim/SAINT_list_input.csv" as input, where empty vector is control
# SAINT analysis performed with https://reprint-apms.org/?q=analysis_front_apms online server 
#Import raw .tsv file of the "raw SAINT results". When downloaded a .output file is given. Just add .tsv file format
data_raw_SAINT <- read.table("data/raw/21262.output.tsv",
                             header = TRUE,
                             check.names = FALSE)
#Import the unfiltered matrix from "Bait-prey heatmap" tab. Downloaded as a .tab file. Just add .tsv file format
data_baitPrey <- read.table(
  "data/raw/21262_baitPrey_data.tab.tsv",
  header = TRUE,
  check.names = FALSE
)

#Merge SAINT files
#Note preys that are found in raw data but not found in bait-prey data will not be merged. E.g. A0A075B6S2
data_SAINT_merged <- merge(data_raw_SAINT, data_baitPrey)
# Save Merged SAINT data
write.csv(data_SAINT_merged, "data/output/SAINT_analysis_merged.csv")



#Plotting unmodified----
#SP vs FC-A all data [NOT IN ARTICLE]
plot_SP_FCA <- ggplot(data = data_SAINT_merged, aes(x = SP, y = FCA, label = PreyGene)) +
  geom_point(aes(color = Bait), show.legend = FALSE) +
  
  # #hihglight sample names
  geom_text_repel(
    data = filter(
      data_SAINT_merged,
      PreyGene == "MTHFR" |
        PreyGene == "MTHFD1" |
        PreyGene == "GCN1" | PreyGene == "PRKDC" | PreyGene == "FASN"
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
  ggtitle("SAINT analysis") +
  theme(plot.title = element_text(hjust = 0.5))

plot_SP_FCA

ggsave(
  "plot_SP_FCA_check.png",
  plot = plot_SP_FCA,
  height = 2000,
  width = 2000,
  units = "px",
  path = "data/output"
)




#Modify data ----
data_SAINT_merged_modified <- data_SAINT_merged



##Kinases -------------------------------------------------------------------------------

#Kinase
#raw data obtained (date: 26 July 2024) from https://amigo.geneontology.org/amigo/term/GO:0004672 with "Gene/product (bioentity_label)"
data_kinases <- read.csv(
  "data/raw/AmiGO_kinase_GO-0004672.csv",
  header = FALSE,
  check.names = FALSE
)
data_SAINT_merged_modified <- data_SAINT_merged_modified %>%
  mutate(kinase = data_SAINT_merged_modified$PreyGene %in% data_kinases$V1)


#Kinase activators
#raw data obtained (date: 3 June 2024) from https://amigo.geneontology.org/amigo/term/GO:0019207
data_kinases_activator <- read.csv(file = "data/raw/AmiGO_kinase_activator.csv",
                                   header = FALSE,
                                   check.names = FALSE)
data_SAINT_merged_modified <- data_SAINT_merged_modified %>%
  mutate(kinase_activator = data_SAINT_merged_modified$PreyGene %in% data_kinases_activator$V2)

# If you want to prioritize "kinase" over "kinase_activator" when both are TRUE
data_SAINT_merged_modified <- data_SAINT_merged_modified %>%
  mutate(
    merged_activity = case_when(
      kinase == TRUE ~ "kinase",
      kinase_activator == TRUE & kinase == FALSE ~ "kinase regulator",
      TRUE ~ "FALSE"
    )
  )




## STRING--------------------------------------------------------------------------------
#Setup, usally use thresholwd above 400
string_db <- STRINGdb$new(
  version = "12.0",
  species = 9606,
  score_threshold = 100,
  network_type = "full",
  input_directory = ""
)

#Gene names to look up in STRING
dat_example = data.frame(list(genes = unique(data_SAINT_merged_modified$PreyGene)))

#Retrieve STRING accession ID for each gene name
dat_examples_mapped = string_db$map(dat_example, 'genes')

data_string <- string_db$get_interactions(dat_examples_mapped$STRING_id) %>%
  merge(dat_examples_mapped %>%
          rename(from_gene = genes),
        by.x = 'from',
        by.y = 'STRING_id') %>%
  merge(dat_examples_mapped %>%
          rename(to_gene = genes),
        by.x = 'to',
        by.y = 'STRING_id')


#Get only interactions with MTHFR
data_STRING_MTHFR <- data_string %>%
  unique() %>%
  filter(from_gene == "MTHFR" | to_gene == "MTHFR") %>%
  dplyr::select(c(combined_score, from_gene, to_gene)) %>%
  # Reshape it to have the interacting gene
  tidyr::gather('direction', 'gene', -combined_score) %>%
  filter(gene != 'MTHFR')




#Merge to input combined scores and direction of interaction
data_SAINT_merged_modified <- merge(
  data_SAINT_merged_modified,
  data_STRING_MTHFR,
  by.x = "PreyGene",
  by.y = "gene",
  all = TRUE
)

#rename combined_score column
data_SAINT_merged_modified <- data_SAINT_merged_modified %>%
  rename(STRING_score = combined_score)


##Crapome----------------------------------------------------------------------------------
#Generated from entering the Prey accession IDs in file "data/interim/crapome_input.csv" 
#Insert into "Workflow 1: Query proteins and retrieve profiles"at https://reprint-apms.org/?q=chooseworkflow
#Date obtained 31st of May 2024
data_crapome <- read.csv("data/raw/Crapome_1716919084_gp.csv",
                         header = TRUE,
                         check.names = FALSE)

# Split the 'Num of Expt. (found/total)' column into two new columns
split_values <- strsplit(as.character(data_crapome$`Num of Expt. (found/total)`), " / ")

# Add the new columns to the data frame and convert to numeric
data_crapome$found <- as.numeric(sapply(split_values, `[`, 1))
data_crapome$total <- as.numeric(sapply(split_values, `[`, 2))

#calculate which of the Prey proteins that appear less than 50% of the crapome experiments
data_crapome <- data_crapome %>%
  mutate(crapome_frequency = found / total)
#Merge the frequency data with the data_SAINT_nerged_nodified
#Remove the extra entries in the crapome file, these comes from presence of pulled down proteins in empty vector.
data_SAINT_merged_modified <- merge(
  data_SAINT_merged_modified,
  data_crapome[, c('Mapped Gene Symbol', "crapome_frequency")],
  by.x = "PreyGene",
  by.y = 'Mapped Gene Symbol',
  all.x = TRUE,
  all.y = FALSE
)

#If data is crapome or not
data_SAINT_merged_modified = data_SAINT_merged_modified %>%
  mutate(is_not_crapome = crapome_frequency <= 0.5 |
           is.na(crapome_frequency))

##Reorder data ----
# Reorder the data frame to highlight kinase, then no crapome and crapome in the back
data_SAINT_merged_modified <- data_SAINT_merged_modified %>%
  arrange(is_not_crapome) %>% # Change this to your desired ordering logic
  arrange(merged_activity)



##Rename data ----
#Bait names with subscript
data_SAINT_merged_modified <- data_SAINT_merged_modified %>%
  mutate(
    Bait_name = case_when(
      Bait == "MTHFR" ~ "MTHFR[WT]",
      Bait == "MTHFR_T34A" ~ "MTHFR[T34A]",
      Bait == "MTHFR_38to656" ~ "MTHFR[38-656]",
      TRUE ~ Bait
    )
  )

#Reorder order of baits based on original sample name
data_SAINT_merged_modified$Bait = factor(data_SAINT_merged_modified$Bait,
                                         levels = c("MTHFR", "MTHFR_T34A", "MTHFR_38to656"))
#Reorder order of baits based on new name
data_SAINT_merged_modified$Bait_name = factor(
  data_SAINT_merged_modified$Bait_name,
  levels = c("MTHFR[WT]", "MTHFR[T34A]", "MTHFR[38-656]")
)



#Plotting ----------------------------------------------------------------------------

##Main Plot ----

# Temporarely to make crapome points go in the back
plot_main <- data_SAINT_merged_modified %>%
  arrange(is_not_crapome) %>% # Change this to your desired ordering logic
  
  
  ggplot(data = ,
         aes(
           x = SP,
           y = FCA,
           label = PreyGene,
           fill = is_not_crapome,
           color = merged_activity
         )) +
  
  
  geom_point(size = 2.5,
             shape = 21,
             stroke = 0.3) +
  
  #Limits for high confidence interactions
  geom_hline(yintercept = 2,
             color = "darkcyan",
             linetype = "dashed") +
  annotate(
    "text",
    x = 0.77,
    y = 2.1 ,
    label = "FCA = 2",
    color = "darkcyan"
  ) +
  geom_vline(xintercept = 0.95,
             color = "darkcyan",
             linetype = "dashed") +
  annotate(
    "text",
    x = 0.945,
    y = 8 ,
    angle = 90,
    label = "SP = 0.95",
    color = "darkcyan"
  ) +
  
  
  scale_color_manual(
    values = c(
      "FALSE" = "white",
      "kinase" = "white",
      "kinase regulator activity" = "white"
    ),
    labels = c("FALSE" = " ")
  ) +
  scale_fill_manual(
    values = c("TRUE" = "black", "FALSE" = "#B38A5E"),
    labels = c("TRUE" = "Not crapome", "FALSE" = "Crapome")
  ) +
  coord_cartesian(xlim = c(0.75, 1), ylim = c(1, 10)) +
  geom_text_repel(
    data = filter(
      data_SAINT_merged_modified,
      FCA > 2 &
        (is_not_crapome) &
        SP > 0.95 & PreyGene != 'MTHFR'
    ),
    aes(label = PreyGene),
    color = 'black',
    max.overlaps = 140,
    size = 5,
    #
    nudge_y = 0.07,
    nudge_x = -0.12,
    #
    box.padding = 0.4,
    segment.colour = "darkgray",
  ) +
  
  geom_text_repel(
    data = filter(data_SAINT_merged_modified, PreyGene == 'MTHFR'),
    aes(label = PreyGene),
    color = 'black',
    max.overlaps = 140,
    size = 5,
    #
    nudge_y = 0.05,
    nudge_x = -0.01,
    #
    box.padding = 0.4,
    segment.colour = "darkgray",
  ) +
  
  facet_grid( ~ Bait_name, labeller = label_parsed) + #label_parsed allows for the subscript to show in graph
  
  scale_y_log10() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 20),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15),
    
    
  ) +
  
  labs(fill = NULL, color = NULL) +  # Remove legend titles
  guides(fill = guide_legend(title = NULL), # Remove the fill legend title
         color = "none")  # Remove the color legend title

plot_main

#Have to separate the legend as this affects the size of the images, and the plots dont align with eachother
legend <- as_ggplot(cowplot::get_legend(plot_main))
grid <- grid.arrange(
  plot_main + theme(legend.position = "none"),
  legend,
  nrow = 1,
  widths = c(8, 1)
)


ggsave(
  "plot_SP_FCA_modified_main.png",
  grid,
  height = 10,
  width = 40,
  units = "cm",
  dpi = 600,
  path = "data/output"
  
)



##Kinases ----
#The protein CAD seems to be incorrectly labled with the kinase GO term and will therfore not be highlighted in the graph
data_SAINT_merged_modified[which(data_SAINT_merged_modified$PreyGene == "CAD"), "kinase"] <-
  FALSE
data_SAINT_merged_modified[which(data_SAINT_merged_modified$PreyGene == "CAD"), "merged_activity"] <-
  FALSE

plot_kinases <- ggplot(
  data = filter(data_SAINT_merged_modified, merged_activity != FALSE),
  aes(
    x = SP,
    y = FCA,
    label = PreyGene,
    fill = is_not_crapome,
    color = merged_activity,
    #     alpha = merged_activity
  )
) +
  
  
  
  
  
  
  geom_point(data = data_SAINT_merged_modified, aes(x = SP, y = FCA), color = "#E5E7E9") +
  geom_point(size = 2.5, shape = 21) +
  
  #Limits for high confidence interactions
  geom_hline(yintercept = 2,
             color = "darkcyan",
             linetype = "dashed") +
  annotate(
    "text",
    x = 0.77,
    y = 2.1 ,
    label = "FCA = 2",
    color = "darkcyan"
  ) +
  geom_vline(xintercept = 0.95,
             color = "darkcyan",
             linetype = "dashed") +
  annotate(
    "text",
    x = 0.945,
    y = 8 ,
    angle = 90,
    label = "SP = 0.95",
    color = "darkcyan"
  ) +
  
  
  scale_color_manual(
    values = c(
      "FALSE" = "white",
      "kinase" = "#D41159",
      "kinase regulator" = "#1A85FF"
    ),
    labels = c("FALSE" = " ")
  ) +
  scale_fill_manual(
    values = c("TRUE" = "black", "FALSE" = "#B38A5E"),
    labels = c("TRUE" = "Not crapome", "FALSE" = "Crapome")
  ) +
  #scale_alpha_manual ( values =c("FALSE" = 0.1), breaks = waiver(), na.value = NA)+
  coord_cartesian(xlim = c(0.75, 1), ylim = c(1, 10)) +
  geom_text_repel(
    data = filter(
      data_SAINT_merged_modified,
      FCA > 2  & SP > 0.95 & merged_activity != FALSE
    ),
    aes(label = PreyGene),
    color = 'black',
    max.overlaps = 140,
    size = 5,
    nudge_y = 0.05,
    nudge_x = -0.05,
    box.padding = 0.4,
    segment.colour = "darkgray",
  ) +
  
  facet_grid( ~ Bait_name, labeller = label_parsed) + #label_parsed allows for the subscript to show in graph
  
  
  scale_y_log10() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 20),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15)
    
  ) +
  labs(fill = NULL, color = NULL) +  # Remove legend titles
  guides(fill = guide_legend(title = NULL), # Remove the fill legend title
         color = guide_legend(title = NULL))  # Remove the color legend title

plot_kinases

#Have to separate the legend as this affects the size of the images, and the plots dont align with eachother
legend <- as_ggplot(cowplot::get_legend(plot_kinases))

grid <- grid.arrange(
  plot_kinases + theme(legend.position = "none"),
  legend,
  nrow = 1,
  widths = c(8, 1)
)

ggsave(
  "plot_SP_FCA_modified_kinase.png",
  grid,
  height = 10,
  width = 40,
  units = "cm",
  dpi = 600,
  path = "data/output"
  
)



##STRING  ----

# Temporarely
plot_STRING <- data_SAINT_merged_modified %>%
  arrange(STRING_score) %>% # Change this to your desired ordering logic
  arrange(merged_activity) %>%
  filter(STRING_score != FALSE) %>%
  
  ggplot(data = ,
         aes(
           x = SP,
           y = FCA,
           label = PreyGene,
           fill = STRING_score,
           color = merged_activity
         )) +
  
  
  geom_point(data = data_SAINT_merged_modified, aes(x = SP, y = FCA), color = "#E5E7E9") +
  geom_point(size = 2.5,
             shape = 21,
             stroke = 0.5) +
  
  #Limits for high confidence interactions
  geom_hline(yintercept = 2,
             color = "darkcyan",
             linetype = "dashed") +
  annotate(
    "text",
    x = 0.77,
    y = 2.1 ,
    label = "FCA = 2",
    color = "darkcyan"
  ) +
  geom_vline(xintercept = 0.95,
             color = "darkcyan",
             linetype = "dashed") +
  annotate(
    "text",
    x = 0.945,
    y = 8 ,
    angle = 90,
    label = "SP = 0.95",
    color = "darkcyan"
  ) +
  
  
  scale_color_manual(
    values = c(
      "FALSE" = "black",
      "kinase" = "black",
      "kinase regulator activity" = "black"
    ),
    labels = c("FALSE" = " ")
  ) +
  scale_fill_gradient(low = "gray",
                      high = "purple",
                      na.value = "black") +
  coord_cartesian(xlim = c(0.75, 1), ylim = c(1, 10)) +
  geom_text_repel(
    data = filter(
      data_SAINT_merged_modified,
      FCA > 2 &
        SP > 0.95 & STRING_score != FALSE
    ),
    aes(label = PreyGene),
    color = 'black',
    max.overlaps = 140,
    size = 5,
    nudge_y = 0.05,
    nudge_x = -0.05,
    box.padding = 0.4,
    segment.colour = "darkgray",
  ) +
  
  facet_grid( ~ Bait_name, labeller = label_parsed) + #label_parsed allows for the subscript to show in graph
  
  
  
  scale_y_log10() +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text.x = element_text(size = 20),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 15)
    
  ) +
  labs(fill = NULL, color = NULL) +  # Remove legend titles
  guides(fill = guide_legend(title = "STRING score"), color = "none")  # Remove the color legend title

plot_STRING

#Have to separate the legend as this affects the size of the images, and the plots dont align with eachother
legend <- as_ggplot(cowplot::get_legend(plot_STRING))
grid <- grid.arrange(
  plot_STRING + theme(legend.position = "none"),
  legend,
  nrow = 1,
  widths = c(8, 1)
)


ggsave(
  "plot_SP_FCA_modified_STRING.png",
  grid,
  height = 10,
  width = 40,
  units = "cm",
  dpi = 600,
  path = "data/output"
  
)

# grid arrange
# Arrange the plots using gridExtra
grid_plots <- grid.arrange(plot_main, plot_kinases, plot_STRING, ncol = 1)
grid_plots

##Crapome [NOT IN ARTICLE] ----


ggplot(
  data = data_SAINT_merged_modified,
  aes(
    x = SP,
    y = FCA,
    label = PreyGene,
    fill = crapome_frequency,
    color = merged_activity
  )
) +
  
  geom_hline(yintercept = 2,
             color = "darkcyan",
             linetype = "dashed") +
  annotate(
    "text",
    x = 0.77,
    y = 2.1 ,
    label = "FCA = 2",
    color = "darkcyan"
  ) +
  geom_vline(xintercept = 0.95,
             color = "darkcyan",
             linetype = "dashed") +
  annotate(
    "text",
    x = 0.94,
    y = 9 ,
    angle = 90,
    label = "SP = 0.95",
    color = "darkcyan"
  ) +
  
  
  facet_grid( ~ Bait) +
  geom_point(size = 2, shape = 21) +
  scale_color_manual(
    values = c(
      "FALSE" = "white",
      "kinase" = "red",
      "kinase regulator activity" = "green"
    ),
    labels = c("FALSE" = " ")
  ) +
  scale_fill_gradient(low = "green",
                      high = "red",
                      na.value = "black") +
  
  coord_cartesian(xlim = c(0.75, 1), ylim = c(1, 10)) +
  # geom_text_repel(
  #   data = filter(data_SAINT_merged_modified, kinase==TRUE),
  #   aes(label = PreyGene),
  #   color = 'black',
  #   max.overlaps = 140,
  #   size = 3,
  #   nudge_y = 0.05,
  #   nudge_x = -0.1,
  #   box.padding = 1,
  
  #) +
  scale_y_log10() +
  theme_bw() +
  labs(fill = NULL, color = NULL) +  # Remove legend titles
  guides(fill = guide_legend(title = NULL), # Remove the fill legend title
         color = guide_legend(title = NULL))  # Remove the color legend title

ggsave(
  "plot_SP_FCA_modified_crapome_gradient.png",
  height = 2000,
  width = 4000,
  units = "px",
  path = "data/output"
  
)


ggplot(
  data = filter(data_SAINT_merged_modified),
  aes(
    x = SP,
    y = FCA,
    label = PreyGene,
    fill = is_not_crapome,
    color = merged_activity
  )
) +
  
  facet_grid( ~ Bait) +
  geom_point(size = 2, shape = 21) +
  
  #Limits for high confidence interactions
  geom_hline(yintercept = 2,
             color = "darkcyan",
             linetype = "dashed") +
  annotate(
    "text",
    x = 0.77,
    y = 2.1 ,
    label = "FCA = 2",
    color = "darkcyan"
  ) +
  geom_vline(xintercept = 0.95,
             color = "darkcyan",
             linetype = "dashed") +
  annotate(
    "text",
    x = 0.945,
    y = 9 ,
    angle = 90,
    label = "SP = 0.95",
    color = "darkcyan"
  ) +
  
  
  scale_color_manual(
    values = c(
      "FALSE" = "white",
      "kinase" = "red",
      "kinase regulator activity" = "green"
    ),
    labels = c("FALSE" = " ")
  ) +
  scale_fill_manual(
    values = c("TRUE" = "black", "FALSE" = "gray"),
    labels = c("TRUE" = "Not crapome", "FALSE" = "Crapome")
  ) +
  coord_cartesian(xlim = c(0.75, 1), ylim = c(1, 10)) +
  
  scale_y_log10() +
  theme_bw() +
  labs(fill = NULL, color = NULL) +  # Remove legend titles
  guides(fill = guide_legend(title = NULL), # Remove the fill legend title
         color = guide_legend(title = NULL))  # Remove the color legend title

ggsave(
  "plot_SP_FCA_modified_crapome.png",
  height = 15,
  width = 40,
  units = "cm",
  dpi = 600,
  path = "data/output"
  
)

##Heatmap [NOT IN ARTICLE]----
#filter out prey genes that are significant according to parameters
data_significant_Prey<-data_SAINT_merged_modified %>%
  filter( FCA>=2, SP > 0.95)
#single out names  
data_significant_Prey<-unique(data_significant_Prey$Prey)

#Set up plot settings
plot_heatmap <- ggplot(data = filter(data_SAINT_merged_modified, Prey %in% data_significant_Prey) ,
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
    axis.title.y = element_blank(),
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
ggsave("plot_heatmaps_modified.png",plot=plot_heatmaps,  height = 2000, width=3000, units = "px", path ="data/output" )



#Retrieve protein names----
#Generate dataframe and retrieve protein names from Uniprot
data_protein_name <- data.frame(Prey = unique(data_SAINT_merged$Prey))
data_protein_name$protein_name <- GetProteinAnnontate(unique(data_SAINT_merged$Prey), c("protein_name"))
#merge data
data_SAINT_merged_modified <- merge(data_SAINT_merged_modified, data_protein_name)


#Save data file ----
# Save Merged SAINT data
write.csv(data_SAINT_merged_modified,
          "data/output/SAINT_analysis_merged_modified.csv")

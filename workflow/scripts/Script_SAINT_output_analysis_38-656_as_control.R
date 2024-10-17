library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(DataCombine)
library(UniProt.ws)
library(grid)
library(gridExtra)
library("ggrepel")
library(UniProt.ws)
library(UniprotR)
library(STRINGdb)
library(tidyverse)
library(scrutiny)
library(ggpubr)


#Read the input file
# Data files ----
getwd()

#SAINT raw data---
#Data from SAINT analysis using "data/interim/crapome_input.csv" as input, where MTHFR38-656 is used as a control. 
# SAINT analysis performed with https://reprint-apms.org/?q=analysis_front_apms online server 
#Import raw .tsv file of the "raw SAINT results". When downloaded a .output file is given. Just add .tsv file format
data_raw_SAINT <- read.table("data/raw/21302.output.tsv",
                             header = TRUE,
                             check.names = FALSE)
#Import the unfiltered matrix from "Bait-prey heatmap" tab. Downloaded as a .tab file. Just add .tsv file format
data_baitPrey <- read.table(
  "data/raw/21302_baitPrey_data.tab.tsv",
  header = TRUE,
  check.names = FALSE
)

#Merge SAINT files
#Note preys that are found in raw data but not found in bait-prey data will not be merged. E.g. A0A075B6S2
data_SAINT_merged <- merge(data_raw_SAINT, data_baitPrey)


#Modify data ----
data_SAINT_merged_modified <- data_SAINT_merged

##Crapome----------------------------------------------------------------------------------
#Generated from entering the Prey accession ID into "Workflow 1: Query proteins and retrieve profiles"at https://reprint-apms.org/?q=chooseworkflow
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




##Rename ----
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

#Reorder order of baits
data_SAINT_merged_modified$Bait = factor(data_SAINT_merged_modified$Bait,
                                         levels = c("MTHFR", "MTHFR_T34A", "MTHFR_38to656"))

data_SAINT_merged_modified$Bait_name = factor(
  data_SAINT_merged_modified$Bait_name,
  levels = c("MTHFR[WT]", "MTHFR[T34A]", "MTHFR[38-656]")
)




#Plotting ----------------------------------------------------------------------------

##Main Plot ----

# Temporarely to make crapome points go in the back
plot_main <- data_SAINT_merged_modified %>%
  arrange(is_not_crapome) %>% # Change this to your desired ordering logic
  
  
  ggplot(data = , aes(
    x = SP,
    y = FCA,
    label = PreyGene,
    fill = is_not_crapome
  )) +
  
  
  geom_point(
    size = 2.5,
    shape = 21,
    stroke = 0.3,
    color = "white"
  ) +
  
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
    x = 0.94,
    y = 3 ,
    angle = 90,
    label = "SP = 0.95",
    color = "darkcyan"
  ) +
  
  
  scale_fill_manual(
    values = c("TRUE" = "black", "FALSE" = "#B38A5E"),
    labels = c("TRUE" = "Not crapome", "FALSE" = "Crapome")
  ) +
  
  
  
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
  
  xlim(0, 1) +
  ylim(0, 10) +
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
  guides(fill = guide_legend(title = NULL))  # Remove the color legend title

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
  "plot_SP_FCA_modified_main_38-656_as_control.png",
  grid,
  height = 10,
  width = 40,
  units = "cm",
  dpi = 600,
  path = "data/output"
  
)

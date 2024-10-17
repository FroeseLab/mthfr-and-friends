if (exists("snakemake")) {
  # input
  fn_saint_merged <- snakemake@input$saint_merged
  fn_crapome <- snakemake@input$crapome
  fn_amigo_kinase <- snakemake@input$amigo_kinase
  fn_amigo_kinase_activator <- snakemake@input$amigo_kinase_activator
  fol_string_db <- snakemake@input$string_db
  # output
  fn_plot_sp_fca <- snakemake@output$plot_sp_fca
  fn_merged_saint_output <- snakemake@output$merged_saint_output
  fn_merged_saint_output_rds <- snakemake@output$merged_saint_output_rds
} else {
  fn_saint_merged <- "data/output/saint_output_merged_main.csv"
  fn_amigo_kinase <- "data/interim/AmiGO_kinase_GO-0004672.csv"
  fn_amigo_kinase_activator <- "data/interim/AmiGO_kinaseactivator_GO-0019207.csv"
  fn_crapome <- "data/interim/crapome_output.csv"
  fn_plot_sp_fca <- snakemake@output$plot_sp_fca
  fol_string_db <- "data/interim/string_db"
  # output
  fn_merged_saint_output <- "results/tmp/merged_saint_output.csv"
  fn_merged_saint_output_rds <- "results/tmp/merged_saint_output.rds"
}

# Read the input file
library(ggplot2)
library(ggrepel)
library(dplyr)
library(readr)
library(stringr)
library(UniProt.ws)
library(UniprotR)
library(STRINGdb)
library(tidyverse)


# SAINT raw data---
data_SAINT_merged <- read.csv(fn_saint_merged)



# Plotting unmodified----
# SP vs FC-A all data [NOT IN ARTICLE]
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
  fn_plot_sp_fca,
  plot = plot_SP_FCA,
  height = 2000,
  width = 2000,
  units = "px",
)




# Modify data ----
data_SAINT_merged_modified <- data_SAINT_merged



## Kinases -------------------------------------------------------------------------------

# Kinase
# raw data obtained (date: 26 July 2024) from https://amigo.geneontology.org/amigo/term/GO:0004672 with "Gene/product (bioentity_label)"
data_kinases <- read.csv(
  fn_amigo_kinase,
  header = FALSE,
  check.names = FALSE
)
data_SAINT_merged_modified <- data_SAINT_merged_modified %>%
  mutate(kinase = data_SAINT_merged_modified$PreyGene %in% data_kinases$V1)


# Kinase activators
# raw data obtained (date: 3 June 2024) from https://amigo.geneontology.org/amigo/term/GO:0019207
data_kinases_activator <- read.csv(
  file = fn_amigo_kinase_activator,
  header = FALSE,
  check.names = FALSE
)
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
# Setup, usally use thresholwd above 400
string_db <- STRINGdb$new(
  version = "12.0",
  species = 9606,
  score_threshold = 100,
  network_type = "full",
  input_directory = fol_string_db
)

# Gene names to look up in STRING
dat_example <- data.frame(list(genes = unique(data_SAINT_merged_modified$PreyGene)))

# Retrieve STRING accession ID for each gene name
dat_examples_mapped <- string_db$map(dat_example, "genes")

data_string <- string_db$get_interactions(dat_examples_mapped$STRING_id) %>%
  merge(
    dat_examples_mapped %>%
      rename(from_gene = genes),
    by.x = "from",
    by.y = "STRING_id"
  ) %>%
  merge(
    dat_examples_mapped %>%
      rename(to_gene = genes),
    by.x = "to",
    by.y = "STRING_id"
  )


# Get only interactions with MTHFR
data_STRING_MTHFR <- data_string %>%
  unique() %>%
  filter(from_gene == "MTHFR" | to_gene == "MTHFR") %>%
  dplyr::select(c(combined_score, from_gene, to_gene)) %>%
  # Reshape it to have the interacting gene
  tidyr::gather("direction", "gene", -combined_score) %>%
  filter(gene != "MTHFR")




# Merge to input combined scores and direction of interaction
data_SAINT_merged_modified <- merge(
  data_SAINT_merged_modified,
  data_STRING_MTHFR,
  by.x = "PreyGene",
  by.y = "gene",
  all = TRUE
)

# rename combined_score column
data_SAINT_merged_modified <- data_SAINT_merged_modified %>%
  rename(STRING_score = combined_score)


## Crapome----------------------------------------------------------------------------------
# Generated from entering the Prey accession IDs in file "data/interim/crapome_input.csv"
# Insert into "Workflow 1: Query proteins and retrieve profiles"at https://reprint-apms.org/?q=chooseworkflow
# Date obtained 31st of May 2024
data_crapome <- read.csv(fn_crapome,
  header = TRUE,
  check.names = FALSE
)

# Split the 'Num of Expt. (found/total)' column into two new columns
split_values <- strsplit(as.character(data_crapome$`Num of Expt. (found/total)`), " / ")

# Add the new columns to the data frame and convert to numeric
data_crapome$found <- as.numeric(sapply(split_values, `[`, 1))
data_crapome$total <- as.numeric(sapply(split_values, `[`, 2))

# calculate which of the Prey proteins that appear less than 50% of the crapome experiments
data_crapome <- data_crapome %>%
  mutate(crapome_frequency = found / total)
# Merge the frequency data with the data_SAINT_nerged_nodified
# Remove the extra entries in the crapome file, these comes from presence of pulled down proteins in empty vector.
data_SAINT_merged_modified <- merge(
  data_SAINT_merged_modified,
  data_crapome[, c("Mapped Gene Symbol", "crapome_frequency")],
  by.x = "PreyGene",
  by.y = "Mapped Gene Symbol",
  all.x = TRUE,
  all.y = FALSE
)

# If data is crapome or not
data_SAINT_merged_modified <- data_SAINT_merged_modified %>%
  mutate(is_not_crapome = crapome_frequency <= 0.5 |
    is.na(crapome_frequency))

## Reorder data ----
# Reorder the data frame to highlight kinase, then no crapome and crapome in the back
data_SAINT_merged_modified <- data_SAINT_merged_modified %>%
  arrange(is_not_crapome) %>% # Change this to your desired ordering logic
  arrange(merged_activity)



## Rename data ----
# Bait names with subscript
data_SAINT_merged_modified <- data_SAINT_merged_modified %>%
  mutate(
    Bait_name = case_when(
      Bait == "MTHFR" ~ "MTHFR[WT]",
      Bait == "MTHFR_T34A" ~ "MTHFR[T34A]",
      Bait == "MTHFR_38to656" ~ "MTHFR[38-656]",
      TRUE ~ Bait
    )
  )

# Reorder order of baits based on original sample name
data_SAINT_merged_modified$Bait <- factor(data_SAINT_merged_modified$Bait,
  levels = c("MTHFR", "MTHFR_T34A", "MTHFR_38to656")
)
# Reorder order of baits based on new name
data_SAINT_merged_modified$Bait_name <- factor(
  data_SAINT_merged_modified$Bait_name,
  levels = c("MTHFR[WT]", "MTHFR[T34A]", "MTHFR[38-656]")
)


# Retrieve protein names----
# Generate dataframe and retrieve protein names from Uniprot
data_protein_name <- data.frame(Prey = unique(data_SAINT_merged$Prey))
data_protein_name$protein_name <- GetProteinAnnontate(unique(data_SAINT_merged$Prey), c("protein_name"))
# merge data
data_SAINT_merged_modified <- merge(data_SAINT_merged_modified,
  data_protein_name,
  on = "Prey", all.x = TRUE
)


# Save data file ----
# Save Merged SAINT data
write.csv(
  data_SAINT_merged_modified,
  fn_merged_saint_output
)

# write out as rds as well
saveRDS(data_SAINT_merged_modified, fn_merged_saint_output_rds)

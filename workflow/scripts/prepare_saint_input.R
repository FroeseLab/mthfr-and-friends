# Read the input file ----

if (exists("snakemake")) {
  # input
  fn_raw_scaffold <- snakemake@input$data_raw_scaffold
  fn_data_rename <- snakemake@input$data_rename

  # output
  fn_saint <- snakemake@output$data_saint
  fn_saint_mthfr38to65 <- snakemake@output$data_saint_MTHFR38to65
  fn_crapome <- snakemake@output$crapome
} else {
  fn_raw_scaffold <- "data/raw/Prepare_SAINT_list_input_file/Proteins Report of merged human files with all 4 conditions (1. Feb, 5. Feb, 18. Feb).txt"
  fn_data_rename <- "data/metadata/SAINT_rename_file.csv"

  fn_saint <- "data/interim/SAINT_list_input.csv"
  fn_saint_mthfr38to65 <- "data/interim/SAINT_list_input_MTHFR38to656_control.csv"
  fn_crapome <- "data/interim/crapome_input.csv"
}

library(dplyr)

# Import raw csv file. Manually converted from output Scaffold .txt file prior.
data_raw_scaffold <- read.csv(fn_raw_scaffold,
  skip = 1, header = TRUE, check.names = FALSE, sep = "\t"
) # it is just the E1 SYPRO 5x from the original raw data from 220420
# Import csv file for renaming of Bait and AP-names
data_rename <- read.csv(fn_data_rename, header = TRUE, check.names = FALSE) # it is just the E1 SYPRO 5x from the original raw data from 220420

# Modify file ----

# Rename the columns
data_saint_input <- data_raw_scaffold %>%
  rename(
    `#BaitName_Condition` = `Biological sample category`,
    `#APName` = `Biological sample name`,
    `#PreyName` = `Protein accession number`,
    `#SpectralCounts` = `Number of total spectra`
  )

# Select only the renamed columns , select suddenly did not work therefore using dplyr::select
data_saint_input <- data_saint_input %>% dplyr::select(`#BaitName_Condition`, `#APName`, `#PreyName`, `#SpectralCounts`)

# Remove duplicate rows
data_saint_input <- data_saint_input %>% distinct()

#  Remove all rows where column "Number of total spectra" == 0
data_saint_input <- data_saint_input %>% filter(`#SpectralCounts` != 0)

# Remove all rows with random AP names, e.g. the column "Protein accession number" contains "|"
data_saint_input <- data_saint_input %>% filter(!grepl("\\|", `#PreyName`))

# Rename entries in #BaitName_Condition and #APName FindReplace
# requirements for SAINT input file in list format following https://reprint-apms.org/?q=inputfileformatting
# Note, the CONTROL should not contain any "_" , the Bait names can not contain any "-"
bait_name_map <- data_rename[, "#BaitName_Condition_new"]
names(bait_name_map) <- data_rename[, "#BaitName_Condition_old"]
ap_name_map <- data_rename[, "#APName_new"]
names(ap_name_map) <- data_rename[, "#APName_old"]

data_saint_input <- data_saint_input %>%
  mutate(
    `#BaitName_Condition` = bait_name_map[`#BaitName_Condition`],
    `#APName` = ap_name_map[`#APName`]
  )

## Save data ----
# Save the cleaned data to a new CSV file
write.csv(data_saint_input, fn_saint, row.names = FALSE, quote = FALSE)

# SAINT input file with MTHFR38to656 as control ----
# Modified input file where Empty vector is removed as a control
data_saint_input_control <- data_saint_input %>%
  filter(`#BaitName_Condition` != "CONTROL")

# Re-name MTHFR to CONTROL
data_saint_input_control$`#BaitName_Condition` <- replace(
  data_saint_input_control$`#BaitName_Condition`,
  data_saint_input_control$`#BaitName_Condition` == "MTHFR_38to656", "CONTROL"
)
## Save data ----
# Save the cleaned data to a new CSV file
write.csv(data_saint_input_control, fn_saint_mthfr38to65, row.names = FALSE, quote = FALSE)


# Crapome list ----
# For manual retrieval of the Prey appearing in the crapome
# Input in "Workflow 1: Query proteins and retrieve profiles"at https://reprint-apms.org/?q=chooseworkflow
data_crapome <- unique(data_saint_input$`#PreyName`)

## Save data ----
# Save the cleaned data to a new CSV file
write.table(data_crapome, file = fn_crapome, row.names = FALSE, col.names = FALSE, quote = FALSE)

# Combine the saint outputs

if (exists("snakemake")) {
  # input
  fn_raw_saint <- snakemake@input$saint_output
  fn_baitprey <- snakemake@input$saint_prey_output

  # output
  fn_saint_merged <- snakemake@output$saint_output
} else {
  fn_raw_saint <- "data/resources/prepared/21262.output.tsv"
  fn_baitprey <- "data/resources/prepared/21262_baitPrey_data.tab.tsv"

  fn_saint_merged <- "data/output/saint_analysis_merged.csv"
}

# saint raw data---
# Data from saint analysis using "data/interim/saint_list_input.csv" as input, where empty vector is control
# saint analysis performed with https://reprint-apms.org/?q=analysis_front_apms online server
# Import raw .tsv file of the "raw saint results". When downloaded a .output file is given. Just add .tsv file format
data_raw_saint <- read.table(fn_raw_saint,
  header = TRUE,
  check.names = FALSE
)
# Import the unfiltered matrix from "Bait-prey heatmap" tab. Downloaded as a .tab file. Just add .tsv file format
data_baitprey <- read.table(
  fn_baitprey,
  header = TRUE,
  check.names = FALSE
)

# Merge saint files
# Note preys that are found in raw data but not found in bait-prey data will not be merged. E.g. A0A075B6S2
data_saint_merged <- merge(data_raw_saint, data_baitprey, all.x = TRUE)

# Save Merged saint data
write.csv(data_saint_merged, fn_saint_merged)

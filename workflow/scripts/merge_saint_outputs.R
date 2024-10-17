# Combine the saint outputs

if (exists("snakemake")) {
  # input
  fn_raw_saint <- snakemake@input$saint_output
  fn_baitprey <- snakemake@input$saint_prey_output

  # output
  fn_saint_merged <- snakemake@output$saint_output
} else {
  fn_raw_saint <- "resources/precomputed/21262.output.tsv"
  fn_baitprey <- "resources/precomputed/21262_baitPrey_data.tab.tsv"

  fn_saint_merged <- "data/output/saint_analysis_merged.csv"
}

# saint raw data---
# Data from saint analysis using "data/interim/saint_list_input.csv" as input, where empty vector is control
# saint analysis performed with https://reprint-apms.org/?q=analysis_front_apms online server
# Import raw .tsv file of the "raw saint results". When downloaded a .output file is given. Just add .tsv file format

#' Read SAINT output file
#'
#' The saint output file repeats the header
#' multiple times, thus this function filters
#' out the repeated headers.
#' @param fn File path to SAINT output file
#'
#' @return Data frame with SAINT output
#'
read_saint_table <- function(fn) {
  lines <- readLines(fn)
  want <- !grepl("^Bait", lines)
  want[1] <- TRUE
  read.table(textConnection(lines[want]), sep = "\t", header = TRUE)
}
data_raw_saint <- read_saint_table(fn_raw_saint)
# Import the unfiltered matrix from "Bait-prey heatmap" tab. Downloaded as a .tab file. Just add .tsv file format
data_baitprey <- read.table(
  fn_baitprey,
  header = TRUE,
  check.names = FALSE
)

# Merge saint files
# Note preys that are found in raw data but not found in bait-prey data will not be merged. E.g. A0A075B6S2
data_saint_merged <- merge(data_raw_saint, data_baitprey,
  all.x = TRUE,
  on = c("Prey", "Bait")
)

# Save Merged saint data
write.csv(data_saint_merged, fn_saint_merged, row.names = FALSE)

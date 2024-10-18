if (exists("snakemake")) {
  # input
  fn_saint_input <- snakemake@input$saint_input
  fn_saint_annotated <- snakemake@input$saint_annotated
  # output
  fn_upsetr_input <- snakemake@output$upsetr_input
  fn_upsetr_output <- snakemake@output$upsetr_output
  fn_upsetr_filtered <- snakemake@output$upsetr_filtered
} else {
  fn_saint_input <- "data/interim/SAINT_list_input.csv"
  fn_saint_annotated <- "data/output/saint_output_merged_annotated_main.rds"
  # output
  fn_upsetr_input <- "results/tmp/upset_plot_SAINT_input.png"
  fn_upsetr_output <- "results/tmp/upset_plot_SAINT_output.png"
  fn_upsetr_filtered <- "results/tmp/upset_SAINT_output_SP_0.95_FCA_2.png"
}
library(UpSetR)
library(tidyr)
library(dplyr)
library(ggplot2)

# Prepare raw Input data ------------------------------------------------------------
# Import raw csv file. Manually converted from output Scaffold .txt file prior.


data_raw_SAINT_input <- read.csv(fn_saint_input, header = TRUE, check.names = FALSE) # it is just the E1 SYPRO 5x from the original raw data from 220420

# Re-name
data_raw_SAINT_input <- data_raw_SAINT_input %>%
  mutate(`#BaitName_Condition` = case_when(
    `#BaitName_Condition` == "MTHFR_38to656" ~ "MTHFR_38-656",
    `#BaitName_Condition` == "MTHFR" ~ "MTHFR_WT",
    TRUE ~ `#BaitName_Condition`
  ))



listInput <- data_raw_SAINT_input %>%
  split(.$`#BaitName_Condition`) %>%
  lapply(function(x) x$`#PreyName`)
length(listInput)


### Generate Diagram of raw data

upset_plot_SAINT_input <- upset(
  fromList(listInput),
  order.by = "freq",
  text.scale = c(1.5, 2, 1.5, 1.5, 1.5, 1.5) * 2,
  # Parameters order: c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
  point.size = 8,
  line.size = 3,
  sets.x.label = "Set Size",
  # Adjust this to change the size of the dots in the matrix
  sets.bar.color = "darkgray",
  # Customize as needed
  main.bar.color = "darkgray",
  # Customize as needed
  matrix.color = "black"
) # Customize as needed)
upset_plot_SAINT_input



# Open a PNG device
png(
  filename = fn_upsetr_input, height = 8,
  width = 13,
  units = "in", res = 300
)
upset_plot_SAINT_input
# Close the device
dev.off()

### Generate Diagram of SAINT output data \[not in article\]


# Prepare SAINT Input data ------------------------------------------------------------
data_raw_SAINT_output <- readRDS(fn_saint_annotated)



listInput_output <- data_raw_SAINT_output %>%
  split(.$`Bait`) %>%
  lapply(function(x) x$Prey)

# Remove the Bait element from the list
listInput_output <- listInput_output[!names(listInput_output) %in% "Bait"]


# Generate Diagaram of raw data
upset_plot_SAINT_output <- upset(
  fromList(listInput_output),
  order.by = "freq",
  text.scale = c(1.5, 2, 1.5, 1.5, 1.5, 1.5) * 2,
  # Parameters order: c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
  point.size = 8,
  line.size = 3,
  sets.x.label = "Set Size",
  # Adjust this to change the size of the dots in the matrix
  sets.bar.color = "gray",
  # Customize as needed
  main.bar.color = "darkgray",
  # Customize as needed
  matrix.color = "black"
) # Customize as needed)
upset_plot_SAINT_output

# Open a PNG device
png(
  filename = fn_upsetr_output, height = 8,
  width = 13,
  units = "in", res = 300
)
upset_plot_SAINT_output
# Close the device
dev.off()


### Generate Diagram of SAINT output data filtered


# Prepare SAINT Input data ------------------------------------------------------------



# Filter SAINT output data based onf F

listInput_output <- data_raw_SAINT_output %>%
  filter(SP > 0.95 & FCA > 2) %>%
  split(.$`Bait`) %>%
  lapply(function(x) x$Prey)




# Remove the Bait element from the list
# listInput_output <- listInput_output[!names(listInput_output) %in% "Bait"]


upset_plot_SAINT_output <- upset(
  fromList(listInput_output),
  order.by = "freq",
  text.scale = c(1.5, 2, 1.5, 1.5, 1.5, 1.5) * 2,
  # Parameters order: c(intersection size title, intersection size tick labels, set size title, set size tick labels, set names, numbers above bars)
  point.size = 8,
  line.size = 3,
  sets.x.label = "Set Size",
  # Adjust this to change the size of the dots in the matrix
  sets.bar.color = "gray", # Customize as needed
  main.bar.color = "darkgray", # Customize as needed
  matrix.color = "black"
) # Customize as needed)
upset_plot_SAINT_output


# Generate Diagaram of raw data upset_plot_SAINT_output \<- upset( fromList(listInput_output), order.by = "freq", text.scale = c(2, 2, 2, 1.5, 2, 1.5)\*2 , point.size = 10, line.size = 5, sets.x.label = "Set Size", \# Adjust this to change the size of the dots in the matrix sets.bar.color = "gray", \# Customize as needed main.bar.color = "darkgray", \# Customize as needed matrix.color = "black" ) \# Customize as needed) upset_plot_SAINT_output


# Open a PNG device
png(
  filename = fn_upsetr_filtered, height = 8,
  width = 13,
  units = "in", res = 300
)
upset_plot_SAINT_output
# Close the device
dev.off()

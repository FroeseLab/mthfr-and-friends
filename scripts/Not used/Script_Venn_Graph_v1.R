library(UpSetR)
library(tidyr)
library(dplyr)
library(ggplot2)
#library(ggupset)
#Read th#Read the input file 
# Data files ----
getwd()

#Prepare raw Input data ------------------------------------------------------------
#Import raw csv file. Manually converted from output Scaffold .txt file prior. 
data_raw_SAINT_input<- read.csv("data/interim/SAINT_list_input.csv",  header = TRUE, check.names = FALSE) #it is just the E1 SYPRO 5x from the original raw data from 220420



listInput = data_raw_SAINT_input %>%
  split(.$`#BaitName_Condition`) %>%
  lapply(function(x) x$`#PreyName`)



#Generate Diagaram of raw data
upset_plot_SAINT_input <- upset(
  fromList(listInput),
  order.by = "freq",
  text.scale = c(2, 2, 2, 1.5, 2, 1.5) * 5,
  point.size = 10,
  line.size = 5,
  sets.x.label = "Set Size",
  # Adjust this to change the size of the dots in the matrix
  sets.bar.color = "gray",
  # Customize as needed
  main.bar.color = "darkgray",
  # Customize as needed
  matrix.color = "black"
) # Customize as needed)
upset_plot_SAINT_input

# Open a PNG device
png(filename = "data/output/upset_plot_SAINT_input.png",   height = 2000,
    width = 4000,
    units = "px")
upset_plot_SAINT_input
# Close the device
dev.off()

#Prepare SAINT Input data ------------------------------------------------------------
data_raw_SAINT_output<- read.table("data/raw/21262.output.tsv",  header = TRUE, check.names = FALSE) 

listInput_output = data_raw_SAINT_output %>%
  split(.$`Bait`) %>%
  lapply(function(x) x$Prey)

# Remove the Bait element from the list
listInput_output <- listInput_output[!names(listInput_output) %in% "Bait"]




#Generate Diagaram of raw data 
upset_plot_SAINT_output <- upset(
  fromList(listInput_output),
  order.by = "freq",
  text.scale = c(2, 2, 2, 1.5, 2, 1.5) * 5,
  point.size = 10,
  line.size = 5,
  sets.x.label = "Set Size",
  # Adjust this to change the size of the dots in the matrix
  sets.bar.color = "gray",
  # Customize as needed
  main.bar.color = "darkgray",
  # Customize as needed
  matrix.color = "black"
) # Customize as needed))
upset_plot_SAINT_output


# Open a PNG device
png(filename = "data/output/upset_SAINT_output.png",  height = 2000,
    width = 4000,
    units = "px")
upset_plot_SAINT_output
# Close the device
dev.off()

#Filter SAINT output data based on FDR

listInput_output = data_raw_SAINT_output %>%
  filter(FDR<=0.1) %>%
  split(.$`Bait`) %>%
  lapply(function(x) x$Prey)

# Remove the Bait element from the list
listInput_output <- listInput_output[!names(listInput_output) %in% "Bait"]


#Generate Diagaram of raw data 
upset_plot_SAINT_output <- upset(
  fromList(listInput_output),
  order.by = "freq",
  text.scale = c(2, 2, 2, 1.5, 2, 1.5) * 5,
  point.size = 10,
  line.size = 5,
  sets.x.label = "Set Size",
  # Adjust this to change the size of the dots in the matrix
  sets.bar.color = "gray",
  # Customize as needed
  main.bar.color = "darkgray",
  # Customize as needed
  matrix.color = "black"
) # Customize as needed))
upset_plot_SAINT_output


# Open a PNG device
png(filename = "data/output/upset_SAINT_output_FDR_0.1.png",  height = 2000,
    width = 4000,
    units = "px")
upset_plot_SAINT_output
# Close the device
dev.off()


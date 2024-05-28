library(dplyr)
library(DataCombine)

#Read the input file ----
getwd()
# Set directory 
setwd(getwd())

#Import csv file 
data_raw_scaffold<- read.csv("Proteins_Report_raw.csv", skip = 1, header = TRUE, check.names = FALSE) #it is just the E1 SYPRO 5x from the original raw data from 220420
#Import csv file 
data_rename<- read.csv("rename_file.csv", header = TRUE, check.names = FALSE) #it is just the E1 SYPRO 5x from the original raw data from 220420

#Modify file ----

# Rename the columns
data_SAINT_input <- data_raw_scaffold %>%
  rename(
    `#BaitName_Condition` = `Biological sample category`,
    `#APName` = `Biological sample name`,
    `#PreyName` = `Protein accession number`,
    `#SpectralCounts` = `Number of total spectra`
  )

# Select only the renamed columns
data_SAINT_input <- data_SAINT_input %>% select(`#BaitName_Condition`, `#APName`, `#PreyName`, `#SpectralCounts`)

# Remove duplicate rows
data_SAINT_input <- data_SAINT_input %>% distinct()

#  Remove all rows where column "Number of total spectra" == 0
data_SAINT_input <- data_SAINT_input %>% filter(`#SpectralCounts` != 0)

# Remove all rows with random AP names, e.g. the column "Protein accession number" contains "|"
data_SAINT_input <- data_SAINT_input %>% filter(!grepl("\\|", `#PreyName`))

# Rename entries in #BaitName_Condition and #APName FindReplace
#requirements for SAINT input file in list format following https://reprint-apms.org/?q=inputfileformatting
#Note, the control should not contain any "_" 
data_SAINT_input<-FindReplace(data_SAINT_input, "#BaitName_Condition", data_rename, from = "#BaitName_Condition_old", to = "#BaitName_Condition_new", exact = TRUE, vector = FALSE)
data_SAINT_input<-FindReplace(data_SAINT_input, "#APName", data_rename, from = "#APName_old", to = "#APName_new", exact = TRUE, vector = FALSE)

# Save data ----
# Save the cleaned data to a new CSV file
write.csv(data_SAINT_input, "SAINT_list_input.csv", row.names = FALSE, quote = FALSE)


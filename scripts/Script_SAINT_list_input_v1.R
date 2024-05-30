library(dplyr)
library(DataCombine)
library(UniProt.ws)

#Read the input file ----
getwd() #Detects file location
# Set directory 
#setwd(getwd()) #sets file location as directory 
setwd("~/Programs/Git/mthfr-and-friends")

#Import raw csv file. Manually converted from output Scaffold .txt file prior. 
data_raw_scaffold<- read.csv("data/raw/Proteins Report of merged human files with all 4 conditions (1. Feb, 5. Feb, 18. Feb).csv", skip = 1, header = TRUE, check.names = FALSE) #it is just the E1 SYPRO 5x from the original raw data from 220420
#Import csv file for renaming of Bait and AP-names 
data_rename<- read.csv("data/metadata/SAINT_rename_file.csv", header = TRUE, check.names = FALSE) #it is just the E1 SYPRO 5x from the original raw data from 220420

#Modify file ----

# Rename the columns
data_SAINT_input <- data_raw_scaffold %>%
  rename(
    `#BaitName_Condition` = `Biological sample category`,
    `#APName` = `Biological sample name`,
    `#PreyName` = `Protein accession number`,
    `#SpectralCounts` = `Number of total spectra`
  )

# Select only the renamed columns , select suddenly did not work therefore using dplyr::select
data_SAINT_input <- data_SAINT_input %>% dplyr::select(`#BaitName_Condition`, `#APName`, `#PreyName`, `#SpectralCounts`)

# Remove duplicate rows
data_SAINT_input <- data_SAINT_input %>% distinct()

#  Remove all rows where column "Number of total spectra" == 0
data_SAINT_input <- data_SAINT_input %>% filter(`#SpectralCounts` != 0)

# Remove all rows with random AP names, e.g. the column "Protein accession number" contains "|"
data_SAINT_input <- data_SAINT_input %>% filter(!grepl("\\|", `#PreyName`))

# Rename entries in #BaitName_Condition and #APName FindReplace
#requirements for SAINT input file in list format following https://reprint-apms.org/?q=inputfileformatting
#Note, the CONTROL should not contain any "_" , the Bait names can not contain any "-" 
data_SAINT_input<-FindReplace(data_SAINT_input, "#BaitName_Condition", data_rename, from = "#BaitName_Condition_old", to = "#BaitName_Condition_new", exact = TRUE, vector = FALSE)
data_SAINT_input<-FindReplace(data_SAINT_input, "#APName", data_rename, from = "#APName_old", to = "#APName_new", exact = TRUE, vector = FALSE)

# Save data ----
# Save the cleaned data to a new CSV file
write.csv(data_SAINT_input, "data/interim/SAINT_list_input.csv", row.names = FALSE, quote = FALSE)




#Second SAINT input file ----

# Modified input file where Empty vector is removed as a control
data_SAINT_input_control <- data_SAINT_input %>%
  filter(`#BaitName_Condition` != "CONTROL")

# Re-name MTHFR to CONTROL
data_SAINT_input_control$`#BaitName_Condition` <- replace(data_SAINT_input_control$`#BaitName_Condition`, 
                                                          data_SAINT_input_control$`#BaitName_Condition` == "MTHFR", "CONTROL")

# Save the cleaned data to a new CSV file
write.csv(data_SAINT_input_control, "data/interim/SAINT_list_input_MTHFR_control.csv", row.names = FALSE, quote = FALSE)


#Third SAINT input file with modified acession numbers ----
#Retrieve gene names using mapUniProt
data_accession_names<-mapUniProt( from="UniProtKB_AC-ID", to='RefSeq_Protein', query=unique(data_SAINT_input$`#PreyName`) )
#allToKeys() #check which types of accession names we can retrieve 

#Add column gene_name for converted Uniprot accession numbers to Gene names
#Note seems like some aceession numers from Uniprot can not be found in RefSeq
data_SAINT_input_RefSeq<-data_SAINT_input %>%
  mutate('#PreyName'=
           FindReplace(data_SAINT_input, '#PreyName', data_accession_names, from = "From", to = "To", exact = TRUE, vector = TRUE))

# Save the cleaned data to a new CSV file
#Note Even with this modification, the SAINT analysis fail with crapome 
write.csv(data_SAINT_input_RefSeq, "data/interim/SAINT_list_input_MTHFR_RefSeq.csv", row.names = FALSE, quote = FALSE)



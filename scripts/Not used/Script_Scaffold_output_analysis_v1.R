library(ggplot2)
library(dplyr)
library(readr)
library(stringr)
library(UniProt.ws)
library(gridExtra)
library("ggrepel")
library(DataCombine)
library(UpSetR)


#Read the input file 
# Data files ----
getwd()

#Import raw csv file. Manually converted from output Scaffold .txt file prior. 
data_raw_scaffold<- read.csv("data/raw/Proteins Report of merged human files with all 4 conditions (1. Feb, 5. Feb, 18. Feb).csv", skip = 1, header = TRUE, check.names = FALSE) #it is just the E1 SYPRO 5x from the original raw data from 220420


# Remove rows with NA in the specified column using dplyr
data_scaffold<- na.omit(data_raw_scaffold)


# Renaming specific values within a column using dplyr
data_scaffold <- data_scaffold %>%
  mutate(`Biological sample category` = recode(`Biological sample category`,
                                               "control" = "CONTROL",
                                               "group 1" = "MTHFR_WT",
                                               "group 2" = "MTHFR_T34A",
                                               "group 3" = "MTHFR_trunc"))

# Use mutate and case_when to create the new column called replicates 
data_scaffold <- data_scaffold %>%
  mutate(replicates = case_when(
    str_detect(`Biological sample name`, "1\\.Feb|1\\. Feb") ~ "1 Feb" ,
    str_detect(`Biological sample name`, "18\\.Feb|18\\. Feb") ~ "18 Feb",
    str_detect(`Biological sample name`, "5\\.Feb|5\\. Feb") ~ "5 Feb"
    ))

#Retrieve gene names using mapUniProt
#data_accession_names<-mapUniProt( from="UniProtKB_AC-ID", to='RefSeq_Protein', query=unique(data_raw_SAINT$Prey) )
data_accession_names<-mapUniProt( from="UniProtKB_AC-ID", to='Gene_Name', query=unique(data_scaffold$`Protein accession number`) )
#allToKeys() #check which types of accession names we can retrieve 

#Add column gene_name for converted Uniprot accession numbers to Gene names 
data_scaffold<-data_scaffold %>%
  mutate(gene_name=
           FindReplace(data_scaffold, 'Protein accession number', data_accession_names, from = "From", to = "To", exact = TRUE, vector = TRUE))




plot_scaffold <-ggplot(data_scaffold,
        aes(x = replicates,
            y = `Number of total spectra`)
 )+

  facet_grid(.~ `Biological sample category`) + 
  geom_point(aes(group=replicates),
             position = position_dodge(width = 0.2)) +  # Using a bar plot; adjust `geom` as needed for your data
  
  
  #highlight sample names
  geom_point(
    data = filter(
      data_scaffold,
      gene_name == "MTHFR" |
        gene_name == "MTHFD1" |
        gene_name == "GCN1" |
        gene_name == "PRKDC" |
        gene_name == "FASN"
    ),
    aes(color = gene_name)
  )+
  

  
  labs(title = "Scaffold data",
        x = "Biological Sample Category",
        y = "Number of Total Spectra") +
   theme_minimal()+
  theme(
    axis.text.x = element_text(angle = -45, hjust = 0)
  )
plot_scaffold
 ggsave("Plot_scaffold_output.png",plot=plot_scaffold,  height = 2000, width=2000, units = "px", path =getwd() )
 
 
 
 #Filter data 
 
 data_scaffold_filter<- data_scaffold %>%
  filter(`Protein accession number` == "P78527") 
   
 plot_scaffold_filter<- ggplot(data= data_scaffold_filter, aes(x = replicates, y = `Number of total spectra`)) +
   facet_grid(. ~ `Biological sample category`) +
   geom_point(aes(group = replicates), position = position_dodge(width = 0.5)) +  # Using a bar plot; adjust `geom` as needed for your data
   labs(title = "Scaffold data", x = "Biological Sample Category", y = "Number of Total Spectra") +
   scale_y_log10() +
   theme_minimal()
 plot_scaffold_filter
   
   ggsave("Plot_scaffold_output_filter.png",plot=plot_scaffold_filter,  height = 2000, width=2000, units = "px", path =getwd() )
   
   
   #Ven----
   upsetSamples(miniACC)
   ??upsetSamples
   
   

library(STRINGdb)
library(dplyr)

string_db <- STRINGdb$new( version="12.0", species=9606,
    score_threshold=100, network_type="full", input_directory="")

dat_example = data.frame(list(genes=c('MTHFR', 'MMADHC', 'MMUT', 'PRKDC')))

dat_examples_mapped = string_db$map(dat_example, 'genes')


string_db$get_interactions(dat_examples_mapped$STRING_id) %>% 
    merge(dat_examples_mapped %>%
        rename(from_gene=genes),
        by.x='from', by.y='STRING_id') %>%
    merge(dat_examples_mapped %>%
        rename(to_gene=genes),
        by.x='to', by.y='STRING_id')


# string_db$get_summary(dat_examples_mapped$STRING_id, required_score=200)



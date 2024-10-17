if (exists("snakemake")) {
  string_dir <- snakemake@output$string_db
}

library(STRINGdb)

dir.create(string_dir)
string_db <- STRINGdb$new(
  version = "12.0",
  species = 9606,
  score_threshold = 100,
  network_type = "full",
  input_directory = string_dir
)

res <- string_db$map(data.frame(list(gene = c("MTHFR", "MMACHC"))), "gene")
string_db$get_interactions(res$STRING_id)

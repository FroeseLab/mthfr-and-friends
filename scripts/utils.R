suppressPackageStartupMessages({
    library(UniProt.ws)
})


#' Translate uniprot ids to gene names
#' 
#' @param query_ids vector of uniprot ids
#' @return vector of gene names
translate_up <- function(query_ids){
    # translate uniprot to gene names
    out = mapUniProt(
        from = "UniProtKB_AC-ID",
        to = "Gene_Name",

        query=as.vector(query_ids)
    )

    genes = out[['To']]
    names(genes) = out[['From']]

    return(genes[query_ids])

}

#' Add gene names to a qlf table
#' 
#' @param qlf qlf object
#' @return qlf object with gene names on table
add_genes <- function(qlf){
    # Add gene names to the table
    qlf$table$gene = translate_up(row.names(qlf$table))
    return(qlf)
}

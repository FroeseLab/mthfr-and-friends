library(httr)
library(readr)


endpoint <- "https://sparql.uniprot.org/"

query <- "PREFIX skos: <http://www.w3.org/2004/02/skos/core#>
PREFIX up: <http://purl.uniprot.org/core/>
PREFIX taxon: <http://purl.uniprot.org/taxonomy/>
PREFIX GO: <http://purl.obolibrary.org/obo/GO_>

SELECT DISTINCT
    ?gene_name (?accession as ?uniprot_accession)
WHERE
{   
    ?protein a up:Protein ;
        up:reviewed true  ;
        # Here we query for 'human' proteins
        up:organism taxon:9606 ;
        # That are classified to have the GO term 'Kinase activity'
        up:classifiedWith|(up:classifiedWith/rdfs:subClassOf) GO:0016301 ;
        up:encodedBy ?gene .
        ?gene skos:prefLabel ?gene_name .
         BIND(SUBSTR(STR(?protein),33) AS ?accession)
}
"

# query the endpoint using a httr get request with the query = query
response <- httr::GET(url = endpoint, query = list(query = query, format = "csv"))
# parse the response content into a dataframe
data <- httr::content(response, "parsed")

# write the dataframe to a csv file
write.csv(data, file = "data/interim/uniprot_kinases.csv", row.names = FALSE)

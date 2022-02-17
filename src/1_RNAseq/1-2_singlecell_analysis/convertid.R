
# functions to convert genes ----------------------------------------------


# converts chicken or anything from ensembl to gene id 

convertid <- function(genelist, organism = "ggallus", supplied_id = "ensembl_gene_id", attribs = c('ensembl_gene_id', 'external_gene_name', 'description', 'mmusculus_homolog_ensembl_gene'), host = "asia.ensembl.org", unique_rows = FALSE){

  require(biomaRt)
  
  speciesmart <- useMart("ensembl", dataset = paste0(organism ,"_gene_ensembl"), host = host)
  
  print("some of the attributes")
  
  listAttributes(speciesmart) %>% head(30) %>% print()
  
  converted_genes <- getBM(attributes=attribs, 
        filters = supplied_id, 
        values = genelist, 
        mart = speciesmart,
        uniqueRows = unique_rows
  )
  
  # fill blanks of external gene name with ensembl
  converted_genes$external_gene_name[converted_genes$external_gene_name == ""] <- converted_genes$ensembl_gene_id[converted_genes$external_gene_name == ""]
  
  
  return(converted_genes)
  
}


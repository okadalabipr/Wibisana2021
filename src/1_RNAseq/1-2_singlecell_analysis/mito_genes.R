# function to get mitochondrial related genes


#connect to annotation hub 
ah <- AnnotationHub()

ahDb <- query(ah, 
              pattern = c("Gallus gallus", "EnsDb"), 
              ignore.case = TRUE)
#get latest annotation files
id <- ahDb %>%
  mcols() %>%
  rownames() %>%
  tail(n = 1)

# Download the appropriate Ensembldb database
edb <- ah[[id]]

# Extract gene-level information from database
annotations <- genes(edb, 
                     return.type = "data.frame")  


# Select annotations of interest
annotations <- annotations %>%
  dplyr::select(gene_id, gene_name, gene_biotype, seq_name, description, entrezid)

# Extract IDs for mitochondrial genes
mt <- annotations %>%
  dplyr::filter(seq_name == "MT") %>%
  dplyr::pull(gene_id)

write.table(mt, "mitochondrial_genes.txt")


library(plotly)
library(cowplot)
library(tidyverse)
library(ggpubr)
library(clusterProfiler)
library(biomaRt)

SE <- read_csv("/SE.csv")
TE <- read_csv("/TE.csv")

marker_genes <- read_csv("activated_markers_cv.csv") %>% dplyr::select(-1)

marker_genes <- marker_genes %>% 
  filter(p_val_adj < 0.05) %>%
  mutate(state = case_when(avg_log2FC > avg_log2FC %>% quantile(0.75) ~ "high",
                           avg_log2FC < avg_log2FC %>% quantile(0.25) ~ "low",
                           avg_log2FC <= avg_log2FC %>% quantile(0.75) & avg_log2FC >= avg_log2FC %>% quantile(0.25) ~ "unchanged"
  ))

intersection_list <- list(
  upregulated = SE_marker_genes_intersection %>% filter(state == "high", atac_state == "atac_gain") ,
  downregulated = SE_marker_genes_intersection %>% filter(state == "low", atac_state == "atac_lost"),  
  unchanged = SE_marker_genes_intersection %>% filter(state == "unchanged", atac_state == "atac_unchanged")   
)

# convert id to mouse gene id
convert_join <- function(x){
  converted_id <- convertid(x$ensembl_gene_id, host = "useast.ensembl.org")
  converted_df <- inner_join(x, converted_id, by = "ensembl_gene_id") %>% 
    distinct(ensembl_gene_id, .keep_all = T) 
  
  return(converted_df) 
}

intersection_list <- intersection_list %>% lapply(convert_join)

# write GO plot
go_strings <- c("_go_simplified")
go_strings_2 <- c("BP")

# function for plot
barplot_cp <- function(x, font_size = 20, n = 20){

  df_clust <- x
    
  total_genes <- df_clust$GeneRatio[1] %>% str_replace_all(".*/", "") %>% as.numeric()  
  df_clust$gene_fraction <- df_clust$Count/total_genes
  
  # calculate log qvalue
  df_clust$log.p.adjust <- -log10(df_clust$p.adjust)
  
  # plot barplot
  ggplot(df_clust %>% head(n = n), aes(x = reorder(Description, log.p.adjust), y = gene_fraction, fill = log.p.adjust)) + 
    geom_bar(stat = "identity") +
    coord_flip() + theme_cowplot(font_size) + cowplot::panel_border("black", size =1) +
    labs(fill = "-log(p-adjust)") +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank()) +
    scale_fill_viridis_c()  
  
}


write_go <- function(x, resdir = "", width = 15, height = 7){
    gene_no <- x %>% pull("mmusculus_homolog_ensembl_gene") %>% str_subset(pattern = "") %>% length()
    # ONT = BP
    go_file <- enrichGO(x %>% pull("mmusculus_homolog_ensembl_gene") %>% str_subset(pattern = ""), 
                        OrgDb = "org.Mm.eg.db",
                        ont = go_strings_2[k],
                        keyType = "ENSEMBL") 
    
    go_plot <- go_file %>% simplify() %>% as_tibble() %>% barplot_cp()
    
    # mkdir
    dir.create(resdir)
    
    ggsave(paste0(resdir, names(x)[i], "_go_", go_strings_2, "_", gene_no, ".svg"), go_plot, width = width, height = height)
    
    write_csv(go_file %>% as_tibble(), paste0(resdir, names(x), "_go_", go_strings_2,"_", gene_no, ".csv"))
    write_csv(go_file %>% simplify() %>% as_tibble(), paste0(resdir, names(x), go_strings,"_", gene_no, ".csv"))
    
}

write_go(intersection_list[[1]])
library(tidyverse)

# read SE and TE related genes
SE <- read_csv("SE.csv")
TE <- read_csv("TE.csv")

marker_genes <- read_csv("activated_markers.csv") %>% dplyr::select(-1)

marker_genes <- marker_genes %>% 
  filter(p_val_adj < 0.05) %>%
  mutate(state = case_when(avg_log2FC > avg_log2FC %>% quantile(0.75) ~ "high",
                           avg_log2FC < avg_log2FC %>% quantile(0.25) ~ "low",
                           avg_log2FC <= avg_log2FC %>% quantile(0.75) & avg_log2FC >= avg_log2FC %>% quantile(0.25) ~ "unchanged"
  ))

# read annotated connections
conns_annotated <- read_tsv("conns_annotated.txt")

conns <- read_tsv("conns.tsv", col_types = cols(chr = col_character()))

conns_annotated <- conns_annotated %>% select(peakid = 1, chr = "Chr", start = "Start", end = "End", ensembl_gene_id = `Entrez ID`, external_gene_name = `Gene Name`)

conns <- inner_join(conns, conns_annotated %>% select(ensembl_gene_id, external_gene_name, peakid))

# calculate number of connections per annotated gene
conns <- left_join(conns, conns %>%
 count(ensembl_gene_id)) %>%
  left_join(conns %>% count(ensembl_gene_id, dist_type) %>%
  pivot_wider(names_from = dist_type, values_from = n))

conns_distinct <- conns %>% distinct(across(ensembl_gene_id), .keep_all = T) %>% select(ensembl_gene_id, n, long, medium, short)

# plot connections
# make function for plotting
corr_plot_conns <- function(df_1, atac_st = "atac_gain", x_axis = "n", dir_out){
   plot <- left_join(df_1 %>% filter(atac_state == atac_st), fano_df %>% select(ensembl_gene_id, fano_10, fano_0) %>% mutate(l2fano = log2(fano_10/fano_0)) %>% na.omit()) %>%
    left_join(conns_distinct) %>%
    filter(ensembl_gene_id %in% (marker_genes %>% pull(ensembl_gene_id))) %>%
    ggplot(aes(x = log2(!!rlang::sym(x_axis)), y = l2fano, color = as.numeric(log2(!!rlang::sym(x_axis))))) +
    geom_point(size = 3, alpha = 1)  + 
    stat_smooth(method = "lm", col = "black") +
    stat_cor(method = "spearman", size = 10) +
    scale_color_gradientn(colours=colorRampPalette(c("#386cb0", "#ef3b2c"))(100))+
    geom_hline(yintercept = 0)+
    guides(color = F)+
    cowplot::theme_cowplot(font_size = 40, line_size = 1)
  
   
  # print_N 
  tot_pts <- left_join(df_1 %>% filter(atac_state == atac_st), fano_df %>% select(ensembl_gene_id, fano_10, fano_0) %>% mutate(l2fano = log2(fano_10/fano_0)) %>% na.omit()) %>%
     left_join(conns_distinct) %>%
     filter(ensembl_gene_id %in% (marker_genes %>% pull(ensembl_gene_id))) %>% pull(!!rlang::sym(x_axis)) %>% na.omit() %>% length() %>% as.character()
  
   dir.create(dir_out)
   ggsave(paste0(dir_out, deparse(substitute(df_1)), "_", x_axis, "_", atac_st,  tot_pts, ".svg"), plot, width = 7, height = 7)
   
}


# calculate fano factor and mRNA change correlation --------------------------



fano_df <- read_csv("baynorm_fano.csv")

fano_df <- fano_df %>% mutate(label = case_when(external_gene_name == "CD83" ~ 1, 
                                     external_gene_name == "NFKBIA" ~ 1,
                                     TRUE ~ 0))

# SE and marker genes
inner_join(fano_df, SE, by = "ensembl_gene_id") %>%
  filter(ensembl_gene_id %in% (marker_genes %>% pull(ensembl_gene_id))) %>%
  ggplot(aes(x = l2fc, y = log2(fano_10/fano_0), color = N))+
  geom_point(size = 3, alpha = 1)  + 
  stat_smooth(method = "lm", col = "black") +
  stat_cor(method = "spearman", size = 10) +
  scale_color_gradientn(colours=colorRampPalette(c("#ef3b2c", "#386cb0"))(100))+
  geom_hline(yintercept = 0)+
  guides(color = F)+
  geom_text_repel(inner_join(fano_df, SE, by = "ensembl_gene_id") %>%
                    filter(ensembl_gene_id %in% (marker_genes %>% pull(ensembl_gene_id))) %>% filter(label !=0), min.segment.length = 0, mapping = aes(x=l2fc, y=log2(fano_10/fano_0), label = external_gene_name.x), size = 8.0, nudge_y = 3, nudge_x = 0, segment.alpha = 0.7, color = "black")+
  cowplot::theme_cowplot(font_size = 40, line_size = 1)

ggsave("SE_marker.svg", width = 7, height = 7)

# TE and marker genes
inner_join(fano_df, TE, by = "ensembl_gene_id") %>%
  filter(ensembl_gene_id %in% (marker_genes %>% pull(ensembl_gene_id))) %>%
  ggplot(aes(x = l2fc, y = log2(fano_10/fano_0), color = N))+
  geom_point(size = 3, alpha = 1)  + 
  stat_smooth(method = "lm", col = "black") +
  stat_cor(method = "spearman", size = 10) +
  scale_color_gradientn(colours=colorRampPalette(c("#ef3b2c", "#386cb0"))(100))+
  geom_hline(yintercept = 0)+
  guides(color = F)+
  cowplot::theme_cowplot(font_size = 40, line_size = 1)


ggsave("TE_marker.svg", width = 7, height = 7)


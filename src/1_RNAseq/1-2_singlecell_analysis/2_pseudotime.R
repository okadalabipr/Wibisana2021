# pseudotime calculation

# read coordinates of cells
umap_cell_embeddings <- read_csv("umap_cell_embeddings.csv")

# add information of dose and clusters to the cell coordinates
cell_clusters <- clean_seurat@meta.data %>% 
    dplyr::select(seurat_clusters, dose = sample) %>% 
    rownames_to_column("cell") %>% 
    as_tibble()


# combine cell clusters and umap coordinates
umap_coord <- inner_join(umap_coord, cell_clusters)

# calculate principal curve for the points in umap space
fit_curve <- princurve::principal_curve(umap_coord[,2:3] %>% as.matrix(), smoother = "lowess")
fit_curve_df <- fit_curve$s %>% as_tibble() %>% dplyr::rename(fitted_x = UMAP_1, fitted_y = UMAP_2)

umap_coord_fitted <- bind_cols(umap_coord, fit_curve_df)

# calculate distance between consecutive points in pseudotime axis
# function for euclidean distance
euclid_dist <- function(x1, x2, y1, y2){sqrt((x1 - x2)^2 + (y1 - y2)^2)}

# arrange cells by fitted_x (x coordinate)
umap_coord_fitted <- arrange(umap_coord_fitted, fitted_y)

# calculate euclidean distance to project 2d line to 1d line based on distance, add 0 as first element
umap_coord_fitted$fitted <- sapply(1:(nrow(umap_coord_fitted)-1), function(i){euclid_dist(umap_coord_fitted$fitted_x[i], umap_coord_fitted$fitted_x[i+1], umap_coord_fitted$fitted_y[i], umap_coord_fitted$fitted_y[i+1])}
) %>% prepend(0) %>% cumsum()

# flip pseudotime
umap_coord_fitted$fitted <- rev(umap_coord_fitted$fitted)

# write pseudotime plot
pseudotime_plot <- ggplot(umap_coord_fitted) + 
  geom_point(aes(x = UMAP_1, y = UMAP_2, fill = seurat_clusters,), shape = 21, size = 20, stroke = NA) + 
  geom_path(aes(x = fitted_x, y = fitted_y, color = fitted), size =4) + 
  scale_fill_manual(values = c("#386cb0","#ef3b2c")) + 
  scale_color_viridis_c(option = "magma") +
  cowplot::theme_cowplot(20, line_size = 1) +
  theme(legend.position = "none") +
  theme(axis.line.x = element_line(size = 2),
        axis.line.y = element_line(size = 2))+
  labs(color = "Pseudotime (a.u)", fill = "Clusters")

pseudotime_plot_dose <- ggplot(umap_coord_fitted) + 
    geom_point(aes(x = UMAP_1, y = UMAP_2, fill = as.factor(dose)), shape = 21, size = 20, stroke = NA) + 
    geom_path(aes(x = fitted_x, y = fitted_y, color = fitted), size =4) + 
    scale_fill_viridis_d() + 
    scale_color_viridis_c(option = "magma") +
    cowplot::theme_cowplot(20, line_size = 1) +
    theme(legend.position = "none") +
    theme(axis.line.x = element_line(size = 2),
        axis.line.y = element_line(size = 2))+
    labs(color = "Pseudotime (a.u)", fill = "Clusters")



# get genes expression over pseudotime ------------------------------------
seurat_normalized_df <- clean_seurat@assays$RNA@data %>% 
    as.data.frame() %>%
    rownames_to_column("ensembl_gene_id") %>%
    as_tibble()

seurat_normalized_df <- inner_join(converted_id[,1:2], seurat_normalized_df) %>%
    as_tibble()

# filter counts by genes acquired using findallmarkers
seurat_normalized_df <- activated_markers %>% inner_join(seurat_normalized_df)
seurat_normalized_df <- activated_markers %>% dplyr::select(10, 5, 8) %>% inner_join(seurat_normalized_df)


# Plot gene expression againts pseudotime
for(i in 1:length(seurat_normalized_df$external_gene_name)){
  
  genename <- seurat_normalized_df$external_gene_name[i]
  
  l2fc <- seurat_normalized_df$avg_log2FC[i] %>% signif(3)
  
  padj <- seurat_normalized_df$p_val_adj[i] %>% signif(3)
  
  pseudotime_gene <- seurat_normalized_df %>% 
    filter(external_gene_name == genename) %>%
    dplyr::select(-1:-4) %>%
    reshape2::melt() %>% 
    rename(cell = "variable", expression = "value") %>%
    inner_join(umap_coord_fitted %>% dplyr::select(cell, fitted, dose, seurat_clusters), by = "cell") %>%
    ggplot() + geom_point(aes(x = fitted, y = expression, color = as.factor(dose)), size = 4) + geom_smooth(aes(x = fitted, y = expression), method = "loess")+ 
    cowplot::theme_cowplot(15) + xlab("Pseudotime") + ylab("Normalized count") + ggtitle(genename, subtitle = paste0("log2fc = ", l2fc, ", p_adj = ", padj)) + labs(color = "M4 (ug/ml)")+
    scale_color_viridis_d()

  ggsave(paste0(genename, ".png"), width = 10, height = 7)
  
}
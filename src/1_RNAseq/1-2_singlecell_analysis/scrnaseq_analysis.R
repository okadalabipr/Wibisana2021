library(bayNorm)
library(Seurat)
library(tidyverse)
library(magrittr)
library(RCurl)
library(gt)

source("../function_compilations/biomart_convert_geneid.R")


# specify directories
counts_dir <- "directory of counted data/counts.txt"


# data preprocessing ------------------------------------------------------
# read counts
counts <- read_csv(counts_dir, skip = 2) 

# extract columns containing counts only
counts <- counts[,-2:-6]

# gene id conversion
converted_id <- convertid(counts$Geneid, attribs = c('ensembl_gene_id', 'external_gene_name', 'description'), host = "uswest.ensembl.org")

# make seurat object ------------------------------------------------------
RNA_seurat <- CreateSeuratObject(counts %>% as.data.frame() %>% column_to_rownames("Geneid"))

# add dose metadata depending on cell naming
RNA_seurat@meta.data <- RNA_seurat@meta.data %>%  
  mutate(sample = case_when(
    orig.ident == "B00" ~ 0,
    orig.ident == "B01" ~ 0.01,
    orig.ident == "B02" ~ 0.1,
    orig.ident == "B03" ~ 1,
    orig.ident == "B04" ~ 10
  ))

# add rownames
rownames(RNA_seurat@meta.data) <- colnames(counts)[-1]

# data qc -----------------------------------------------------------------
# add number of genes per non-normalized values for each cell for better comparison between samples
RNA_seurat$log10GenesPerUMI <- log10(RNA_seurat$nFeature_RNA) / log10(RNA_seurat$nCount_RNA)

metadata <- RNA_seurat@meta.data

# Rename columns
metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)


# compute mitochondrial genes ---------------------------------------------

# get mito genes
source("mito_genes.R")

mt <- read.table("../20191106_scrnaseqclustering/mitochondrial_genes.txt", stringsAsFactors = F)$x

# count UMIs assigned to mitochondrial genes
mtUMI <- rbind(
  Matrix::colSums(counts[,-1][which(counts$Geneid %in% mt),], na.rm = T) %>% as.data.frame()
)

metadata$mtUMI <- mtUMI$.

# Calculate mitochondrial gene ratio per cell
metadata$mitoRatio <- metadata$mtUMI/metadata$nUMI

# return metadata
RNA_seurat@meta.data <- metadata

#  filter cells -----------------------------------------------------------

# Filter out low quality reads using selected thresholds - these will change with experiment
filtered_seurat <- subset(RNA_seurat, 
                          subset = nUMI >= 1500000 & 
                            nGene >= 8500 & 
                            #genes have low novelty
                            #(log10GenesPerUMI > 0) & 
                            mitoRatio < 0.04)


metadata_filtered <- filtered_seurat@meta.data

# clean genes -------------------------------------------------------------

# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts_logical <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts_logical > 0L

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- rowSums(as.matrix(nonzero)) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts_logical[keep_genes, ]

# Create a new Seurat object
clean_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

metadata_clean <- clean_seurat@meta.data

#save filtered ensembl id
filtered_ensembl_id <- clean_seurat@assays$RNA@counts %>% rownames()

# function for clustering -------------------------------------------------

# backup seurat data
clean_seurat <- filtered_seurat

# cell cycle scoring ------------------------------------------------------

# Normalize the counts
seurat_phase <- NormalizeData(clean_seurat)

# Cell cycle scoring using curated list (Tirosh et al., 2016)
cc_file <- getURL("https://raw.githubusercontent.com/hbc/tinyatlas/master/cell_cycle/Homo_sapiens.csv") 
cell_cycle_genes <- read.csv(text = cc_file)

# Basic function to convert mouse to human gene names
library(biomaRt)

chicken_ccgenes <- convertid(cell_cycle_genes$geneID, organism = "hsapiens", attribs = c("ensembl_gene_id", "ggallus_homolog_ensembl_gene"), host = "asia.ensembl.org")

cell_cycle_genes <- cell_cycle_genes %>% rename(ensembl_gene_id = "geneID")

cell_cycle_genes <- inner_join(cell_cycle_genes, chicken_ccgenes, by = "ensembl_gene_id")

cell_cycle_genes <- cell_cycle_genes[!cell_cycle_genes$ggallus_homolog_ensembl_gene == "",]

# Acquire the S phase genes
s_genes <- cell_cycle_genes %>%
  dplyr::filter(phase == "S") %>%
  pull("ggallus_homolog_ensembl_gene")

# Acquire the G2M phase genes        
g2m_genes <- cell_cycle_genes %>%
  dplyr::filter(phase == "G2/M") %>%
  pull("ggallus_homolog_ensembl_gene")


# Perform cell cycle scoring
seurat_phase <- CellCycleScoring(seurat_phase,
                                 g2m.features = g2m_genes,
                                 s.features = s_genes)

# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                                     selection.method = "vst",
                                     nfeatures = 2000, 
                                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
cellcycleplot <- DimPlot(seurat_phase,
                         reduction = "pca",
                         group.by= "Phase",
                         split.by = "Phase")

cellcycleplot_dose <- DimPlot(seurat_phase,
                              reduction = "pca",
                              group.by= "sample",
                              split.by = "Phase")

# cell cycle scoring ------------------------------------------------------

# perform cell cycle scoring on all samples
clean_seurat <- NormalizeData(clean_seurat, verbose = TRUE)
clean_seurat <- FindVariableFeatures(clean_seurat, selection.method = "vst", nfeatures = 2000)
clean_seurat <- ScaleData(clean_seurat)
clean_seurat <- CellCycleScoring(clean_seurat, g2m.features=g2m_genes, s.features=s_genes)

# scale data based on mitoratio and cell cycle scoring
clean_seurat <- ScaleData(clean_seurat, vars.to.regress = c("mitoRatio", "S.Score", "G2M.Score"))
clean_seurat <- RunPCA(clean_seurat)

DimPlot(clean_seurat, group.by = "Phase")

# perform data normalization using baynorm ---------------------------------------------------------
# get the processed cells filtered using seurat
filtered_cells <- rownames(clean_seurat@meta.data)
counts_filtered <- counts %>% dplyr::select(1, all_of(filtered_cells))

# filter features
counts_filtered <- counts_filtered %>% filter(Geneid %in% filtered_ensembl_id)

bay_out <- bayNorm(counts_filtered %>% as.data.frame() %>% column_to_rownames("Geneid"), mean_version = TRUE, UMI_sffl = 10)

bay_out_counts <- bay_out$Bay_out %>% as.data.frame() %>% rownames_to_column("Geneid") %>% as_tibble()

# create baynorm seurat object 
baynorm_seurat <- CreateSeuratObject(counts = bay_out$Bay_out, assay = "bayNorm")

# combine assays ----------------------------------------------------------
# store normalized seurat data in misc
clean_seurat@misc[["seurat_data"]] <- as.matrix(x = clean_seurat@misc)

# change RNA matrix to baynorm normalized data
clean_seurat@assays$RNA@data <- baynorm_seurat@assays$bayNorm@data
clean_seurat@assays$RNA@counts <- baynorm_seurat@assays$bayNorm@data

# data normalization ------------------------------------------------------
clean_seurat <- NormalizeData(clean_seurat, verbose = TRUE)
clean_seurat <- FindVariableFeatures(clean_seurat, selection.method = "vst", nfeatures = 2000)

# scale data based on mitoratio and cell cycle scoring
clean_seurat <- ScaleData(clean_seurat, vars.to.regress = c("mitoRatio", "S.Score", "G2M.Score"))

# clustering cells based on top PCs ---------------------------------------
# run PCA
clean_seurat %<>% RunPCA()

# Explore heatmap of PCs
PC_heatmap <- DimHeatmap(clean_seurat, 
                         dims = 1:9, 
                         cells = 500, 
                         balanced = TRUE) 

PC_heatmap

# Printing out the most variable genes driving PCs
print(x = clean_seurat[["pca"]], 
      dims = 1:10, 
      nfeatures = 10) 

# plot feature selection
top20_genes <- convertid(head(VariableFeatures(clean_seurat), 20), attribs = c("ensembl_gene_id", "external_gene_name")) %>% pull(external_gene_name)

var_features_plot <- VariableFeaturePlot(clean_seurat)
var_features_plot_label <- LabelPoints(var_features_plot, points = head(VariableFeatures(clean_seurat), 20), labels = top20_genes, repel = TRUE)

var_features_plot_label

# Determine percent of variation associated with each PC
pct <- clean_seurat[["pca"]]@stdev / sum(clean_seurat[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.
print(paste("last 0.1% change point" , co2))

# Minimum of the two calculation
pcs <- min(co1, co2)

#cutoffpoint[[i]] <- pcs

# Create a dataframe with values
plot_df <- data.frame(pct = pct, 
                      cumu = cumu, 
                      rank = 1:length(pct))

# Elbow plot to visualize 
elbowplot  <- ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  #    geom_vline(xintercept = 90, color = "grey") + 
  #    geom_hline(yintercept = min(pct[pct > pcs]), color = "grey") +
  cowplot::theme_cowplot()

elbowplot

# cell clustering ---------------------------------------------------------

# Determine the K-nearest neighbor graph
clean_seurat <- FindNeighbors(object = clean_seurat, 
                              dims = 1:pcs)

# Determine the clusters for various resolutions                                
clean_seurat <- FindClusters(object = clean_seurat,
                             resolution = c(0.2))


# Calculation of dimensionality reduction
clean_seurat <- RunUMAP(clean_seurat, 
                        reduction = "pca", 
                        dims = 1:pcs)



# plot
umap_clust <- DimPlot(clean_seurat,
                      reduction = "umap",
                      label = T,
                      label.size = 0,
                      pt.size = 6,
                      #group.by = "sample",
                      cols = c("#386cb0", "#ef3b2c"),) + theme(axis.line = element_line(size = 1))
umap_clust

umap_dose <- DimPlot(clean_seurat,
                     reduction = "umap",
                     label = F,
#                     label.size = 0,
                     pt.size = 6,
                     group.by = "sample",
                     cols = c(viridis::viridis(5))) + theme(axis.line = element_line(size = 1))

umap_dose 

umap_phase <- DimPlot(clean_seurat,
                      reduction = "umap",
                      label = T,
                      label.size = 6,
                      pt.size = 4,
                      group.by = "Phase") + theme(axis.line = element_line(size = 1))


umap_phase 

# save umap coords
write.csv(clean_seurat@reductions$umap@cell.embeddings, "umap_cell_embeddings.csv")

# save csv for cell metadata
write.csv(clean_seurat@meta.data, "metadata_final.csv")

# find DE genes -----------------------------------------------------------
all_markers <- FindAllMarkers(clean_seurat, logfc.threshold = 0)

marker_clust_1 <- all_markers %>% rename(ensembl_gene_id = "gene") %>% inner_join(converted_id, by = "ensembl_gene_id") %>% filter(cluster == 1) %>% arrange(-avg_logFC)

# convert lnfc to log2fc
marker_clust_1$avg_log2FC <- exp(marker_clust_1$avg_logFC) %>% log2()

write_csv(marker_clust_1, "activated_markers.csv")

# save seurat object as rds object ----------------------------------------

saveRDS(clean_seurat, "seurat_object.RDS")

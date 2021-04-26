
library(monocle3)
library(cicero)
library(tidyverse)
library(rtracklayer)


# specify directories

# peak matrix directory
dir_peakmatrix <- "peak matrix RDS directory/peak_matrix.RDS"
dir_10x <- "cellranger outs folder directory/singelecell.csv"
# refer to chrom_size.sh to calculate chromosome size
dir_genomesizes <- "directory of /genomesize.tsv"
dir_bed <- "directory of bed file /enhancers.bed"
dir_gtf <- "directory of gtf file/genome.gtf"

# set visualization regions c("chromosome num", start, end)
region_to_vis <- c(2, 60658235, 60710751)

# loading 10x data --------------------------------------------------------
# load cellranger filtered cells
cell_barcode <- read_csv(dir_10x) %>% 
    filter(is__cell_barcode == 1) %>%
    pull(barcode)

# read peak matrix file
indata <- readRDS(dir_peakmatrix)

# filter cells
indata <- indata$all_peaks_unstitched[,colnames(indata$all_peaks_unstitched) %in% cell_barcode]

# binarize matrix
indata@x[indata@x > 0] <- 1

# make CDS
input_cds <-  suppressWarnings(new_cell_data_set(indata,
                                                 cell_metadata = NULL,
                                                 gene_metadata = NULL))

input_cds <- detect_genes(input_cds)

# Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 

# constructing cis regulatory networks ------------------------------------

# perform dimensionality reduction
set.seed(2017)
# input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', 
                              preprocess_method = "LSI")

# make cicero cds
umap_coords <- reducedDims(input_cds)$UMAP
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)


# read chromosome lengths
chrom_length <- read_tsv(dir_genomesizes, col_names = F)

# run cicero alternatively
distance_parameters <- estimate_distance_parameter(cicero_cds,
                                                   sample_num=100,
                                                   genomic_coords = chrom_length,
                                                   max_sample_windows = 1000)

model_output <- generate_cicero_models(cicero_cds,
                                       distance_parameter = distance_parameters %>% median(),
                                       genomic_coords = chrom_length, max_elements = 500)

conns <- assemble_connections(model_output)

# paste chromosome
conns[,1:2] <- conns[,1:2] %>% apply(2, function(x){paste0("chr", x)})

# saveRDS(conns, file="allcells_conns.RDS")


# visualization -----------------------------------------------------------
conns <-readRDS("allcells_conns.RDS")

# visualization
# export annotation
gene_anno <- readGFF(dir_gtf)

# rename some columns to match requirements
gene_anno$chromosome <- paste0("chr", gene_anno$seqid)
gene_anno$gene <- gene_anno$gene_id
gene_anno$transcript <- gene_anno$transcript_id
gene_anno$symbol <- gene_anno$gene_name

# plot list ---------------------------------------------------------------
# add enhancer stitched track
SE_bed <- read_tsv(dir_bed, col_names =  F) %>% mutate(X1 = paste0("chr", X1))

SE_bed <- makeGRangesFromDataFrame(SE_bed, start.field = "X2", end.field = "X3", seqnames.field = "X1")

dTrack3_SE <- AnnotationTrack(SE_bed,
)

# plots
plot_conns <- function(x, filename, start, end, chr, sizes = c(2,1,1,1)){
  svg(paste0(filename, ".svg"), width = 20, height = 14)
  Gviz::plotTracks(x,  
                   sizes = sizes,
                   from = start, to = end, chromosome = chr, 
                   transcriptAnnotation = "symbol",
                   col.axis = "black", 
                   fontsize = 40,
                   fontsize.group = 40,
                   fontcolor.legend = "black",
                   lwd=.3,
                   margin = 100, 
                   innerMargin = 50,
                   title.width = .5,
                   background.title = "transparent", 
                   col.border.title = "transparent")
  dev.off()
}

res_vis <- plot_connections(conns, paste0(region_to_vis[1]), region_to_vis[2], region_to_vis[3],
                 alpha_by_coaccess = TRUE,
                 gene_model = gene_anno,
                 coaccess_cutoff = 0.1,
                 connection_width = 1,
                 collapseTranscripts = "longest", 
                 return_as_list = T)


plot_conns(res_vis, "filename", 60658235, 60710751, "chr2")


# calculating number of tags per peak -------------------------------------
library(tidyverse)
library(Signac)
library(rtracklayer)

# peaks directory
dir_peaks <- "specify peaks directory here/peaks.bed"
# 10x fragment directory 
dir_fragments <- "fragments directory/fragments.tsv.gz"

# output directory
dir_out <- "output directory"

# read bed file as granges
bed_granges <- import(dir_peaks, format = "BED")

# assign fragments to matrix
  peakmatrix <- FeatureMatrix(CreateFragmentObject(dir_fragments),
                                   features = bed_granges, sep = c("_", "_"))
  
# save as RDS
saveRDS(peakmatrix, paste0(dir_out, "peak_matrix.RDS"))











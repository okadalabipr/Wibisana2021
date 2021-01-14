# Convert motifs from homer to fimo file -----------------------------------------------
​library(universalmotif)
library(tidyverse)
​
# input directory contains homer motif files
inputdir <- "directory of homer motif files"
dir.create("output directory", recursive = T)
outputdir <- "output directory"
​
# list input motif files
input_files <- list.files(inputdir) %>%
  str_subset("motif$") %>%
  str_remove_all("\\.motif$")
​
​# read homer motif files and write meme motif files
for(i in input_files){
  motif <- read_homer(paste0(inputdir, i, ".motif"))
  
  write_meme(motif, paste0(outputdir, i, ".meme"), version = 4, overwrite = T)
}


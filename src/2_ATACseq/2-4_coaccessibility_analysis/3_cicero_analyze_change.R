library(cicero)
library(tidyverse)
library(rtracklayer)

conns_stim <- readRDS("stim_conns.RDS")
conns_unstim <- readRDS("unstim_conns.RDS")

conns_stim <- conns_stim %>% filter(conns_stim$coaccess >= 0.1)

conns_joined <- left_join(conns_stim, conns_unstim, by = c("Peak1", "Peak2"))

conns_joined <- conns_joined %>% transmute(Peak1 = Peak1,
                                           Peak2 = Peak2,
                                           coaccess = coaccess.x-coaccess.y)

# filter by 0.05
conns_joined <- conns_joined %>% filter(coaccess >= 0.05)

# calculate distance between connections
conns_joined


conns_joined <- conns_joined %>%
  mutate(center_1 = str_split(Peak1, "_")[[1]] %>% .[2:3] %>% as.numeric() %>% mean(),
         center_2 = str_split(Peak2, "_")[[1]] %>% .[2:3] %>% as.numeric() %>% mean())

for (i in 1:nrow(conns_joined)) {
  conns_joined[i,"center_1"] <- conns_joined[i,"Peak1"] %>% str_split("_") %>% .[[1]] %>% .[2:3] %>% as.numeric() %>% mean()
  conns_joined[i,"center_2"] <- conns_joined[i,"Peak2"] %>% str_split("_") %>% .[[1]] %>% .[2:3] %>% as.numeric() %>% mean()

 
}

conns_joined <- conns_joined %>% mutate(dist = abs(center_1-center_2)) %>% as_tibble()

conns_bed <- conns_joined %>% select(Peak1, coaccess, dist)


conns_bed <- conns_bed %>%
  separate(c("Peak1"), c("chr", "start", "end"), "_")

# edit columns to match homer wanted file
conns_bed <- conns_bed %>% mutate(strand = "+", peakid = rownames(conns_bed)) %>% select(1:3, "peakid", "coaccess", "strand", "dist")

conns_bed <- conns_bed %>% mutate(chr = str_remove(chr, "chr"))

# make bed files of one of the peaks and the coaccessibility score
dir.create("conns_bed")
conns_bed %>% write_tsv("conns.bed", col_names = F)




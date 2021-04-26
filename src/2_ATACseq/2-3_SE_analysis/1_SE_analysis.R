library(tidyverse)
library(gridExtra)
library(VennDiagram)

TE <- read_tsv("TE_5000_combined_annotated_intensity.txt") %>% dplyr::select(1:4, 12, 16, 20:21)
SE <- read_tsv("SE_5000_combined_annotated_intensity.txt") %>% dplyr::select(1:4, 12, 16, 20:21)

colnames(TE) <- c("peakid", "chr", "start", "end", "ensembl_gene_id", "external_gene_name", "stim", "unstim")
colnames(SE) <- c("peakid", "chr", "start", "end", "ensembl_gene_id", "external_gene_name", "stim", "unstim")

# calculate l2fc
SE <- SE %>% mutate(l2fc = log2(stim/unstim))
TE <- TE %>% mutate(l2fc = log2(stim/unstim))

# add rank
SE <- SE %>% arrange(desc(l2fc)) %>% mutate(N = seq.int(nrow(.)))
TE <- TE %>% arrange(desc(l2fc)) %>% mutate(N = seq.int(nrow(.)))


# gained and lost SE and TE graphs
ggplot(SE,aes(x = N , y = l2fc, fill = -N))+
  geom_bar(stat = "identity", width = 1)+
  geom_hline(yintercept = c(SE %>% pull(l2fc) %>% quantile(0.25),SE %>% pull(l2fc) %>% quantile(0.75)), linetype="dashed", size = 2)+
  coord_flip()+
  scale_x_reverse()+
  guides(fill = F)+
  cowplot::theme_cowplot(40)+
  cowplot::panel_border("black", size = 2)+
  scale_fill_viridis_c()+
  ylab("") +
  xlab("") +
  ylab(element_blank())+
  xlab(element_blank())

ggsave("SE.svg", height = 7, width = 7)

ggplot(TE,aes(x = N , y = l2fc, fill = -N))+
  geom_bar(stat = "identity", width = 1)+
  geom_hline(yintercept = c(TE %>% pull(l2fc) %>% quantile(0.25),TE %>% pull(l2fc) %>% quantile(0.75)), linetype="dashed", size = 2)+
  coord_flip()+
  scale_x_reverse()+
  guides(fill = F)+
  cowplot::theme_cowplot(40)+
  cowplot::panel_border("black", size = 2)+
  scale_fill_viridis_c()+
  ylab(element_blank())+
  xlab(element_blank())


ggsave("TE.svg", height = 7, width = 7)

# calculate TE SE gain and lost
SE <- SE %>% 
  mutate(atac_state = case_when(l2fc > l2fc %>% quantile(0.75) ~ "atac_gain",
                                l2fc < l2fc %>% quantile(0.25) ~ "atac_lost",
                                l2fc <= l2fc %>% quantile(0.75) & l2fc >= l2fc %>% quantile(0.25) ~ "atac_unchanged"
  ))

TE <- TE %>% 
  mutate(atac_state = case_when(l2fc > l2fc %>% quantile(0.75) ~ "atac_gain",
                                l2fc < l2fc %>% quantile(0.25) ~ "atac_lost",
                                l2fc <= l2fc %>% quantile(0.75) & l2fc >= l2fc %>% quantile(0.25) ~ "atac_unchanged"
  ))

# remove TE in SE
TE <- TE %>% filter(!(ensembl_gene_id %in% SE$ensembl_gene_id))

# save SE and TE file
write_csv(SE, "SE.csv")
write_csv(TE, "TE.csv")

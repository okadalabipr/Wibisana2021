library(tidyverse)
library(cowplot)

# motif analysis for PU.1 and NF-kB, input file is from Homer output file
motif_list <- list.files("raw_data/") %>% str_subset("NFkB|PU1")

df_list <- list()

for (i in 1:length(motif_list)) {
  df_list[[i]] <- read_tsv(paste0("raw_data/", motif_list[i]))
  names(df_list)[i] <- motif_list[i] %>% str_remove(".txt")  
  df_list[[i]] <- df_list[[i]] %>% mutate(motif = motif_list[i] %>% str_extract("[^_]+"),
                                          atac = motif_list[i]  %>% str_remove("[^_]+") %>% str_remove(".txt") %>% str_remove("_"))
}

df_list <- bind_rows(df_list)
df_list <- df_list %>% separate(atac, into = c("Type", "atac_state") , sep = "_")
df_list <- df_list %>% group_split(motif,Type,atac_state) %>% as.list()

# calculate the number of motif occurences
df_list <- df_list %>% lapply(function(x){
  x %>% group_by(PositionID) %>% 
    summarise(occurences = n(), motif = motif, Type = Type, atac_state = atac_state) %>% 
    slice_head(n = 1) %>% ungroup()
})

df_list <- df_list %>% bind_rows() 

# add regions with no occurences
SE_TE_df <- bind_rows(
  SE_merged_raw %>% mutate(Type = "SE", atac_state = str_remove(atac_state, "atac_"), length=end-start) %>% dplyr::select(PositionID = peakid, Type, atac_state, length),
  TE_merged_raw %>% mutate(Type = "TE", atac_state = str_remove(atac_state, "atac_"), length=end-start) %>% dplyr::select(PositionID = peakid, Type, atac_state, length)
  
)

SE_TE_df <- bind_rows(
  SE_TE_df %>% mutate(motif = "NFkB"),
  SE_TE_df %>% mutate(motif = "PU1")
)

SE_TE_df <- SE_TE_df %>% mutate(occurences_0 = 0)

SE_TE_df <- SE_TE_df %>% mutate(new_id = paste(PositionID, motif, Type, atac_state, sep = "_"))

# join by new id
df_list <- df_list %>% mutate(new_id = paste(PositionID, motif, Type, atac_state, sep = "_"))

#sanity check if all true is ok
df_list$new_id %in% (left_join(SE_TE_df, df_list) %>% pull(new_id)) %>% table

df_list <- left_join(SE_TE_df, df_list)

# change NA to 0
df_list <- df_list %>% mutate_all(~replace(., is.na(.), 0))

# save as csv
df_list %>% mutate(occurencesbp = occurences/length*10000) %>% write_csv("results/motif_df.csv")

df_list %>% ggplot(aes(x = Type, y = occurences, fill = Type)) + geom_boxplot(outlier.shape = NA) +
  facet_grid(cols = vars(motif), scales="fixed")+ cowplot::theme_cowplot(25) + cowplot::panel_border(color = "black")+
  scale_fill_manual(values = c("#FCBBB6", "#80DFE2")) + 
  coord_cartesian(ylim = c(0, 125))+
  theme(legend.position = "none", strip.background = element_blank(), strip.text.y = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

dir.create("results/")
ggsave("results/motif_SE_TE.svg", height = 5, width = 7)

df_list %>% ggplot(aes(x = Type, y = occurences/length*10000, fill = Type)) + geom_boxplot(outlier.shape = NA) +
  facet_grid(cols = vars(motif), scales="fixed")+ cowplot::theme_cowplot(25) + cowplot::panel_border(color = "black")+
  scale_fill_manual(values = c("#FCBBB6", "#80DFE2")) + 
  coord_cartesian(ylim = c(0, 40))+
  theme(legend.position = "none", strip.background = element_blank(), strip.text.y = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

dir.create("results/")
ggsave("results/motif_SE_TE_perbp.svg", height = 5, width = 7)

# plot SE and TE difference
boxplot_motif <- function(df, motif_1, type, y1, y2){
  df <- df %>% mutate(atac_state = str_replace(atac_state, "gain", "gained")) %>%
    filter(motif == motif_1, Type == type) 
  
  df %>%
    ggplot(aes(x = atac_state, y = occurences, fill = atac_state)) + geom_boxplot(outlier.shape = NA) + 
    cowplot::theme_cowplot(25) + cowplot::panel_border(color = "black")+
    scale_fill_manual(values = c("#FCBBB6", "#80DFE2", "grey")) + 
    coord_cartesian(ylim = c(y1, y2))+
    theme(legend.position = "none", strip.background = element_blank(), strip.text.y = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
}


df_list %>% boxplot_motif("NFkB", "SE", 0, 60)
ggsave("results/nfkb_SE.svg", width = 5, height = 5)
df_list %>% boxplot_motif("NFkB", "TE", 0, 8)
ggsave("results/nfkb_TE.svg", width = 5, height = 5)
df_list %>% boxplot_motif("PU1", "SE", 0, 160)
ggsave("results/PU1_SE.svg", width = 5, height = 5)
df_list %>% boxplot_motif("PU1", "TE", 0, 30)
ggsave("results/PU1_TE.svg", width = 5, height = 5)

df_list %>%
  group_by(atac_state, motif, Type, ) %>% 
  summarise(median = median(occurences))

df_list %>% group_by(motif, Type, ) %>% 
  summarise(median = median(occurences))



#### t test and anova
minsample <- df_list %>% filter(motif == "NFkB") %>% dplyr::select(atac_state, Type) %>% table() %>% min()

aov(occurences ~ atac_state, data = df_list %>% filter(Type == "SE", motif == "NFkB") %>% group_by(atac_state) %>% slice_sample(n = minsample) %>% ungroup()) %>% summary

aov(occurences ~ atac_state, data = df_list %>% filter(Type == "TE", motif == "NFkB") %>% group_by(atac_state) %>% slice_sample(n = minsample) %>% ungroup()) %>% summary

# both
t.test(df_list %>% filter(Type == "SE", motif == "NFkB") %>% group_by(atac_state) %>% slice_sample(n = minsample) %>% ungroup() %>% pull(occurences),
       df_list %>% filter(Type == "TE", motif == "NFkB") %>% group_by(atac_state) %>% slice_sample(n = minsample) %>% ungroup() %>% pull(occurences),
       alternative = "two.sided", var.equal = FALSE)$p.value

t.test(df_list %>% filter(Type == "SE", motif == "PU1") %>% group_by(atac_state) %>% slice_sample(n = minsample) %>% ungroup() %>% pull(occurences),
       df_list %>% filter(Type == "TE", motif == "PU1") %>% group_by(atac_state) %>% slice_sample(n = minsample) %>% ungroup() %>% pull(occurences),
       alternative = "two.sided", var.equal = FALSE)$p.value


aov(occurences ~ atac_state, data = df_list %>% filter(Type == "SE", motif == "PU1") %>% group_by(atac_state) %>% slice_sample(n = minsample) %>% ungroup()) %>% summary
aov(occurences ~ atac_state, data = df_list %>% filter(Type == "TE", motif == "PU1") %>% group_by(atac_state) %>% slice_sample(n = minsample) %>% ungroup()) %>% summary
aov(occurences ~ atac_state, data = df_list %>% filter(Type == "SE", motif == "NFkB") %>% group_by(atac_state) %>% slice_sample(n = minsample) %>% ungroup()) %>% summary
aov(occurences ~ atac_state, data = df_list %>% filter(Type == "TE", motif == "NFkB") %>% group_by(atac_state) %>% slice_sample(n = minsample) %>% ungroup()) %>% summary



# motif per length --------------------------------------------------------
boxplot_motif_length <- function(df, motif_1, type, y1, y2){
  df <- df %>% mutate(atac_state = str_replace(atac_state, "gain", "gained")) %>% 
    mutate(atac_state = factor(atac_state, levels = c("gained", "unchanged", "lost"))) %>%
    filter(motif == motif_1, Type == type) 
  
  df %>%
    ggplot(aes(x = atac_state, y = occurences/length*10000, fill = atac_state)) + geom_boxplot(outlier.shape = NA) + 
    cowplot::theme_cowplot(25) + cowplot::panel_border(color = "black")+
    scale_fill_manual(values = c("#FCBBB6", "white", "#80DFE2")) + 
    coord_cartesian(ylim = c(y1, y2))+
    theme(legend.position = "none", strip.background = element_blank(), strip.text.y = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
}


df_list %>% boxplot_motif_length("NFkB", "SE", 0, 15)
ggsave("results/nfkb_SE_bp.svg", width = 5, height = 5)
df_list %>% boxplot_motif_length("NFkB", "TE", 0, 10)
ggsave("results/nfkb_TE_bp.svg", width = 5, height = 5)
df_list %>% boxplot_motif_length("PU1", "SE", 0, 25)
ggsave("results/PU1_SE_bp.svg", width = 5, height = 5)
df_list %>% boxplot_motif_length("PU1", "TE", 0, 50)
ggsave("results/PU1_TE_bp.svg", width = 5, height = 5)


df_list %>%
  group_by(atac_state, motif, Type, ) %>% 
  summarise(median = median(occurences/length*10000))

df_list %>% group_by(motif, Type, ) %>% 
  summarise(median = median(occurences/length*10000))

# t test and aov per length
df_list %>% mutate(occurencesbp = occurences/length*10000) %>% filter(Type == "SE", motif == "NFkB") %>% group_by(atac_state) %>% slice_sample(n = minsample) %>% ungroup() %>% pull(occurencesbp) %>% 
  mean()
df_list %>% mutate(occurencesbp = occurences/length*10000) %>% filter(Type == "TE", motif == "NFkB") %>% group_by(atac_state) %>% slice_sample(n = minsample) %>% ungroup() %>% pull(occurencesbp) %>% 
  mean()

t.test(df_list %>% mutate(occurencesbp = occurences/length) %>% filter(Type == "SE", motif == "NFkB") %>% group_by(atac_state) %>% slice_sample(n = minsample) %>% ungroup() %>% pull(occurencesbp),
       df_list %>% mutate(occurencesbp = occurences/length) %>% filter(Type == "TE", motif == "NFkB") %>% group_by(atac_state) %>% slice_sample(n = minsample) %>% ungroup() %>% pull(occurencesbp),
       alternative = "two.sided", var.equal = FALSE)$p.value

t.test(df_list %>% mutate(occurencesbp = occurences/length) %>% filter(Type == "SE", motif == "PU1") %>% group_by(atac_state) %>% slice_sample(n = minsample) %>% ungroup() %>% pull(occurencesbp),
       df_list %>% mutate(occurencesbp = occurences/length) %>% filter(Type == "TE", motif == "PU1") %>% group_by(atac_state) %>% slice_sample(n = minsample) %>% ungroup() %>% pull(occurencesbp),
       alternative = "two.sided", var.equal = FALSE)$p.value


aov(occurencesbp ~ atac_state, data = df_list %>% mutate(occurencesbp = occurences/length) %>% filter(Type == "SE", motif == "PU1") %>% group_by(atac_state) %>% slice_sample(n = minsample) %>% ungroup()) %>% summary
aov(occurencesbp ~ atac_state, data = df_list %>% mutate(occurencesbp = occurences/length) %>% filter(Type == "TE", motif == "PU1") %>% group_by(atac_state) %>% slice_sample(n = minsample) %>% ungroup()) %>% summary
aov(occurencesbp ~ atac_state, data = df_list %>% mutate(occurencesbp = occurences/length) %>% filter(Type == "SE", motif == "NFkB") %>% group_by(atac_state) %>% slice_sample(n = minsample) %>% ungroup()) %>% summary
aov(occurencesbp ~ atac_state, data = df_list %>% mutate(occurencesbp = occurences/length) %>% filter(Type == "TE", motif == "NFkB") %>% group_by(atac_state) %>% slice_sample(n = minsample) %>% ungroup()) %>% summary



# motif counting for NFKBIA and CD83 --------------------------------------
merged_geneid <- bind_rows(SE_merged_raw %>% mutate(motif = "NFkB"),
          SE_merged_raw %>% mutate(motif = "PU1")) %>% mutate(new_id = paste(peakid, motif, "SE", str_remove(atac_state, "atac_"), sep = "_"))


merged_geneid %>% inner_join(df_list, by = "new_id") %>% filter(external_gene_name %in% c("CD83", "NFKBIA")) %>% ggplot(aes(x = external_gene_name, y = (occurences/length)*10000)) +
  geom_bar(stat = "identity") + 
  cowplot::theme_cowplot(25) + cowplot::panel_border(color = "black")+
  scale_fill_manual(values = c("#FCBBB6", "#80DFE2", "grey")) + 
  facet_grid(cols = vars(motif.x)) + 
  theme(legend.position = "none", strip.background = element_blank(), strip.text.y = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())
  

ggsave("results/NFKBIA_CD83.svg", width = 7, height = 5)


merged_geneid %>% inner_join(df_list, by = "new_id") %>% filter(external_gene_name %in% c("CD83", "NFKBIA")) %>% ggplot(aes(x = motif.x, y = (occurences/length)*10000)) +
  geom_bar(stat = "identity") + 
  cowplot::theme_cowplot(25) + cowplot::panel_border(color = "black")+
  scale_fill_manual(values = c("#FCBBB6", "#80DFE2", "grey")) + 
  facet_grid(cols = vars(external_gene_name)) + 
  theme(legend.position = "none", strip.background = element_blank(), strip.text.y = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())


ggsave("results/NFKBIA_CD83_2.svg", width = 7, height = 5)


merged_geneid %>% inner_join(df_list, by = "new_id") %>% filter(external_gene_name %in% c("CD83", "NFKBIA")) %>% ggplot(aes(x = external_gene_name, y = occurences)) +
  geom_bar(stat = "identity") + 
  cowplot::theme_cowplot(25) + cowplot::panel_border(color = "black")+
  scale_fill_manual(values = c("#FCBBB6", "#80DFE2", "grey")) + 
  facet_grid(cols = vars(motif.x)) + 
  theme(legend.position = "none", strip.background = element_blank(), strip.text.y = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

ggsave("results/NFKBIA_CD83_2_not_normalized.svg", width = 7, height = 5)


library(tidyverse)

# read raw data -----------------------------------------------------------
# read counts
counts <- read_csv("baynorm_normalized.csv")

# read SE and TE
SE <- read_csv("SE.csv") %>% mutate(Type = "SE")
TE <- read_csv("TE.csv")%>% 
  filter(!ensembl_gene_id %in% SE$ensembl_gene_id) %>% mutate(Type = "TE")

SE_TE <- bind_rows(SE, TE)

# read DEGs
activated_markers <- read_csv("activated_markers.csv") %>% filter(p_val_adj < 0.05)

# calculate CV ------------------------------------------------------------
# convert to long df
counts_long <- counts %>% pivot_longer(starts_with("B")) %>% mutate(dose = case_when(str_detect(name, "^B00") == T ~ 0,
                                                                              str_detect(name, "^B01") == T ~ 0.01,
                                                                              str_detect(name, "^B02") == T ~ 0.1,
                                                                              str_detect(name, "^B03") == T ~ 1,
                                                                              str_detect(name, "^B04") == T ~ 10

                                                                               ))

# summarize cv by identifier
counts_long <- counts_long %>%
  group_by(Geneid, dose) %>%
  summarise(cv = sd(value)/mean(value), mean = mean(value), fano = var(value)/mean(value), sd = sd(value), var = var(value)) 
counts_long <- counts_long %>% pivot_wider(names_from = dose, values_from = c(cv, mean, fano, sd, var))
counts_long <- counts_long %>% na.omit()

# calculate fano ratio and cv ratio
counts_long <- counts_long %>% mutate(cv_ratio = cv_10/cv_0, fano_ratio = fano_10/fano_0)

# compare with fano factor ------------------------------------------------
counts_long %>% ggplot(aes(x = mean_10, y = fano_10)) + geom_point() + ggpubr::stat_cor(method = "spearman")
counts_long %>% ggplot(aes(x = fano_10, y = cv_10)) + geom_point() + ggpubr::stat_cor(method = "spearman")

# compare CV and Fano by using correlation --------------------------------
# correlation CV and fano at 10 ug/ml (DEGs)
counts_long %>% filter(Geneid %in% activated_markers$ensembl_gene_id) %>% 
  ggplot(aes(x = log2(cv_10), y = log2(fano_10))) + geom_point(size = 3, alpha = 0.5) +
  stat_smooth(method = "lm", col = "black") +
  ggpubr::stat_cor(method = "spearman", size = 10, label.y = 14)+
  cowplot::theme_cowplot(font_size = 40, line_size = 1)

# in all genes
counts_long %>% 
  ggplot(aes(x = log2(cv_10), y = log2(fano_10))) + geom_point(size = 3, alpha = 0.5) +
  stat_smooth(method = "lm", col = "black") +
  ggpubr::stat_cor(method = "spearman", size = 10, label.y = 14)+
  cowplot::theme_cowplot(font_size = 40, line_size = 1)


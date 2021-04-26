# Fano factor calculations

# read normalized read counts and convert id to external gene name
counts_normalized <- read_csv("baynorm_normalized.csv") %>% rename(ensembl_gene_id = "Geneid")

converted_id <- read_csv("converted_id.csv") %>% select(1:2)

counts_normalized <- inner_join(converted_id, counts_normalized, by = "ensembl_gene_id")

# read metadata
metadata <- read_csv("metadata_final.csv")

# extract inactivated and activated cells
counts_inactivated <- counts %>% select(external_gene_name, metadata %>% filter(sample == 0) %>% pull(X1))
counts_activated <- counts %>% select(external_gene_name, metadata %>% filter(sample == 10) %>% pull(X1))

counts_list <- list(counts_inactivated = counts_inactivated,
                    counts_activated = counts_activated)

# extract NFKBIA and CD83
counts_list <- lapply(counts_list, function(x){x %>% filter(external_gene_name %in% c("CD83", "NFKBIA"))})

counts_list <- list(CD83_act = counts_list$counts_activated %>% filter(external_gene_name == "CD83") %>% select(-1) %>% as.numeric(),
                    CD83_inact = counts_list$counts_inactivated %>% filter(external_gene_name == "CD83") %>% select(-1) %>% as.numeric(),
                    NFKBIA_act = counts_list$counts_activated %>% filter(external_gene_name == "NFKBIA") %>% select(-1) %>% as.numeric(),
                    NFKBIA_inact = counts_list$counts_inactivated %>% filter(external_gene_name == "NFKBIA") %>% select(-1) %>% as.numeric())


# read rna-fish data
CD83_fish <- read_csv("CD83.csv")
NFKBIA_fish <- read_csv("NFKBIA.csv")

# filter dose 10 and 0
rnafish_list <- list(CD83_act = CD83_fish %>% filter(dose == 10) %>% pull(N_Q570),
                 CD83_inact = CD83_fish %>% filter(dose == 0) %>% pull(N_Q570),
                 NFKBIA_act = NFKBIA_fish %>% filter(dose == 10) %>% pull(N_Q570),
                 NFKBIA_inact = NFKBIA_fish %>% filter(dose == 0) %>% pull(N_Q570))




# calculate CV for CD83 and NFKBIA ----------------------------------------
mean_list <- lapply(counts_list, mean)
sd_list <- lapply(counts_list, sd)

cv_df <- bind_cols(mean_list %>% as_tibble() %>% t() %>% as.data.frame(), sd_list %>% as_tibble() %>% t() %>% as.data.frame()) %>% rownames_to_column("data")

colnames(cv_df) <- c("data", "mean", "sd")

cv_df <- cv_df %>% mutate(cv = sd/mean)

# calculate CV for RNA-FISH
rnafish_mean_list <- lapply(rnafish_list, mean)
rnafish_sd_list <- lapply(rnafish_list, sd)

rnafish_cv_df <- bind_cols(rnafish_mean_list %>% as_tibble() %>% t() %>% as.data.frame(), rnafish_sd_list %>% as_tibble() %>% t() %>% as.data.frame()) %>% rownames_to_column("data")

colnames(rnafish_cv_df) <- c("data", "mean", "sd")

rnafish_cv_df <- rnafish_cv_df %>% mutate(cv = sd/mean)

# calculate Fano factor (var/mean)
cv_df <- cv_df %>% mutate(noise_strength = sd^2 / mean)
rnafish_cv_df <- rnafish_cv_df %>% mutate(noise_strength = sd^2 / mean)

# bind_cols(cv_df, rnafish_cv_df) %>% write_csv("dose_0_10.csv")

# normalization and comparison of distribution ----------------------------
# create DF for plotting
counts_df <- c()

for(i in 1:length(counts_list)){
  counts_df[[i]] <- counts_list[[i]] %>%
    data.frame() %>%
    mutate(state = names(counts_list)[i]) %>%
    mutate(data = "RNA-seq")
}

# calculate z score per gene
counts_df <- lapply(counts_df, function(x){
  scaled_values <- x[,1] %>% scale() %>% as.vector()
  
  x$scaled_value <- scaled_values  
  
  return(x)
  
})


counts_df <- do.call(bind_rows, counts_df) %>% as_tibble()

colnames(counts_df)[1] <- "value"


# calculate z score and make new dataframe for RNA-FISH data
rnafish_df <- c()

for(i in 1:length(rnafish_list)){
  rnafish_df[[i]] <- rnafish_list[[i]] %>%
    data.frame() %>%
    mutate(state = names(rnafish_list)[i]) %>%
    mutate(data = "RNA-FISH")
}

# calculate z score for RNA-FISH
rnafish_df <- lapply(rnafish_df, function(x){
  scaled_values <- x[,1] %>% scale() %>% as.vector()
  
  x$scaled_value <- scaled_values  
  
  return(x)
  
})


rnafish_df <- do.call(bind_rows, rnafish_df) %>% as_tibble()

colnames(rnafish_df)[1] <- "value"


# combine both data frames
all_df <- bind_rows(counts_df, rnafish_df)

# change labels
all_df <- all_df %>% mutate(state = case_when(state == "NFKBIA_inact" ~ "NFKBIA/0",
                                    state == "NFKBIA_act" ~ "NFKBIA/10",
                                    state == "CD83_act" ~ "CD83/10",
                                    state == "CD83_inact" ~ "CD83/0",
                                    ))


# calculate z score per gene
z_score_per_gene_1 <- all_df %>% filter(state == "CD83/10" | state == "CD83/0") %>% filter(data == "RNA-seq") %>% pull(value) %>% scale() %>% as.vector()
z_score_per_gene_2 <- all_df %>% filter(state == "NFKBIA/10" | state == "NFKBIA/0") %>% filter(data == "RNA-seq") %>% pull(value) %>% scale() %>% as.vector()
z_score_per_gene_3 <- all_df %>% filter(state == "CD83/10" | state == "CD83/0") %>% filter(data == "RNA-FISH") %>% pull(value) %>% scale() %>% as.vector()
z_score_per_gene_4 <- all_df %>% filter(state == "NFKBIA/10" | state == "NFKBIA/0") %>% filter(data == "RNA-FISH") %>% pull(value) %>% scale() %>% as.vector()


# return to data frame
all_df$scaled_value_per_gene <- c(z_score_per_gene_1, z_score_per_gene_2, z_score_per_gene_3, z_score_per_gene_4)

# plot ridgeplot
all_df$data <- all_df$data %>% as.factor() %>% relevel("RNA-seq")

ggplot(all_df, aes(x = scaled_value_per_gene, y = state,
                   fill = fct_rev(data)))+
  geom_density_ridges(alpha = .5)+
  #  scale_fill_manual(values = data, name = 'Party',labels = c('Democrat','Republican'))+
  labs(y = element_blank(),
       x = element_blank(),
       title = element_blank()) + cowplot::theme_cowplot(font_size = 30, line_size = 2) + theme(legend.position="top", legend.title = element_blank())

ggsave("ridgeplot_scaled_per_gene.svg", width = 8, height = 12)



# fano factor and CV calculation for all data at 10 and 0 ug/ml (baynorm normalized)
count_fano <- function(cell_10_df, cell_0_df){
  
  # calculate cv
  cell_10_df <- cell_10_df %>% mutate(mean = rowMeans(.[-1:-2]), sd = rowSds(as.matrix(.[-1:-2]))) %>% mutate(cv = (sd/mean)*100, fano = (sd^2/mean)) %>% dplyr::select(ensembl_gene_id, external_gene_name, mean_10 = mean, sd_10 = sd, cv_10 = cv, fano_10 = fano)
  cell_0_df <- cell_0_df %>% mutate(mean = rowMeans(.[-1:-2]), sd = rowSds(as.matrix(.[-1:-2]))) %>% mutate(cv = (sd/mean)*100, fano = (sd^2/mean)) %>% dplyr::select(ensembl_gene_id, external_gene_name, mean_0 = mean, sd_0 = sd, cv_0 = cv, fano_0 = fano)
  
  # combine onto one data frame
  cv_df <- inner_join(cell_10_df, cell_0_df, by = c("external_gene_name", "ensembl_gene_id"))
  
  return(cv_df)
}

bay_out_counts_dose <- list(
  dose_10 = bay_out_counts_wide %>% select(1:2, contains("B04_")) ,
  dose_0 = bay_out_counts_wide %>% select(1:2, contains("B00_"))
)

bay_out_counts_dose_fano <- count_fano(bay_out_counts_dose$dose_10,
           bay_out_counts_dose$dose_0) %>% na.omit()


write_csv(bay_out_counts_dose_fano, "baynorm_fano.csv")
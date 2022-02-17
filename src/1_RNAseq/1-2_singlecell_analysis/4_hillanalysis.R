library(tidyverse)
library(optimx)

# functions
hillfunc <-function(X){
  (param$p1)*X^param$p2 / (param$p3^param$p2 + X^param$p2)
}

fitting_func <- function(X, Y, n_predict = 2, change_predict = 0.1){
  require(optimx)
  require(data.table)
  require(tidyverse)
  
  X <- X
  Y <- Y
  
  fun_estimate <- function (Y, X) { # Y = gene expression data of each dose or time-point
    resid <- function (param) {
      yhat <- (param[1])*X^param[2] / (param[3]^param[2] + X^param[2]) # Hill equation
      sum((Y-yhat)^2)	# return Residual sum of squares (RSS)
    }
    
    #omptimize
    initialp <- c(min(Y), max(Y), n_predict, change_predict)  #initial parameteres of parameter, 3rd argument = n, 4th argument = where the change is expected to happen
    RESULT <- optimx(par = initialp, fn = resid, control = list(all.methods = TRUE)) #you can constrain range of parameteres lower = or upper =
    RESULT %>% data.table(keep.rownames=TRUE) %>% 
      dplyr::filter(convcode==0) %>% # filter(min(p1,p2,p3,p4)>0) %>% fr(p4<12) filtering of the results if you want
      dplyr::filter(value == min(value)) -> RESULT_OPT
    return(RESULT_OPT)
  }
  
  # make hill function
  hillfunc <-function(X, param){
    (param$p1)*X^param$p2 / (param$p3^param$p2 + X^param$p2)
  }
  
  # perform parameter estimation
  res <- fun_estimate(Y, X)
  res <- res[1]
  
  # extract only one if there are 2 best 
  fitted_df <- tibble(X = X, Y = Y, fitted_Y = hillfunc(X, res[1]))
  
  # extract parameters and create df
  res_param <- tibble(param = c("n", "ka", "kd"), value = c(res$p2, res$p3, res$p3^res$p2))
  
  list(parameters = res,
       params_short = res_param,
       fitted = fitted_df) %>% return()
  
}

#  read files -------------------------------------------------------------
# normalized table
normalized_counts <- read_csv("baynorm normalized counts")
# remove counts below 10
normalized_counts <- normalized_counts[(normalized_counts %>% select(-1) %>% rowSums()) > 10,]
# read fano calculations
fano_raw <- read_csv("matrix containing fano factor calculation for every genes per dose")
# converted id list
converted_id <- read_csv("csv of converted ID (ensembl id, external gene name)")
# enhancer list
SE <- read_csv("csv of SE list")
TE <- read_csv("csv of TE list")

# preprocess data --------------------------------------------------------------
# group by the cell identifier for each dose
expression <- normalized_counts %>% dplyr::rename(ensembl_gene_id = "Geneid") %>%
  pivot_longer(cols = starts_with(c("B0"))) %>%
  dplyr::select(name, value, ensembl_gene_id)

# add dose based on the column "name"
expression <- expression %>% mutate(dose = case_when(str_detect(name, "B00") == T ~ 0,
                                                     str_detect(name, "B01") == T ~ 0.01,
                                                     str_detect(name, "B02") == T ~ 0.1,
                                                     str_detect(name, "B03") == T ~ 1,
                                                     str_detect(name, "B04") == T ~ 10
))

# summarise to get the mean expression per dose
expression_mean <- expression %>% group_by(ensembl_gene_id, dose) %>% 
  summarise(mean = mean(value)) %>%
  ungroup()

# calculate fold change compared to dose 0
expression_mean_list <- expression_mean %>% group_split(ensembl_gene_id)

expression_mean_list <- lapply(expression_mean_list, function(x){
  dose0_mean <- x %>% filter(dose == 0) %>% pull(mean)
  df <- x %>% mutate(fc = (mean - dose0_mean)/dose0_mean)
  return(df)  
})

expression_mean <- expression_mean_list %>% bind_rows

# remove genes without gained TE or SE
expression_mean <- expression_mean %>% 
  filter(ensembl_gene_id %in% c(TE %>% filter(atac_state == "atac_gain") %>% pull(ensembl_gene_id), 
                                SE %>% filter(atac_state == "atac_gain") %>% pull(ensembl_gene_id)))

# perform fitting per genes -----------------------------------------------
# get unique gene set
gene_set <- expression_mean %>% pull(ensembl_gene_id) %>% unique()

# perform fitting and compute parameters for each gene 
opt_res <- list()

for(i in 1:length(gene_set)){
  tryCatch({
    exp <- expression_mean %>%
      filter(ensembl_gene_id == gene_set[i]) %>%
      arrange(dose) %>% pull(fc)
        
    fitted_vals <- fitting_func(X = c(0, 0.01, 0.1, 1, 10),
                                Y = exp,
                                change_predict = 0.1, n_predict = 2)
     
    opt_res[[i]] <- fitted_vals
    
    names(opt_res)[i] <- gene_set[i]
    
    print(paste0(i, " out of ", length(gene_set), " done"))
        
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  
}

# save plots
for (i in 1:length(opt_res)){
  param <- opt_res[[i]]$parameters
  
  scaleFUN <- function(x) {x}
  
  df <- expression_mean %>% filter(ensembl_gene_id == gene_set[i]) 
  
  df$dose[df$dose == 0] <- 0.001
  
  df %>%  ggplot(aes(x = dose, y = fc)) +  
    stat_function(fun = hillfunc, size=3, color= "black") +
    geom_point(color = "red", size = 4)+
    # geom_point(data = df %>% group_by(dose) %>% summarise(mean = mean(fc)), aes(x = dose, y = fc), color = "red", size = 4) +
    scale_x_log10(labels = scaleFUN) + 
    cowplot::theme_cowplot(30) +
    theme(legend.position = "none") +
    labs(title=element_blank(), x = element_blank(), y = element_blank())
  
  ggsave(paste0("output/", opt_res[[i]]$parameters[[3]],
                "_", converted_id %>% filter(ensembl_gene_id == gene_set[i]) %>%
                  pull(external_gene_name), ".svg"), width = 7,height=7)
  
}


# store parameters in a dataframe ----------------------------------------------------

# extract n ka kd only
opt_res_params <- lapply(opt_res, function(x){
  x$params_short
})

# add another column of gene id
for (i in 1:length(opt_res_params)) {
  opt_res_params[[i]]$ensembl_gene_id <- names(opt_res_params)[i]
}

opt_res_params <- bind_rows(opt_res_params) %>% left_join(converted_id) %>% left_join(fano_clusters %>% select(1,2))
opt_res_params %>% write_csv("optimized_params.csv")

# calculate mean per cluster
opt_res_params_summary <- opt_res_params %>% filter(param == "n") %>% group_by(cluster) %>% 
  summarise(mean = mean(value), sd = sd(value)) %>% ungroup

# make df of fitted values
opt_res_fitted <- c()
for (i in 1:length(opt_res)) {
  opt_res_fitted[[i]] <- opt_res[[i]]$fitted
  opt_res_fitted[[i]]$ensembl_gene_id <- names(opt_res)[i]
  opt_res_fitted[[i]]$val <- "fitted"
}

opt_res_fitted <- opt_res_fitted %>% bind_rows %>% left_join(converted_id)

# extract all parameters
opt_res_params_full <- c()
for (i in 1:length(opt_res)) {
  opt_res_params_full[[i]] <- opt_res[[i]]$parameters
  opt_res_params_full[[i]]$ensembl_gene_id <- names(opt_res)[i]
}

opt_res_params_full <- bind_rows(opt_res_params_full) %>% filter(convcode == 0)

opt_res_params_full %>% write_csv("optimized_params_full.csv")

# remove non meaningful genes -------------------------------------------------
# genes with p1 < 0 are decreasing and thus not meaningful
# lfc > 0.2 (dose 0 and 10)
# remove genes with p1 < 0
# remove hill outliers p2 < 9, p2 > 0.3

filtered_genes <- opt_res_params_full %>% filter(p1 > 0) %>% filter(p2 < 9 & p2 > 0.3) %>%
  filter(ensembl_gene_id %in% (expression_mean %>% filter(dose == 10, fc > 0.2) %>% pull(ensembl_gene_id) )) %>%
  inner_join(converted_id) 

filtered_genes %>% write_csv("filtered_params.csv")

# boxplot -------------------------------------------------
filtered_params <- filtered_genes %>% as_tibble() %>% mutate(atac = case_when(ensembl_gene_id %in% TE$ensembl_gene_id ~ "TE_gained",
                                                           ensembl_gene_id %in% SE$ensembl_gene_id ~ "SE_gained",
                                                           T ~ "none"
                                                           )) %>% filter(atac != "none")

filtered_params %>% ggplot(aes(x = atac, y = p2, fill = atac)) + geom_boxplot(outlier.shape = NA) +
  cowplot::theme_cowplot(25) + cowplot::panel_border(color = "black")+
  scale_fill_manual(values = c("#FCBBB6", "#80DFE2")) + 
  coord_cartesian(ylim = c(0, 10))+
  theme(legend.position = "none", strip.background = element_blank(), strip.text.y = element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank())

# ggsave("boxplot.svg", height = 5, width = 5)


# t test
#t.test(filtered_params %>% filter(atac == "SE_gained") %>% pull(p2),
#       filtered_params %>% filter(atac == "TE_gained") %>% pull(p2),
#       alternative = "two.sided", var.equal = FALSE)$p.value

# separate hill coefficient
filtered_params <- filtered_params  %>% mutate(hill_cat = case_when(p2 <= 1 ~ "low",
                                                                   p2 > 1 & p2 < 5 ~ "high",
                                                                   p2 >= 5 ~ "1_high"
                                                                   ))

# make bar plot
filtered_params %>% ggplot(aes(x = atac, y = 100, fill = hill_cat)) + 
  geom_bar(stat="identity", position = "fill") + cowplot::theme_cowplot(25) + cowplot::panel_border(color = "black") +
  theme(legend.position = "none", axis.title.x=element_blank(), axis.title.y=element_blank())

# ggsave("barplot.svg", height = 5, width = 5)





















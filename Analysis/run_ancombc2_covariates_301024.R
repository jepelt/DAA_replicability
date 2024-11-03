library(tidyverse)

#Helper functions---------------------------------------------------------------
make_physeq <- function(counts, meta){
  otu <- phyloseq::otu_table(t(counts), taxa_are_rows = TRUE)
  meta_data <- phyloseq::sample_data(meta)
  phyloseq::phyloseq(otu, meta_data)
}


#ANCOM-BC2----------------------------------------------------------------------
run_ancombc2 <- function(counts, meta, fm = ~ group, p_adj = 'BH'){
  
  physeq <- make_physeq(counts, meta)
  fmc <- paste(as.character(fm)[-1], collapse = ' + ')
  
  obj_p <- NULL
  try(obj_p <- ANCOMBC::ancombc2(data = physeq, fix_formula = fmc,
                                 prv_cut = 0, p_adj_method = 'none',
                                 verbose = F, alpha = .05))
  
  obj_q01 <- NULL
  try(obj_q01 <- ANCOMBC::ancombc2(data = physeq, fix_formula = fmc,
                                   prv_cut = 0, p_adj_method = p_adj,
                                   verbose = F, alpha = .01))
  
  obj_q05 <- NULL
  try(obj_q05 <- ANCOMBC::ancombc2(data = physeq, fix_formula = fmc,
                                   prv_cut = 0, p_adj_method = p_adj,
                                   verbose = F, alpha = .05))
  
  obj_q10 <- NULL
  try(obj_q10 <- ANCOMBC::ancombc2(data = physeq, fix_formula = fmc,
                                   prv_cut = 0, p_adj_method = p_adj,
                                   verbose = F, alpha = .10))
  
  obj_q20 <- NULL
  try(obj_q20 <- ANCOMBC::ancombc2(data = physeq, fix_formula = fmc,
                                   prv_cut = 0, p_adj_method = p_adj,
                                   verbose = F, alpha = .20))
  
  if(!is.null(obj_p) & !is.null(obj_q01) & !is.null(obj_q05) & 
     !is.null(obj_q10) & !is.null(obj_q20)){
    
    r <- obj_p$res
    res <- tibble(taxon = r$taxon,
                  est = r$lfc_groupcase,
                  se = r$se_groupcase,
                  df = nrow(counts) - length(all.vars(fm)) - 1,
                  p0 = r$p_groupcase,
                  ss_pass_p = r$passed_ss_groupcase,
                  q01 = obj_q01$res$q_groupcase,
                  ss_pass_q01 = obj_q01$res$passed_ss_groupcase,
                  q05 = obj_q05$res$q_groupcase,
                  ss_pass_q05 = obj_q05$res$passed_ss_groupcase,
                  q10 = obj_q10$res$q_groupcase,
                  ss_pass_q10 = obj_q10$res$passed_ss_groupcase,
                  q20 = obj_q20$res$q_groupcase,
                  ss_pass_q20 = obj_q20$res$passed_ss_groupcase,
                  method = paste0('ANCOM-BC2 (', p_adj, ')')) %>% 
      mutate(p = ifelse(p0 < 0.05 & !ss_pass_p, runif(length(p0), 0.05, 1), p0),
             q01 = ifelse(ss_pass_q01, q01, 1),
             q05 = ifelse(ss_pass_q05, q05, 1),
             q10 = ifelse(ss_pass_q10, q10, 1),
             q20 = ifelse(ss_pass_q20, q20, 1)
      ) %>% 
      select(taxon, est, se, df, p, method, q01, q05, q10, q20,
             ss_pass_p, ss_pass_q01, ss_pass_q05, ss_pass_q10, ss_pass_q20)
    
  }else{
    res <- tibble(taxon = colnames(counts),
                  method = paste0('ANCOM-BC2 (', p_adj, ')'))
  }
  
  return(res %>% mutate(data_id = meta$data_id[1]))
}

#Run methods--------------------------------------------------------------------
load('data_171023.rds')

create_formula <- function(d){
  bmi <- ifelse(d$meta$is_bmi[1], '+ bmi', '')
  age <- ifelse(d$meta$is_age[1], '+ age', '')
  sex <- ifelse(d$meta$is_sex[1], '+ sex', '')
  
  fm <- paste0('~ group', bmi, age, sex) %>% as.formula
  return(fm)
}

datasets <- keep(d_all, ~ create_formula(.) != ~ group)

lcounts <- map(datasets, ~ .$counts)
lmeta <- map(datasets, ~ .$meta)
lformula <- map(datasets, create_formula)

future::plan(future::multisession, workers = 7)

print(Sys.time())
resl_ancombc2 <- furrr::future_pmap(list(lcounts, lmeta, lformula),
                                    ~ run_ancombc2(..1, ..2, ..3),
                                    .options = furrr::furrr_options(seed = 1),
                                    .progress = T)
print(Sys.time())

res_cov_ancombc2 <- bind_rows(resl_ancombc2)
#save(res_cov_ancombc2, file = 'res_ancombc2_cov_301024.rds')

library(tidyverse)
library(data.table)

#Load files and prepare tibble of results---------------------------------------
load('res_main.rds')
load('res_ancombc_251024.rds')
load('data_meta_041023.rds')
load('data_prevl_061023.rds')


#Split-data-analyses------------------------------------------------------------
res_all <- bind_rows(res_list) %>% 
  bind_rows(., bind_rows(resl_ancombc2)) %>% 
  left_join(., d_meta, by = 'data_id') %>% 
  left_join(., d_prevl, by = c('data_id', 'taxon')) %>% 
  as.data.table() %>% 
  .[method %in% c('ALDEx2 (glm)',
                  'ANCOM-BC2 (BH)',
                  'corncob (LRT, ev = TRUE)',
                  'DESeq2 (Wald, Default)',
                  'edgeR (QL, TMM)',
                  'fastANCOM (0.5, 0.05, 0.05)',
                  'LDM (CLR = FALSE)',
                  'limma-voom (TMM)',
                  'LinDA (Adaptive, 0.5)',
                  'LogR (Firth)',
                  'MaAsLin2 (LOG, TSS)',
                  'mSeq (CSS)',
                  'ORMW (Score, TSS)',
                  'ZicoSeq (Default)')] %>% 
  .[, method := ifelse(method == 'mSeq (CSS)', 'metagenomeSeq', method)] %>% 
  .[, method := ifelse(method == 'ORMW (Score, TSS)', 'ORM/Wilcoxon', method)] %>% 
  .[, method := ifelse(method == 'MaAsLin2 (LOG, TSS)', 'MaAsLin2/t-test', method)] %>% 
  .[, method := ifelse(method == 'LogR (Firth)', 'LogR', method)] %>% 
  .[, method := tstrsplit(method, " ", fixed = T, keep = 1)] %>% 
  .[, analysis := tstrsplit(data_id, " ", fixed = T, keep = 1)] %>% 
  .[, q := p.adjust(p, method = 'BH'), by = c('data_id', 'method')] %>% 
  .[, q := ifelse(grepl('DESeq2|LDM|Zico', method, perl = T), p_fdr, q)] %>%
  .[, t91.7 := qt(.834 + (1 - .834) / 2, df = df)] %>% 
  .[, lwr := ifelse(is.na(lwr), est - t91.7 * se, lwr)] %>% 
  .[, upr := ifelse(is.na(upr), est + t91.7 * se, upr)] %>% 
  .[, disease := ifelse(disease == 'EDD', 'CDI', disease)] %>% 
  filter(site == 'stool')


#Summarize results--------------------------------------------------------------
alphas <- c(.05)
res_summary <- list()

for(i in 1:length(alphas)){
  
  res_exploratory <- res_all %>% 
    filter(analysis == 'half' & half == 1) %>% 
    mutate(q = case_when(method == 'ANCOM-BC2' & alphas[i] == .01 ~ q01,
                         method == 'ANCOM-BC2' & alphas[i] == .05 ~ q05,
                         method == 'ANCOM-BC2' & alphas[i] == .10 ~ q10,
                         method == 'ANCOM-BC2' & alphas[i] == .20 ~ q20,
                         T ~ q)) %>% 
    filter(q < alphas[i]) %>% 
    select(method, type, study, iter, taxon, est1 = est, q)
  
  res_validation <- res_all %>% 
    filter(analysis == 'half' & half == 2) %>% 
    mutate(is2 = T) %>% 
    select(method, study, iter, taxon, is2, est2 = est, p)
  
  res_combined <- left_join(res_exploratory, res_validation,
                            by = c('method', 'study', 'iter', 'taxon')) %>% 
    filter(!is.na(is2)) %>% 
    mutate(sgn1 = T,
           d1 = ifelse(is.na(est1) | is.nan(est1), 0, sign(round(est1, 8))),
           s1 = if_else(sgn1 & d1 != 0, T, F, F),
           sgn2 = p < .05,
           d2 = ifelse(is.na(est2) | is.nan(est2), 0, sign(round(est2, 8))),
           s2 = if_else(sgn2 & d2 != 0, T, F, F),
           nores = d2 == 0,
           repl = s1 & s2 & d1 == d2,
           same = s1 & !s2 & d1 == d2 & !nores,
           opposite = s1 & !s2 & d1 != d2 & !nores,
           anti = s1 & s2 & d1 != d2) %>%
    mutate_if(is.logical, as.numeric) %>% 
    mutate(check = repl + same + opposite + anti + nores) %>%
    select(method, type, study, iter, taxon,
           repl, same, opposite, anti, nores, check)
  
  stopifnot(all(res_combined$check == 1))
  
  hits <- res_exploratory %>% group_by(method) %>% summarize(r_nhit = n())
  
  res_summary[[i]] <- res_combined %>% 
    group_by(method) %>% 
    summarize_at(vars(repl, same, opposite, anti, nores),
                 ~ sum(., na.rm = T)) %>% 
    ungroup() %>% 
    left_join(., hits, 'method') %>% 
    mutate(n = repl + same + opposite + anti + nores,
           r_confl = anti / n * 100,
           r_repl = repl / n * 100,
           alpha = alphas[i]) %>% 
    select(method, r_confl, r_repl, r_nhit, alpha)
}

res_split <- bind_rows(res_summary) %>% mutate(analysis = 'split')


#Separate study analyses--------------------------------------------------------
res_all <- bind_rows(res_list) %>% 
  bind_rows(., bind_rows(resl_ancombc2)) %>% 
  left_join(., d_meta, by = 'data_id') %>% 
  left_join(., d_prevl, by = c('data_id', 'taxon')) %>% 
  as.data.table() %>% 
  .[method %in% c('ALDEx2 (glm)',
                  'ANCOM-BC2 (BH)',
                  'corncob (LRT, ev = TRUE)',
                  'DESeq2 (Wald, Default)',
                  'edgeR (QL, TMM)',
                  'fastANCOM (0.5, 0.05, 0.05)',
                  'LDM (CLR = FALSE)',
                  'limma-voom (TMM)',
                  'LinDA (Adaptive, 0.5)',
                  'LogR (Firth)',
                  'MaAsLin2 (LOG)',
                  'mSeq (CSS)',
                  'ORM (Score, TSS)',
                  'ZicoSeq (Default)')] %>% 
  .[, method := ifelse(method == 'mSeq (CSS)', 'metagenomeSeq', method)] %>% 
  .[, method := ifelse(method == 'ORM (Score, TSS)', 'ORM/Wilcoxon', method)] %>% 
  .[, method := ifelse(method == 'MaAsLin2 (LOG)', 'MaAsLin2/t-test', method)] %>% 
  .[, method := ifelse(method == 'LogR (Firth)', 'LogR', method)] %>% 
  .[, method := tstrsplit(method, " ", fixed = T, keep = 1)] %>% 
  .[, analysis := tstrsplit(data_id, " ", fixed = T, keep = 1)] %>% 
  .[, q := p.adjust(p, method = 'BH'), by = c('data_id', 'method')] %>% 
  .[, q := ifelse(grepl('DESeq2|LDM|Zico', method, perl = T), p_fdr, q)] %>%
  .[, t91.7 := qt(.834 + (1 - .834) / 2, df = df)] %>% 
  .[, lwr := ifelse(is.na(lwr), est - t91.7 * se, lwr)] %>% 
  .[, upr := ifelse(is.na(upr), est + t91.7 * se, upr)] %>% 
  .[, disease := ifelse(disease == 'EDD', 'CDI', disease)] %>% 
  filter(site == 'stool') %>% 
  mutate(n0 = n + avg_ls / 10 ^ 10,
         is2 = T) %>% 
  filter(analysis == 'whole' & n_case >= 10 & n_control >= 10) %>% 
  filter(!disease %in% c('HIV', 'schizophrenia', 'NASH'))


#summarize results--------------------------------------------------------------
alphas <- c(.05)
res_summary <- list()

expl_studies <- left_join(
  res_all %>% select(type, disease, taxon, study, n1 = n0, data_id),
  res_all %>% select(type, disease, taxon, n2 = n0, is2),
  by = join_by(type, disease, taxon, n1 < n2)) %>% 
  filter(!is.na(is2)) %>% 
  distinct(study, .keep_all = T) %>% 
  pull(data_id)

for(i in 1:length(alphas)){
  
  res_exploratory <- res_all %>% 
    filter(data_id %in% expl_studies) %>% 
    mutate(q = case_when(method == 'ANCOM-BC2' & alphas[i] == .01 ~ q01,
                         method == 'ANCOM-BC2' & alphas[i] == .05 ~ q05,
                         method == 'ANCOM-BC2' & alphas[i] == .10 ~ q10,
                         method == 'ANCOM-BC2' & alphas[i] == .20 ~ q20,
                         T ~ q)) %>% 
    filter(q < alphas[i]) %>% 
    select(method, type, disease, taxon, study1 = study, n1 = n0, est1 = est, q)
    
  res_validation <- res_all %>% 
    filter(analysis == 'whole') %>% 
    select(method, type, disease, taxon, study2 = study, n2 = n0, is2,
           est2 = est, p)
  
  res_combined <- left_join(res_exploratory, res_validation,
                            by = join_by(method, type, disease, taxon, n1 < n2)) %>% 
    filter(!is.na(is2)) %>% 
    mutate(sgn1 = T,
           d1 = ifelse(is.na(est1) | is.nan(est1), 0, sign(round(est1, 8))),
           s1 = if_else(sgn1 & d1 != 0, T, F, F),
           sgn2 = p < .05,
           d2 = ifelse(is.na(est2) | is.nan(est2), 0, sign(round(est2, 8))),
           s2 = if_else(sgn2 & d2 != 0, T, F, F),
           nores = d2 == 0,
           repl = s1 & s2 & d1 == d2,
           same = s1 & !s2 & d1 == d2 & !nores,
           opposite = s1 & !s2 & d1 != d2 & !nores,
           anti = s1 & s2 & d1 != d2) %>%
    mutate_if(is.logical, as.numeric) %>% 
    mutate(check = repl + same + opposite + anti + nores) %>%
    select(method, type, disease, study1, study2, taxon,
           repl, same, opposite, anti, nores, check)
  
  stopifnot(all(res_combined$check == 1))
  stopifnot(nrow(res_combined %>% distinct(study1)) == 37)
  stopifnot(nrow(res_combined %>% distinct(study2)) == 37)
  
  n_studies <- length(expl_studies)
  
  hits <- res_exploratory %>% 
    group_by(method) %>% 
    summarize(r_nhit = n())
  
  res_summary[[i]] <- res_combined %>% 
    group_by(method, study1, taxon) %>% 
    summarize_at(vars(repl, same, opposite, anti, nores), ~ sum(., na.rm = T)) %>% 
    mutate(n = repl + same + opposite + anti + nores) %>% 
    mutate_at(vars(repl, same, opposite, anti, nores), ~ . / n) %>% 
    group_by(method) %>% 
    summarize_at(vars(repl, same, opposite, anti, nores), ~ sum(., na.rm = T)) %>% 
    ungroup() %>% 
    left_join(., hits, 'method') %>% 
    mutate(n = repl + same + opposite + anti + nores,
           r_confl = anti / n * 100,
           r_repl = repl / n * 100,
           alpha = alphas[i]) %>% 
    select(method, r_confl, r_repl, r_nhit, alpha)
}

date <- format(Sys.Date(), "%d%m%y")
round2 <- function(x, decimals = 2){
  format(round(x, decimals), nsmall = decimals)
}

res_sep <- bind_rows(res_summary) %>% mutate(analysis = 'sep')


#Create figure------------------------------------------------------------------
res_both <- bind_rows(res_split, res_sep) %>% 
  select(analysis, method, r_confl, r_nhit) %>% 
  mutate(analysis = factor(analysis, levels = c('split', 'sep'),
                           labels = c('Split-data analyses',
                                      'Separate study analyses')))

dline <- tibble(analysis = as_factor(levels(res_both$analysis)),
                yint = c(.125, NA))

ggrep <- ggrepel::geom_label_repel
(pnc <- ggplot(res_both, aes(r_nhit, r_confl)) +
    geom_hline(yintercept = 0, linetype = 'solid') +
    geom_hline(data = dline, aes(yintercept = yint), linetype = '44',
               color = '#E03E18') +
    geom_point(alpha = .5) +
    ggrep(aes(label = method), 
          box.padding   = 0.7, 
          point.padding = 0.0,
          min.segment.length = unit(0, 'lines'),
          size = 2.3,
          max.overlaps = 20,
          alpha = .7,
          seed = 1,
          max.time = 2) +
    scale_y_continuous(limits = c(-1.8, 9.5),
                       labels = scales::percent_format(scale = 1)) +
    labs(x = 'Number of significant taxa (NHits)',
         y = 'Percentage of conflicting results (Conflict%)',
         title = 'Figure A13: Sensitivity and percentage of conflicting results') +
    facet_grid(. ~ analysis, scales = 'free_x') +
    theme_light() +
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 8),
          strip.text = element_text(size = 10),
          plot.title = element_text(face = 'bold', size = 12),
          plot.title.position = 'plot')
)

ggsave(filename = paste0('fig_A13_sens_confl_', date,'.png'),
       plot = pnc, width = 170, height = 100, dpi = 300, units = 'mm',
       bg = 'white')

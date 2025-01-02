library(tidyverse)
library(data.table)

#Load files and prepare tibble of results---------------------------------------
load('res_main.rds')
load('res_ancombc_251024.rds')
load('data_meta_041023.rds')
load('data_prevl_061023.rds')

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
  .[, method := ifelse(method == 'LogR (Firth)', 'Logistic\nregression', method)] %>% 
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
  filter(!disease %in% c('HIV', 'schizophrenia', 'NASH')) %>% 
  filter(n >= 50)


#Summarize results--------------------------------------------------------------
alphas <- c(.01, .05, .10, .20)
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
  # stopifnot(nrow(res_combined %>% distinct(study1)) == 37)
  # stopifnot(nrow(res_combined %>% distinct(study2)) == 37)
  
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


#Create Figure------------------------------------------------------------------
date <- format(Sys.Date(), "%d%m%y")
round2 <- function(x, decimals = 2){
  format(round(x, decimals), nsmall = decimals)
}

res_figure <- bind_rows(res_summary) %>%
  group_by(alpha) %>% 
  mutate(s_confl = -as.numeric(scale(sqrt(r_confl))),
         s_repl = as.numeric(scale(r_repl)),
         s_nhit = as.numeric(scale(r_nhit))) %>% 
  ungroup() %>% 
  select(method, alpha, r_confl, r_repl, r_nhit, s_confl, s_repl, s_nhit) %>% 
  pivot_longer(cols = -c(method, alpha)) %>% 
  separate(name, into = c('type', 'var')) %>% 
  pivot_wider(names_from = type, values_from = value) %>% 
  mutate(name = factor(var,
                       levels = c('confl','repl', 'nhit'),
                       labels = c('Conflict%', 'Replication%', 'NHits'))) %>% 
  group_by(method) %>% 
  mutate(total = mean(s)) %>% 
  ungroup() %>% 
  mutate(method = fct_reorder(method, total),
         res = case_when(var == 'confl' ~ paste0(round2(r, 1), '%'),
                         var == 'repl' ~ paste0(round2(r, 0), '%'),
                         var == 'nhit' ~ round2(r, 0) %>% as.character()),
         ff = ifelse(alpha == .05, 'bold', 'plain'))


alphas2 <- rep(c('.01', '.05', '.10', '.20'), 3)
ffs <- rep(c('plain', 'bold', 'plain', 'plain'), 3)
x_lab <- expression(paste("Nominal FDR level in exploratory datasets (", alpha, ")"))

(p1 <- ggplot(res_figure, aes(interaction(alpha, name), method, fill = s)) +
    geom_tile() +
    geom_text(aes(label = res, fontface = ff),
              color = 'grey30', size = 2.2, hjust = 0.55) +
    geom_vline(xintercept = c(4.5, 8.5), color = 'white',
               linewidth = 1) +
    annotate(geom = 'text', x = 2.5, y = 15.6, vjust = 1, hjust = .5, color = 'grey0',
             label = 'Percentage of conflicting results\n(Conflict%)', size = 2.5) +
    annotate(geom = 'text', x = 6.5, y = 15.6, vjust = 1, hjust = .5, color = 'grey0',
             label = 'Replication percentage\n(Replication%)', size = 2.5) +
    annotate(geom = 'text', x = 10.5, y = 15.6, vjust = 1, hjust = .5, color = 'grey0',
             label = 'Number of significant taxa\n(NHits)', size = 2.5) +
    annotate(geom = 'text', x = 1:12, y = 0, vjust = 0, hjust = .5, color = 'grey20',
             label = alphas2, size = 2.5, fontface = ffs) +
    scale_fill_gradient2(low = "#f7AD19", mid = "white", high = "#70aed1",
                         midpoint = 0, na.value = 'white',
                         breaks = -2:2,
                         labels = c('-2 (Worse)', '-1', '0', '1', '2 (Better)'),
                         name = 'Standardized\nvalue') +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(#title = 'Figure 4: Results of the separate analyses',
         x = x_lab) +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 9, color = 'grey0'),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 8),
          legend.title = element_text(size = 9, color = 'grey0'),
          legend.text = element_text(size = 8),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0, size = 12, face = 'bold',
                                    margin = margin(0, 0, 10, 0)),
          plot.title.position = 'plot')
)

ggsave(filename = paste0('fig_A18_filter_sample_size_sep_', date,'.png'),
       plot = p1, width = 170, height = 120, dpi = 300, unit = 'mm',
       bg = 'white')

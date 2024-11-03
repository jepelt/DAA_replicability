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
  filter(!disease %in% c('HIV', 'schizophrenia', 'NASH'))


#Summarize results--------------------------------------------------------------
expl_studies <- left_join(
  res_all %>% select(type, disease, taxon, study, n1 = n0, data_id),
  res_all %>% select(type, disease, taxon, n2 = n0, is2),
  by = join_by(type, disease, taxon, n1 < n2)) %>% 
  filter(!is.na(is2)) %>% 
  distinct(study, .keep_all = T) %>% 
  pull(data_id)

res_exploratory <- res_all %>% 
    filter(data_id %in% expl_studies) %>% 
    mutate(q = case_when(method == 'ANCOM-BC2' ~ q05, T ~ q)) %>% 
    filter(q < .05) %>% 
    select(method, type, disease, taxon, study1 = study, n1 = n0, est1 = est, q)
  
res_validation <- res_all %>% 
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

hits <- res_all %>% 
  filter(data_id %in% expl_studies) %>% 
  mutate(q = ifelse(method == 'ANCOM-BC2', q05, q)) %>% 
  group_by(study, method) %>%
  summarize(r_nhit = sum(q < .05, na.rm = T)) %>% 
  rename(study1 = study)

res_bydata <- res_combined %>% 
  group_by(study1, study2, method) %>% 
  summarize_at(vars(repl, same, opposite, anti, nores),
               ~ sum(., na.rm = T)) %>% 
  ungroup() %>% 
  left_join(hits, ., c('study1', 'method')) %>% 
  mutate(n = repl + same + opposite + anti + nores,
         r_confl = anti / n * 100,
         r_repl = repl / n * 100) %>% 
  select(study1, study2, method, n_cand = n, r_confl, r_repl, r_nhit)

res_overall <- res_combined %>% 
  group_by(method, study1, taxon) %>% 
  summarize(n = n(),
            confl = sum(anti) / n,
            repl = sum(repl) / n) %>% 
  ungroup() %>% 
  group_by(method) %>% 
  summarize(confl = 100 * mean(confl),
            repl = 100 * mean(repl))

hits_overall <- res_all %>% 
  filter(data_id %in% expl_studies) %>% 
  mutate(q = ifelse(method == 'ANCOM-BC2', q05, q)) %>% 
  group_by(method) %>%
  summarize(nhit = sum(q < .05, na.rm = T) / 37)

#Create Figure------------------------------------------------------------------
date <- format(Sys.Date(), "%d%m%y")
round2 <- function(x, decimals = 2){
  format(round(x, decimals), nsmall = decimals)
}

res_figure <- res_bydata %>% 
  left_join(res_overall, 'method') %>% 
  left_join(hits_overall, 'method') %>% 
  mutate(m_confl = fct_reorder(method, -confl),
         m_repl = fct_reorder(method, repl),
         m_nhit = fct_reorder(method, nhit))

res_mean <- res_figure %>% distinct(method, .keep_all = T)

geom_qr <- ggbeeswarm::geom_quasirandom
(p1 <- ggplot(res_figure%>% filter(!is.na(r_confl)),
              aes(m_confl, r_confl, size = n_cand)) +
    geom_qr(alpha = .07, width = .2) +
    geom_point(data = res_mean, aes(y = confl), shape = 124, size = 7, color = '#d61800') +
    scale_y_continuous(limits = c(0, 52), breaks = c(0, 10, 20, 30, 40),
                       labels = paste0(c(0, 10, 20, 30, 40), '%')) +
    scale_size_continuous(breaks = c(1, 10, 50),
                          name = 'Number of\ncandidate\ntaxa') +
    coord_flip() +
    labs(#title = 'Figure A4.1: Distribution of the percentage of conflicting results\nin the split-data analyses',
      y = 'Percentage of conflicting results (Conflict%)') +
    theme_light() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(size = 9, color = 'grey0'),
          axis.text.x = element_text(size = 8),
          legend.title = element_text(size = 9, color = 'grey0'),
          legend.text = element_text(size = 8),
          plot.title = element_text(hjust = 0, size = 12, face = 'bold',
                                    margin = margin(0, 0, 10, 0)),
          plot.title.position = 'plot')
)

ggsave(filename = paste0('fig_A10.1_distr_confl_', date,'.png'),
       plot = p1, width = 170, height = 150, dpi = 300, unit = 'mm',
       bg = 'white')


(p2 <- ggplot(res_figure %>% filter(!is.na(r_repl)),
              aes(m_repl, r_repl, size = n_cand)) +
    geom_qr(alpha = .07, width = .2) +
    geom_point(data = res_mean, aes(y = repl), shape = 124, size = 7, color = '#d61800') +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),
                       labels = paste0(seq(0, 100, 20), '%')) +
    scale_size_continuous(breaks = c(1, 10, 50),
                          name = 'Number of\ncandidate\ntaxa') +
    coord_flip() +
    labs(#title = 'Figure A4.2: Distribution of the percentage of replicated results\nin the split-data analyses',
      y = 'Replication percentage (Replication%)') +
    theme_light() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(size = 9, color = 'grey0'),
          axis.text.x = element_text(size = 8),
          legend.title = element_text(size = 9, color = 'grey0'),
          legend.text = element_text(size = 8),
          plot.title = element_text(hjust = 0, size = 12, face = 'bold',
                                    margin = margin(0, 0, 10, 0)),
          plot.title.position = 'plot')
)

ggsave(filename = paste0('fig_A10.2_distr_repl_', date,'.png'),
       plot = p2, width = 170, height = 150, dpi = 300, unit = 'mm',
       bg = 'white')


(p3 <- ggplot(res_figure, aes(m_nhit, r_nhit)) +
    geom_qr(alpha = .07, width = .2) +
    geom_point(data = res_mean, aes(y = nhit), shape = 124, size = 7, color = '#d61800') +
    labs(#title = 'Figure A4.3: Distribution of the number of significant taxa detected\nin the split-data analyses',
      y = 'Number of significant taxa (NHits)') +
    coord_flip() +
    theme_light() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(size = 9, color = 'grey0'),
          axis.text.x = element_text(size = 8),
          legend.title = element_text(size = 9, color = 'grey0'),
          legend.text = element_text(size = 8),
          plot.title = element_text(hjust = 0, size = 12, face = 'bold',
                                    margin = margin(0, 0, 10, 0)),
          plot.title.position = 'plot')
)

ggsave(filename = paste0('fig_A10.3_distr_nhit_', date,'.png'),
       plot = p3, width = 170, height = 150, dpi = 300, unit = 'mm',
       bg = 'white')





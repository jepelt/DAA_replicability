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
res_exploratory <- res_all %>% 
  filter(analysis == 'half' & half == 1) %>% 
  mutate(q = ifelse(method == 'ANCOM-BC2', q05, q)) %>% 
  filter(q < .05) %>% 
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

hits <- res_all %>% 
  filter(analysis == 'half' & half == 1) %>% 
  mutate(q = ifelse(method == 'ANCOM-BC2', q05, q)) %>% 
  group_by(study, iter, method) %>%
  summarize(r_nhit = sum(q < .05, na.rm = T))

res_bydata <- res_combined %>% 
  group_by(study, iter, method) %>% 
  summarize_at(vars(repl, same, opposite, anti, nores),
               ~ sum(., na.rm = T)) %>% 
  ungroup() %>% 
  left_join(hits, ., c('study', 'iter', 'method')) %>% 
  mutate(n = repl + same + opposite + anti + nores,
         r_confl = anti / n * 100,
         r_repl = repl / n * 100) %>% 
  select(study, iter, method, n_cand = n, r_confl, r_repl, r_nhit) %>% 
  group_by(method) %>% 
  mutate(confl = sum(r_confl * n_cand, na.rm = T) / sum(n_cand, na.rm = T),
         repl = sum(r_repl * n_cand, na.rm = T) / sum(n_cand, na.rm = T),
         nhit = sum(r_nhit, na.rm = T) / 285) %>% 
  ungroup() %>% 
  mutate(m_confl = fct_reorder(method, -confl),
         m_repl = fct_reorder(method, repl),
         m_nhit = fct_reorder(method, nhit))


#Create Figure------------------------------------------------------------------
date <- format(Sys.Date(), "%d%m%y")
round2 <- function(x, decimals = 2){
  format(round(x, decimals), nsmall = decimals)
}

res_mean <- res_bydata %>% distinct(method, .keep_all = T)

geom_qr <- ggbeeswarm::geom_quasirandom
(p1 <- ggplot(res_bydata %>% filter(!is.na(r_confl)),
              aes(m_confl, r_confl, size = n_cand)) +
    geom_qr(alpha = .07, width = .2) +
    geom_point(data = res_mean, aes(y = confl), shape = 124, size = 7, color = '#d61800') +
    scale_y_continuous(limits = c(0, 40), breaks = c(0, 10, 20, 30, 40),
                       labels = paste0(c(0, 10, 20, 30, 40), '%')) +
    scale_size_continuous(breaks = c(1, 10, 100),
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

ggsave(filename = paste0('fig_A4.1_distr_confl_', date,'.png'),
       plot = p1, width = 170, height = 150, dpi = 300, unit = 'mm',
       bg = 'white')


(p2 <- ggplot(res_bydata %>% filter(!is.na(r_repl)),
              aes(m_repl, r_repl, size = n_cand)) +
    geom_qr(alpha = .07, width = .2) +
    geom_point(data = res_mean, aes(y = repl), shape = 124, size = 7, color = '#d61800') +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),
                       labels = paste0(seq(0, 100, 20), '%')) +
    scale_size_continuous(breaks = c(1, 10, 100),
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

ggsave(filename = paste0('fig_A4.2_distr_repl_', date,'.png'),
       plot = p2, width = 170, height = 150, dpi = 300, unit = 'mm',
       bg = 'white')


(p3 <- ggplot(res_bydata, aes(m_nhit, r_nhit)) +
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

ggsave(filename = paste0('fig_A4.3_distr_nhit_', date,'.png'),
       plot = p3, width = 170, height = 150, dpi = 300, unit = 'mm',
       bg = 'white')




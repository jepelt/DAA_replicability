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
  .[, method := tstrsplit(method, " ", fixed = T, keep = 1)] %>% 
  .[, analysis := tstrsplit(data_id, " ", fixed = T, keep = 1)] %>% 
  .[, q := p.adjust(p, method = 'BH'), by = c('data_id', 'method')] %>% 
  .[, q := ifelse(grepl('DESeq2|LDM|Zico', method, perl = T), p_fdr, q)] %>%
  .[, t91.7 := qt(.834 + (1 - .834) / 2, df = df)] %>% 
  .[, lwr := ifelse(is.na(lwr), est - t91.7 * se, lwr)] %>% 
  .[, upr := ifelse(is.na(upr), est + t91.7 * se, upr)] %>% 
  .[, disease := ifelse(disease == 'EDD', 'CDI', disease)] %>% 
  filter(site == 'stool') %>% 
  filter(!(method == 'LogR' & prevl_control == 1 & prevl_case == 1))


#Summarize results--------------------------------------------------------------
res_exploratory <- res_all %>% 
  filter(analysis == 'half' & half == 1) %>% 
  select(method, type, study, iter, taxon, est1 = est)

res_validation <- res_all %>% 
  filter(analysis == 'half' & half == 2) %>% 
  select(method, study, iter, taxon, est2 = est)

res_combined <- inner_join(res_exploratory, res_validation, 
                           by = c('method', 'study', 'iter', 'taxon'))

#Create Figure------------------------------------------------------------------
res_figure <- res_combined %>% 
  group_by(type, study, method, iter) %>% 
  summarize(r = cor(est1, est2, method = 'spearman', use = 'pairwise.complete'),
            n = n()) %>% 
  group_by(method) %>% 
  mutate(mean_r = tanh(mean(atanh(r), na.rm = T))) %>% 
  ungroup() %>% 
  mutate(method = fct_reorder(method, mean_r))

res_mean <- res_figure %>% distinct(method, .keep_all = T)

geom_qr <- ggbeeswarm::geom_quasirandom
(p1 <- ggplot(res_figure, aes(method, r)) + 
    geom_hline(yintercept = 0, linetype = '44', linewidth = 1, color = 'grey40') +
    geom_qr(size = 2, alpha = .07, width = .2) +
    geom_point(data = res_mean, aes(method, mean_r), shape = 124, size = 7,
               color = '#d61800') +
    #scale_x_discrete(limits = rev) +
    coord_flip() +
    labs(#title = 'Figure A5: Correlation of DA estimates between exploratory and validation datasets\nin the split-data analyses (285 pairs of datasets)',
         y = 'Spearman correlation of DA estimates') +
    theme_light() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(size = 8, color = 'grey0'),
          axis.text.x = element_text(size = 8),
          legend.title = element_text(size = 9, color = 'grey0'),
          legend.text = element_text(size = 8),
          plot.title = element_text(hjust = 0, size = 12, face = 'bold',
                                    margin = margin(0, 0, 10, 0)),
          plot.title.position = 'plot')
)

date <- format(Sys.Date(), "%d%m%y")
ggsave(filename = paste0('fig_A5_cor_split_', date,'.png'),
       plot = p1, width = 170, height = 150, dpi = 300, unit = 'mm',
       bg = 'white')

#writexl::write_xlsx(res_mean, 'table_est_corr_split_201024.xlsx')

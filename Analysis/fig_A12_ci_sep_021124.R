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
  filter(!method %in% c('edgeR', 'LDM', 'ZicoSeq')) %>% 
  filter(!(method == 'ANCOM-BC2' & !ss_pass_p)) %>% 
  mutate(n0 = n + avg_ls / 10 ^ 10,
         is2 = T) %>% 
  filter(analysis == 'whole' & n_case >= 10 & n_control >= 10) %>% 
  filter(!disease %in% c('HIV', 'schizophrenia', 'NASH'))


#Calculate the metrics for overlap of CIs---------------------------------------
expl_studies <- left_join(
  res_all %>% select(type, disease, taxon, study, n1 = n0, data_id),
  res_all %>% select(type, disease, taxon, n2 = n0, is2),
  by = join_by(type, disease, taxon, n1 < n2)) %>% 
  filter(!is.na(is2)) %>% 
  distinct(study, .keep_all = T) %>% 
  pull(data_id)

res_exploratory <- res_all %>% 
  filter(data_id %in% expl_studies) %>% 
  mutate(q = ifelse(method == 'ANCOM-BC2', q05, q)) %>% 
  filter(q < .05) %>% 
  filter(prevl_control >= .10 & prevl_case >= .10) %>% 
  select(method, type, disease, taxon, study1 = study, n1 = n0,
         q, lwr1 = lwr, upr1 = upr)
  
res_validation <- res_all %>% 
    filter(prevl_control >= .10 & prevl_case >= .10) %>% 
    select(method, type, disease, taxon, study2 = study, n2 = n0, is2,
           lwr2 = lwr, upr2 = upr)
  
res_combined <- left_join(res_exploratory, res_validation,
                          by = join_by(method, type, disease, taxon, n1 < n2)) %>% 
  filter(!is.na(lwr1) & !is.na(upr1) & !is.na(lwr2) & !is.na(upr2)) %>% 
  mutate(ol_start = pmax(lwr1, lwr2),
         ol_end = pmin(upr1, upr2),
         overlap = ol_start <= ol_end)


#Create figure------------------------------------------------------------------
res_mean <- res_combined %>% 
  group_by(method, study1, taxon) %>% 
  summarize(ol = mean(overlap)) %>% 
  group_by(method) %>% 
  summarize(m = mean(ol)) %>% 
  ungroup() %>% 
  mutate(method = fct_reorder(method, m))
  
res_figure <- res_combined %>% 
  group_by(method, study1, study2) %>% 
  summarize(ol = mean(overlap),
            n = n()) %>% 
  mutate(method = factor(method, levels(res_mean$method)))
  
geom_qr <- ggbeeswarm::geom_quasirandom
(p <- ggplot(res_figure, aes(method, ol, size = n)) +
    geom_hline(yintercept = 0:1, color = 'grey', linewidth = .5) +
    #geom_hline(yintercept = .95, linetype = '33', color = 'grey', linewidth = .5) +
    geom_qr(alpha = .07, width = .2) +
    geom_point(data = res_mean, aes(method, m), shape = 124, size = 7, color = '#d61800') +
    scale_y_continuous(breaks = seq(0, 100, 10) / 100,
                       labels = paste0(seq(0, 100, 10), '%'),
                       limits = c(0, 1.00)) +
    scale_size_continuous(breaks = c(1, 10, 40),
                          name = 'Number of\ncandidate\ntaxa') +
    coord_flip() +
    labs(#title = 'Figure A6: Overlap percentage of 83.4% confidence intervals in the split-data analyses',
      y = 'Overlap percentage of confidence intervals (CI%)') +
    theme_light() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(size = 8, color = 'grey0'),
          axis.text.x = element_text(size = 8),
          legend.title = element_text(size = 9, color = 'grey0'),
          legend.text = element_text(size = 8),
          plot.title = element_text(hjust = 0, size = 11, face = 'bold',
                                    margin = margin(0, 0, 10, 0)),
          plot.title.position = 'plot')
)

date <- format(Sys.Date(), "%d%m%y")
ggsave(filename = paste0('fig_A12_ci_split_', date,'.png'),
       plot = p, width = 170, height = 150, dpi = 300, unit = 'mm',
       bg = 'white')


res_table <- res_mean %>% mutate(ol = round(100 * m, 0))

#writexl::write_xlsx(res_table, 'table_ci_sep_301024.xlsx')

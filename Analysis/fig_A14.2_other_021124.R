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
  .[method %in% c(#'ALDEx2 (glm)',
    'ALDEx2 (t-test0.5)',
    'ALDEx2 (Wilcoxon)',
    #'ANCOM-BC2 (pc = 0)',
    #'corncob (LRT, ev = TRUE)',
    'corncob (LRT, ev = FALSE)',
    #'DESeq2 (Wald, Default)',
    'DESeq2 (LRT, Default)',
    #'edgeR (QL, TMM)',
    #'fastANCOM (0.5, 0.05, 0.05)',
    #'LDM (CLR = FALSE)',
    'LDM (CLR = TRUE)',
    #'limma-voom (TMM)',
    #'LinDA (Adaptive, 0.5)',
    #'LogR (Firth)',
    #'MaAsLin2 (LOG, TSS)',
    'MaAsLin2 (AST, TSS)',
    'MaAsLin2 (LOG, CSS)',
    'MaAsLin2 (LOG, TMM)',
    'MaAsLin2 (NONE, CLR)',
    #'mSeq (CSS)',
    'NB',
    #'ORMW (Score, TSS)',
    'ORMW (Score, GMPR)',
    'ORMW (Score, Wrench)',
    'radEmu (Wald)'
    #'ZicoSeq (Default)'
  )] %>% 
  .[, method := ifelse(method == 'ALDEx2 (Wilcoxon)', 'ALDEx2-Wilcox', method)] %>% 
  .[, method := ifelse(method == 'ALDEx2 (t-test0.5)', 'ALDEx2-scale', method)] %>% 
  .[, method := ifelse(method == 'corncob (LRT, ev = FALSE)', 'corncob-UEV', method)] %>% 
  .[, method := ifelse(method == 'DESeq2 (LRT, Default)', 'DESeq2-LRT', method)] %>% 
  .[, method := ifelse(method == 'LDM (CLR = TRUE)', 'LDM-CLR', method)] %>% 
  .[, method := ifelse(method == 'MaAsLin2 (AST, TSS)', 'MaAsLin2-AST', method)] %>% 
  .[, method := ifelse(method == 'MaAsLin2 (LOG, CSS)', 'MaAsLin2-CSS', method)] %>% 
  .[, method := ifelse(method == 'MaAsLin2 (LOG, TMM)', 'MaAsLin2-TMM', method)] %>% 
  .[, method := ifelse(method == 'MaAsLin2 (NONE, CLR)', 'MaAsLin2-CLR', method)] %>% 
  .[, method := ifelse(method == 'NB', 'NegBin', method)] %>% 
  .[, method := ifelse(method == 'ORMW (Score, TSS)', 'ORM/W', method)] %>% 
  .[, method := ifelse(method == 'ORMW (Score, GMPR)', 'ORM/W-GMPR', method)] %>% 
  .[, method := ifelse(method == 'ORMW (Score, Wrench)', 'ORM/W-Wrench', method)] %>% 
  .[, method := ifelse(method == 'mSeq (CSS)', 'metagenomeSeq', method)] %>%  
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

datasets <- res_all %>% 
  distinct(study, .keep_all = T) %>% 
  select(study, type, disease, n, mean_count)

res_bydata <- res_combined %>% 
  group_by(method, study, taxon) %>% 
  summarize_at(vars(repl, anti), ~ sum(., na.rm = T)) %>% 
  ungroup() %>% 
  group_by(method, study) %>% 
  summarize_at(vars(repl, anti), ~ sum(. > 0.5)) %>% 
  left_join(res_exploratory %>% distinct(study), ., 'study') %>% 
  mutate(method = ifelse(is.na(method), 'ALDEx2', method)) %>% 
  pivot_wider(names_from = method, values_from = c(anti, repl)) %>% 
  pivot_longer(cols = -study) %>% 
  left_join(., datasets, 'study') %>% 
  mutate(value = ifelse(is.na(value), 0, value)) %>% 
  separate(name , into = c('var', 'method'), sep = '_') %>% 
  mutate(last = 0)


#Create figure------------------------------------------------------------------
load('res_beta_141024.rds')
round2 <- function(x, decimals = 2){
  format(round(x, decimals), nsmall = decimals)
}
std01 <- function(x){(x - min(x)) / (max(x) - min(x))}
date <- format(Sys.Date(), "%d%m%y")

res_sum <- res_bydata %>% 
  group_by(method, var) %>% 
  summarize(value = sum(value)) %>% 
  mutate(last = 1, type = '', n = NULL, disease = '', study = 'tot')

res_figure <- bind_rows(res_bydata, res_sum) %>% 
  group_by(study, var) %>% 
  mutate(std = std01(value),
         std = ifelse(is.na(std), 0, std),
         n = ceiling(n / 2),
         disease = ifelse(disease == 'adenoma', 'Ade.', disease),
         disease = ifelse(disease == 'cirrhosis', 'Cirr.', disease),
         disease = ifelse(disease == 'CIRR', 'Cirr.', disease),
         disease = ifelse(disease == 'parkinson', 'PD', disease),
         disease = ifelse(disease == 'diarrhea', 'Diarr.', disease),
         disease = ifelse(disease == 'parkinson', 'PD', disease),
         disease = ifelse(disease == 'hypertension', 'HT', disease),
         disease = ifelse(disease == 'pre-hypertension', 'PHT', disease),
         disease = ifelse(disease == 'asthma', 'Asth.', disease),
         disease = ifelse(disease == 'schizophrenia', 'Schiz.', disease),
         disease = ifelse(disease == 'migraine', 'Migr.', disease),
         disease = ifelse(disease == 'cephalosporins', 'Ceph.', disease),
         disease = ifelse(disease == 'ME/CFS', 'ME', disease),
         disease = ifelse(disease == 'EDD', 'CDI', disease),
         type = ifelse(type == 'Shotgun', 'SG', type),
         xc = tolower(paste(method, var))) %>% 
  ungroup() %>% 
  arrange(-last, desc(type), desc(tolower(disease)), n) %>% 
  mutate(study = as_factor(study)) %>% 
  arrange(xc) %>% 
  mutate(xn = dense_rank(xc))

types <- res_figure %>% group_by(study) %>% slice(1) %>% pull(type)
diseases <- res_figure %>% group_by(study) %>% slice(1) %>% pull(disease)
ns <- res_figure %>% group_by(study) %>% slice(1) %>% pull(n)

beta_divs <- betad %>% 
  filter(str_detect(data_id, 'half')) %>% 
  mutate(study = str_replace(data_id, 'half 16S |half sg ', ''),
         study = str_sub(study, 1, -5)) %>% 
  group_by(study) %>% 
  summarize(r2 = mean(r2)) %>% 
  mutate(r2 = str_sub(round2(r2, 2), 2, 4)) %>% 
  left_join(res_figure %>% distinct(study), ., 'study') %>% 
  pull(r2)
  
methods <- res_figure %>% distinct(method) %>% arrange(tolower(method)) %>% 
  pull(method)
vlines <- tibble(x = .5 + seq(0, 26, 2), xend = x, y = 0.5, yend = 58.5)

ann_size <- 2.2
(p1 <- ggplot(res_figure, aes(xn, study)) +
    geom_tile(aes(alpha = std, fill = var)) +
    geom_text(aes(label = value), size = 2.1, vjust = .4) +
    geom_segment(data = vlines, aes(x = x, xend = xend, y = y, yend = yend),
                 linewidth = .4, color = 'black') +
    geom_hline(yintercept = c(1.55), linewidth = .4) +
    geom_hline(yintercept = c(34.55), linewidth = .4, linetype = '33') +
    annotate(geom = 'point', x = 10, y = 60.5, color = 'white') +
    annotate(geom = 'text', x = -2.5 - 2, y = 1:58, vjust = .4, hjust = .5, 
             color = 'grey20', label = types, size = ann_size) +
    annotate(geom = 'text', x = -1.7 - 2, y = 1:58, vjust = .4, hjust = 0,
             color = 'grey20', label = diseases, size = ann_size) +
    annotate(geom = 'text', x = 1.0 - 2, y = 1:58, vjust = .4, hjust = 1,
             color = 'grey20', label = beta_divs, size = ann_size) +
    annotate(geom = 'text', x = 2.3 - 2, y = 1:58, vjust = .4, hjust = 1,
             color = 'grey20', label = ns, size = ann_size) +
    annotate(geom = 'text', x = c(-3.0, -1.8, 0.0, 1.7) - 2, y = 59, vjust = 0,
             hjust = 0, color = 'grey20', label = c('Seq', 'Cond', 'Beta', 'N'),
             size = ann_size, fontface = 'bold', angle = 0) +
    annotate(geom = 'text', x = -3.0, y = 1, vjust = 0.4, hjust = 0.5,
             color = 'grey20', label = 'Total', size = ann_size, fontface = 'bold',
             angle = 0) +
    scale_x_continuous(breaks = seq(1.5, 2 * 13 - .5, 2), labels = methods,
                       expand = c(0, 0), limits = c(-5.1, .5 + 2 * 13)) +
    scale_alpha(breaks = seq(0, 1, .25), labels = c('Min', '', '', '', 'Max'),
                name = 'Value\nrelative\nto other\nmethods',
                guide = guide_legend(reverse = T)) +
    scale_fill_manual(name = 'Result', labels = c('Conflicting', 'Replicated'),
                      values = c("#f7AD19", "#70aed1")) +
    labs(#title = 'Figure 3: Number of conflicting and replicated results on split datasets',
         y = 'Split dataset') +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                     size = 8),
          axis.title.y = element_text(size = 10),
          axis.text.y = element_blank(),
          legend.key.size = unit(.4, 'cm'),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          panel.grid = element_blank(),
          plot.title = element_text(face = 'bold', size = 12),
          plot.title.position = 'plot')
)

ggsave(filename = paste0('fig_A14.2_other_', date,'.png'),
       plot = p1, width = 170, height = 180, dpi = 300, unit = 'mm',
       bg = 'white')

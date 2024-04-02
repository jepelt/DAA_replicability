library(tidyverse)
library(data.table)

#Helper functions---------------------------------------------------------------
geom_qr <- ggbeeswarm::geom_quasirandom
ggarr <- ggpubr::ggarrange
pblank <- ggplot() + theme_void()
std01 <- function(x){(x - min(x)) / (max(x) - min(x))}
round2 <- function(x, decimals = 2){
  format(round(x, decimals), nsmall = decimals)
}
date <- format(Sys.Date(), "%d%m%y")


#Load files and prepare tibble of results---------------------------------------
load('res_main.rds')
load('data_meta_041023.rds')
load('data_prevl_061023.rds')

res_all <- bind_rows(res_list) %>% 
  left_join(., d_meta, by = 'data_id') %>% 
  left_join(., d_prevl, by = c('data_id', 'taxon')) %>% 
  as.data.table() %>% 
  .[method %in% c('ALDEx2 (glm)',
                  'ANCOM-BC2 (pc = 0)',
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
                  'ORMW (Score, TSS)')] %>% 
  .[, method := ifelse(method == 'mSeq (CSS)', 'metagenomeSeq', method)] %>% 
  .[, method := ifelse(method == 'ORMW (Score, TSS)', 'ORM/Wilcox', method)] %>% 
  .[, method := tstrsplit(method, " ", fixed = T, keep = 1)] %>% 
  .[, analysis := tstrsplit(data_id, " ", fixed = T, keep = 1)] %>% 
  .[, q := p.adjust(p, method = 'BH'),
    by = c('data_id', 'method')] %>% 
  .[, q := ifelse(grepl('DESeq2|LDM', method, perl = T), 
                  p_fdr, q)] %>%
  .[, t91.7 := qt(.834 + (1 - .834) / 2, df = df)] %>% 
  .[, lwr := ifelse(is.na(lwr), est - t91.7 * se, lwr)] %>% 
  .[, upr := ifelse(is.na(upr), est + t91.7 * se, upr)] %>% 
  .[, disease := ifelse(disease == 'EDD', 'CDI', disease)] %>% 
  filter(site == 'stool')


#Calculate number of conflicting and replicated results-------------------------
sgn_level1 <- sgn_level2 <- .05

ol1h <- res_all %>% 
  filter(analysis == 'half' & half == 1) %>% 
  select(method, type, study, iter, taxon,
         pco1 = prevl_control, pca1 = prevl_case,
         est1 = est, q, lwr1 = lwr, upr1 = upr)

ol2h <- res_all %>% 
  filter(analysis == 'half' & half == 2) %>% 
  mutate(is2 = T) %>% 
  select(method, type, study, iter, taxon, is2,
         pco2 = prevl_control, pca2 = prevl_case,
         est2 = est, p, lwr2 = lwr, upr2 = upr)

olh <- left_join(ol1h, ol2h, by = c('method', 'study', 'iter', 'taxon')) %>% 
  filter(q < sgn_level1 & !is.na(is2)) %>% 
  mutate(sgn1 = q < sgn_level1,
         d1 = ifelse(is.na(est1) | is.nan(est1), 0, sign(round(est1, 8))),
         s1 = if_else(sgn1 & d1 != 0, T, F, F),
         sgn2 = p < sgn_level2,
         d2 = ifelse(is.na(est2) | is.nan(est2), 0, sign(round(est2, 8))),
         s2 = if_else(sgn2 & d2 != 0, T, F, F),
         nores = d2 == 0,
         repl = s1 & s2 & d1 == d2,
         same = s1 & !s2 & d1 == d2 & !nores,
         opposite = s1 & !s2 & d1 != d2 & !nores,
         anti = s1 & s2 & d1 != d2,
         is_ci = pco1 >= .10 & pca1 >= .10 & pco1 >= .10 & pca1 >= .10,
         cvg = lwr1 >= lwr2 & lwr1 <= upr2 | lwr2 >= lwr1 & lwr2 <= upr1,
         cvg = ifelse(is.na(cvg), F, cvg),
         cvg = ifelse(is_ci, cvg, NA)) %>%
  mutate_if(is.logical, as.numeric) %>% 
  mutate(check = repl + same + opposite + anti + nores) %>%
  select(method, study, iter, taxon,
         est1, q, est2, p, sgn1, d1, s1, sgn2, d2, s2,
         repl, same, opposite, anti, nores, is_ci, cvg, check)

stopifnot(sum(olh$check != 1) == 0)


datasets <- res_all %>% 
  distinct(study, .keep_all = T) %>% 
  select(study, type, disease, n, mean_count)

dtb_h <- olh %>% 
  group_by(method, study, taxon) %>% 
  summarize_at(vars(repl, anti), ~ sum(., na.rm = T)) %>% 
  ungroup() %>% 
  group_by(method, study) %>% 
  summarize_at(vars(repl, anti), ~ sum(. > 0.5)) %>% 
  left_join(ol1h %>% distinct(study), ., 'study') %>% 
  mutate(method = ifelse(is.na(method), 'ALDEx2', method)) %>% 
  pivot_wider(names_from = method, values_from = c(anti, repl)) %>% 
  pivot_longer(cols = -study) %>% 
  left_join(., datasets, 'study') %>% 
  mutate(value = ifelse(is.na(value), 0, value)) %>% 
  separate(name , into = c('var', 'method'), sep = '_') %>% 
  mutate(last = 0)

dsh <- dtb_h %>% 
  group_by(method, var) %>% 
  summarize(value = sum(value)) %>% 
  mutate(last = 1, type = '', n = NULL, disease = '', study = 'tot')

dtb_h2 <- bind_rows(dtb_h, dsh) %>% 
  group_by(study, var) %>% 
  mutate(std = std01(value),
         std = ifelse(is.na(std), 0, std),
         n = ceiling(n / 2),
         disease = ifelse(disease == 'adenoma', 'Ade.', disease),
         disease = ifelse(disease == 'NASH', 'NAS', disease),
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


#Create Figure 3----------------------------------------------------------------
ns <- dtb_h2 %>% group_by(study) %>% slice(1) %>% pull(n)
diseases <- dtb_h2 %>% group_by(study) %>% slice(1) %>% pull(disease)
typs <- dtb_h2 %>% group_by(study) %>% slice(1) %>% pull(type)
methods <- dtb_h2 %>% distinct(method) %>% arrange(tolower(method)) %>% 
  pull(method)
vlines <- tibble(x = .5 + seq(0, 26, 2), xend = x, y = 0.5, yend = 58.5)

(pdh <- ggplot(dtb_h2, aes(xn, study)) +
    geom_tile(aes(alpha = std, fill = var)) +
    geom_text(aes(label = value), size = 2.5, vjust = .4) +
    geom_segment(data = vlines, aes(x = x, xend = xend, y = y, yend = yend),
                 linewidth = .4, color = 'black') +
    geom_hline(yintercept = c(1.55), linewidth = .4) +
    geom_hline(yintercept = c(34.55), linewidth = .4, linetype = '33') +
    annotate(geom = 'point', x = 10, y = 60.5, color = 'white') +
    annotate(geom = 'text', x = -3.5, y = 1:58, vjust = .4, hjust = .5, 
             color = 'grey20', label = typs, size = 2.5) +
    annotate(geom = 'text', x = -2.4, y = 1:58, vjust = .4, hjust = 0,
             color = 'grey20', label = diseases, size = 2.5) +
    annotate(geom = 'text', x = .3, y = 1:58, vjust = .4, hjust = 1,
             color = 'grey20', label = ns, size = 2.5) +
    annotate(geom = 'text', x = c(-4.1, -.2, -2.5), y = 59, vjust = 0,
             hjust = 0, color = 'grey20', label = c('Seq.', 'N', 'Cond.'),
             size = 2.5, fontface = 'bold', angle = 0) +
    annotate(geom = 'text', x = -2.0, y = 1, vjust = 0.4, hjust = 0.5,
             color = 'grey20', label = 'Total', size = 2.5, fontface = 'bold',
             angle = 0) +
    scale_x_continuous(breaks = seq(1.5, 25.5, 2), labels = methods,
                       expand = c(0, 0), limits = c(-4.1, 26.5)) +
    scale_alpha(breaks = seq(0, 1, .25), labels = c('Min', '', '', '', 'Max'),
                name = 'Value\nrelative\nto other\nmethods',
                guide = guide_legend(reverse = T)) +
    scale_fill_manual(name = 'Result', labels = c('Conflicting', 'Replicated'),
                      values = c("#f7AD19", "#70aed1")) +
    labs(y = 'Original dataset',
         title = 'Figure 3: Number of conflicting and replicated results on split datasets') +
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

ggsave(filename = paste0('fig_3_', date,'.png'),
       plot = pdh, width = 170, height = 180, dpi = 300, unit = 'mm',
       bg = 'white')

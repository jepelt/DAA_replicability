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
                  'ALDEx2 (Wilcoxon)',
                  #'ANCOM-BC2 (pc = 0)',
                  'corncob (LRT, ev = TRUE)',
                  'corncob (LRT, ev = FALSE)',
                  'DESeq2 (Wald, Default)',
                  'DESeq2 (LRT, Default)',
                  'edgeR (QL, TMM)',
                  'fastANCOM (0.5, 0.05, 0.05)',
                  'LDM (CLR = FALSE)',
                  'LDM (CLR = TRUE)',
                  'limma-voom (TMM)',
                  'LinDA (Adaptive, 0.5)',
                  'LogR (Firth)',
                  'MaAsLin2 (LOG, TSS)',
                  'MaAsLin2 (AST, TSS)',
                  'MaAsLin2 (LOG, CSS)',
                  'MaAsLin2 (LOG, TMM)',
                  'MaAsLin2 (NONE, CLR)',
                  'mSeq (CSS)',
                  'NB',
                  'ORMW (Score, TSS)',
                  'ORMW (Score, GMPR)',
                  'ORMW (Score, Wrench)')] %>% 
  .[, method := ifelse(method == 'ALDEx2 (Wilcoxon)', 'ALDEx2-Wilcox', method)] %>% 
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
  .[, q := p.adjust(p, method = 'BH'),
    by = c('data_id', 'method')] %>% 
  .[, q := ifelse(grepl('DESeq2|LDM', method, perl = T), 
                  p_fdr, q)] %>%
  .[, t91.7 := qt(.834 + (1 - .834) / 2, df = df)] %>% 
  .[, lwr := ifelse(is.na(lwr), est - t91.7 * se, lwr)] %>% 
  .[, upr := ifelse(is.na(upr), est + t91.7 * se, upr)] %>% 
  .[, disease := ifelse(disease == 'EDD', 'CDI', disease)] %>% 
  filter(site == 'stool') %>% 
  mutate(n0 = n + avg_ls / 10 ^ 10,
         is2 = T) %>% 
  filter(analysis == 'whole' & n_case >= 10 & n_control >= 10) %>% 
  filter(!disease %in% c('HIV', 'schizophrenia', 'NASH'))


#Cluster methods based on significant taxa (53 datasets)------------------------
dz <- res_all %>% 
  filter(analysis == 'whole' & n_case >= 10 & n_control >= 10 & 
           disease != 'HIV') %>% 
  
  filter(!str_detect(method, 'AST|LDM-CLR|UEV|LRT|Wilcox|Pseudo|edgeR|DESeq2|NegBin|metag')) %>% 
  filter(prevl_control >= .10 & prevl_case >= .10) %>% 
  mutate(sgn = as.numeric(q < .05)) %>% 
  select(method, taxon, study, sgn) %>% 
  pivot_wider(names_from = method, values_from = sgn) %>% 
  select(-c(study, taxon)) %>% 
  as.matrix() %>% 
  t()

dists <- dist(dz, method = 'binary')  
mds <- cmdscale(dists, eig = TRUE, k = 5)

dmds <- tibble(x = mds$points[, 1],
               y = mds$points[, 2],
               method = rownames(mds$points))
eigs <- round2(100 * mds$eig / sum(mds$eig), 1)


ggrep <- ggrepel::geom_label_repel
(pclust <- ggplot(dmds, aes(x, y)) +
  geom_point() +
  ggrep(aes(label = method), 
        box.padding   = 0.5, 
        point.padding = 0.5,
        min.segment.length = unit(0, 'lines'),
        size = 3,
        max.overlaps = 20,
        alpha = .7,
        seed = 1,
        max.time = 2) +
  labs(title = 'a) 50 datasets') +
  labs(x = paste0('MDS 1 (', eigs[1], '%)'),
       y = paste0('MDS 2 (', eigs[2], '%)')) +
  theme_light() +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 8),
        plot.title = element_text(face = 'bold', size = 12),
        plot.title.position = 'plot'))


#Cluster methods based on significant taxa (CDI datasets)-----------------------
dz2 <- res_all %>% 
  filter(analysis == 'whole' & n_case >= 10 & n_control >= 10 & 
           disease != 'HIV') %>% 
  
  filter(!str_detect(method, 'AST|LDM-CLR|UEV|LRT|Wilcox|Pseudo|edgeR|DESeq2|NegBin|metag')) %>% 
  filter(str_detect(disease, 'CDI')) %>% 
  filter(prevl_control >= .10 & prevl_case >= .10) %>% 
  mutate(sgn = as.numeric(q < .05)) %>% 
  select(method, taxon, study, sgn) %>% 
  pivot_wider(names_from = method, values_from = sgn) %>% 
  select(-c(study, taxon)) %>% 
  as.matrix() %>% 
  t()

dists2 <- dist(dz2, method = 'binary')  
mds2 <- cmdscale(dists2, eig = TRUE, k = 5)

dmds2 <- tibble(x = mds2$points[, 1],
               y = mds2$points[, 2],
               method = rownames(mds2$points))
eigs <- round2(100 * mds2$eig / sum(mds2$eig), 1)


ggrep <- ggrepel::geom_label_repel
(pclust2 <- ggplot(dmds2, aes(x, y)) +
    geom_point() +
    ggrep(aes(label = method), 
          box.padding   = 0.5, 
          point.padding = 0.5,
          min.segment.length = unit(0, 'lines'),
          size = 3,
          max.overlaps = 20,
          alpha = .7,
          seed = 1,
          max.time = 2) +
    labs(title = 'b) Three CDI datasets') +
    labs(x = paste0('MDS 1 (', eigs[1], '%)'),
         y = paste0('MDS 2 (', eigs[2], '%)')) +
    theme_light(base_size = 10) +
    theme(plot.title = element_text(face = 'bold'),
          plot.title.position = 'plot'))

(p <- ggpubr::ggarrange(pclust, pclust2, nrow = 2))

ggsave(filename = paste0('figure_A16_', date,'.png'),
       plot = p, width = 170, height = 150, dpi = 300, units = 'mm',
       bg = 'white')


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
load('res_covariates.rds')
load('data_meta_041023.rds')
load('data_prevl_061023.rds')

t <- res_all %>% distinct(method)

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
                  'MaAsLin2 (LOG)',
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


#Calculate values of the evaluation metrics-------------------------------------
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


n_studies <- nrow(olh %>% distinct(study, iter))

hits_h <- ol1h %>% 
  filter(q < sgn_level1) %>% 
  group_by(method) %>% 
  summarize(nhit = n())

res_h <- olh %>% 
  group_by(method) %>% 
  summarize_at(vars(repl, same, opposite, anti, nores, cvg, is_ci),
               ~ sum(., na.rm = T)) %>% 
  ungroup() %>% 
  left_join(., hits_h, 'method') %>% 
  mutate(n = repl + same + opposite + anti + nores,
         panti = anti / n * 100,
         prepl = repl / n * 100,
         popp = (anti + opposite + .5 * nores) / n * 100,
         pcvg = cvg / is_ci * 100,
         method_an = fct_reorder(method, -(anti +.001) / n),
         method_re = fct_reorder(method, prepl),
         method_op = fct_reorder(method, -popp),
         method_ci = fct_reorder(method, pcvg),
         method_nh = fct_reorder(method, nhit))


#Create Figure A5---------------------------------------------------------------
title_size <- 12
axis_title_size <- 10
axis_text_size <- 8

acc_lvl <- sgn_level1 * sgn_level2 * 50
(panh <- ggplot(res_h, aes(method_an, panti)) +
    geom_bar(stat = 'identity', width = .8, fill = '#f7AD19') +
    geom_hline(yintercept = 0, linewidth = 1, color = 'grey60') +
    geom_hline(yintercept = acc_lvl, linewidth = 1, linetype = '21',
               color = '#E03E18') +
    geom_text(aes(y = pmax(panti + .15, acc_lvl + .15), label = method_an), 
              vjust = .4, hjust = 1, alpha = 1,
              color = 'black', fontface = 'bold', size = 2.7) +
    coord_flip() +
    scale_y_reverse(expand = c(0, 0), limits = c(19.5, 0),
                    breaks = seq(0, 20, 2),
                    labels = paste0(seq(0, 20, 2), '%')) +
    labs(y = 'Percentage of conflicting results (Conflict%)') +
    theme_minimal() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size = axis_title_size),
          axis.text.x = element_text(size = axis_text_size),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.margin = margin(t = 5.5, r = 0, b = 5.5, l = 15))
)

(pnrh <- ggplot(res_h, aes(method_re, prepl)) +
    geom_bar(stat = 'identity', fill = '#70aed1', width = .8) +
    geom_hline(yintercept = 0, linewidth = 1, color = 'grey60') +
    geom_text(aes(y = prepl - 1, label = method_re),
              vjust = .4, hjust = 1,
              color = 'white', fontface = 'bold', size = 2.7) +
    coord_flip() +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 20),
                       labels = paste0(seq(0, 100, 20), '%'),
                       limits = c(0, 100)) +
    labs(y = 'Replication percentage (Replication%)') +
    theme_minimal() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size = axis_title_size),
          axis.text.x = element_text(size = axis_text_size),
          legend.position = 'none',
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.margin = margin(t = 5.5, r = 15, b = 5.5, l = 30))
)


(poph <- ggplot(res_h, aes(method_op, popp)) +
    geom_bar(stat = 'identity', fill = '#f7AD19', width = .8) +
    geom_hline(yintercept = 0, linewidth = 1, color = 'grey60') +
    geom_text(aes(y = pmax(popp + .4, sgn_level1 * 50 + .4), label = method_op), 
              vjust = .4, hjust = 1,
              color = 'black', fontface = 'bold', size = 2.1) +
    geom_hline(yintercept = sgn_level1 * 50, linewidth = 1, linetype = '21',
               color = '#E03E18') +
    coord_flip() +
    scale_y_reverse(expand = c(0, 0), breaks = seq(0, 50, 5),
                    labels = paste0(seq(0, 50, 5), '%'),
                    limits = c(49, 0)) +
    labs(y = 'Percentage of opposite estimates (Opposite%)') +
    theme_minimal() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size = axis_title_size),
          axis.text.x = element_text(size = axis_text_size),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.margin = margin(t = 5.5, r = 0, b = 5.5, l = 15))
)

(pcih <- ggplot(res_h %>% filter(pcvg > 0), aes(method_ci, pcvg)) +
    geom_bar(stat = 'identity', fill = '#70aed1', width = .8) +
    geom_hline(yintercept = 0, linewidth = 1, color = 'grey60') +
    geom_text(aes(y = pmin(pcvg - 1, 88.5), label = method_ci),
              vjust = .4, hjust = 1,
              color = 'white', fontface = 'bold', size = 2.1) +
    geom_hline(yintercept = 90, linewidth = 1, linetype = '21',
               color = '#E03E18', alpha = .50) +
    geom_hline(yintercept = 95, linewidth = 1, linetype = '21',
               color = '#E03E18') +
    coord_flip() +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 100, 20),
                       labels = paste0(seq(0, 100, 20), '%'),
                       limits = c(0, 100)) +
    labs(y = 'Overlap percentage of CIs (CI%)') +
    theme_minimal() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size = axis_title_size),
          axis.text.x = element_text(size = axis_text_size),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.margin = margin(t = 5.5, r = 15, b = 5.5, l = 30))
)

(pnhh <- ggplot(res_h, aes(method_nh, nhit)) +
    geom_bar(stat = 'identity', fill = '#70aed1', width = .8) +
    geom_hline(yintercept = 0, linewidth = 1, color = 'grey60') +
    geom_text(aes(y = nhit - 30, label = method_nh),
              vjust = .4, hjust = 1,
              color = 'white', fontface = 'bold', size = 2.1) +
    coord_flip() +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 10000, 1000),
                       limits = c(0, max(hits_h$nhit))) +
    labs(y = 'Number of significant taxa (NHits)') +
    theme_minimal() +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size = axis_title_size),
          axis.text.x = element_text(size = axis_text_size),
          panel.grid.minor.y = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.margin = margin(t = 5.5, r = 15, b = 5.5, l = 30))
)

pmainh <- ggplot(tibble(x = 0:1)) +
  labs(title = paste0('Figure A5: Results of the split-data analyses with covariates')) +
  theme_void() +
  theme(plot.title = element_text(hjust = 0, face = 'bold', size = title_size))


(pmajorh <- ggarr(panh, pblank, pnrh, ncol = 3, widths = c(1, .1, 1.1),
                  labels = c('a)', '', 'b)', ''),
                  font.label = list(size = title_size)))

(pminorh <- ggarr(poph, pblank, pcih, ncol = 3, widths = c(1, .1, 1.1),
                  labels = c('c)', '', 'd)', ''),
                  font.label = list(size = title_size)))

psensh <- ggarr(pnhh, ncol = 1, labels = 'e)',
                font.label = list(size = title_size))

(ph <- ggarr(pmainh, pmajorh, pblank, pminorh, pblank, psensh, nrow = 6,
             heights = c(.1, .6, .05, .4, .05, .4)))

ggsave(filename = paste0('fig_A5_', date,'.png'),
       plot = ph, width = 170, height = 180, dpi = 300, unit = 'mm',
       bg = 'white')

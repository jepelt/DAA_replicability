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
                  'ANCOM-BC2 (pc = 0)',
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
                  'MaAsLin2 (LOG)',
                  'MaAsLin2 (AST)',
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
  .[, method := ifelse(method == 'MaAsLin2 (AST)', 'MaAsLin2-AST', method)] %>% 
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
  filter(site == 'stool')


#Calculate values of the evaluation metrics-------------------------------------
sgn_level2s <- .05

ress <- list()
alphas <- c(.01, .05, .10)
for(i in 1:length(alphas)){
  sgn_level1s <- alphas[i]
  
  ol1s <- res_all %>% 
    filter(analysis == 'whole' & n_case >= 10 & n_control >= 10 & 
             disease != 'HIV') %>% 
    mutate(n1 = n + avg_ls / 10 ^ 10) %>% 
    select(method, type, disease, taxon, study1 = study, n1,
           pco1 = prevl_control, pca1 = prevl_case,
           est1 = est, q, lwr1 = lwr, upr1 = upr)
  
  ol2s <- res_all %>% 
    filter(analysis == 'whole') %>% 
    mutate(is2 = T,
           n2 = n + avg_ls / 10 ^ 10) %>% 
    select(type, disease, taxon, method, study2 = study, n2, is2,
           pco2 = prevl_control, pca2 = prevl_case,
           est2 = est, p, lwr2 = lwr, upr2 = upr)
  
  ols <- left_join(ol1s, ol2s, 
                   by = join_by(method, type, disease, taxon, n1 < n2)) %>% 
    filter(q < sgn_level1s & !is.na(is2)) %>% 
    mutate(sgn1 = q < sgn_level1s,
           d1 = ifelse(is.na(est1) | is.nan(est1), 0, sign(round(est1, 8))),
           s1 = if_else(sgn1 & d1 != 0, T, F, F),
           sgn2 = p < sgn_level2s,
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
    select(method, type, disease, study1, study2, taxon,
           est1, q, est2, p, sgn1, d1, s1, sgn2, d2, s2,
           repl, same, opposite, anti, nores, is_ci, cvg, check) %>% 
    filter(!(str_detect(study1, 'Thomas') & str_detect(study2, 'Thomas')))
  
  stopifnot(sum(ols$check != 1) == 0)
  
  n_studiess <- nrow(ols %>% distinct(study1))
  
  hits_s <- left_join(ol1s, ol2s, 
                      by = join_by(method, type, disease, taxon, n1 < n2)) %>% 
    distinct(method, study1, taxon, .keep_all = T) %>% 
    filter(q < sgn_level1s) %>% 
    group_by(method) %>% 
    summarize(nhit = n())
  
  res_s <- ols %>% 
    group_by(method, study1, taxon) %>% 
    summarize_at(vars(repl, same, opposite, anti, nores, cvg, is_ci),
                 ~ sum(., na.rm = T)) %>% 
    mutate(n = repl + same + opposite + anti + nores,
           cvg = cvg / is_ci,
           is_ci = ifelse(is_ci > 0.5, 1, 0)) %>% 
    mutate_at(vars(repl, same, opposite, anti, nores), ~ . / n) %>% 
    group_by(method) %>% 
    summarize_at(vars(repl, same, opposite, anti, nores, cvg, is_ci),
                 ~ sum(., na.rm = T)) %>% 
    ungroup() %>% 
    left_join(., hits_s, 'method') %>% 
    mutate(n = repl + same + opposite + anti + nores,
           panti = anti / n * 100,
           prepl = repl / n * 100,
           popp = (anti + opposite + .5 * nores) / n * 100,
           pcvg = cvg / is_ci * 100,
           method_an = fct_reorder(method, -panti),
           method_re = fct_reorder(method, prepl),
           method_op = fct_reorder(method, -popp),
           method_ci = fct_reorder(method, pcvg),
           method_nh = fct_reorder(method, nhit))
  
  
  res_concl <- res_s %>% select(method, panti, prepl, popp, pcvg, nhit) %>% 
    mutate_at(vars(pcvg), ~ ifelse(. == 0, NA, .)) %>% 
    mutate(across(c(panti, prepl, popp, pcvg, nhit),
                  .fns = list(std = ~ as.numeric(scale(.))))) %>% 
    rename_at(vars(panti:nhit), ~ paste0(., '_res')) %>% 
    mutate_at(vars(panti_std, popp_std), ~ - .) %>% 
    pivot_longer(cols = -method) %>% 
    separate(name, into = c('var', 'type')) %>% 
    pivot_wider(names_from = type, values_from = value) %>% 
    mutate(name = factor(var,
                         levels = c('panti','prepl', 'popp', 'pcvg', 'nhit'),
                         labels = c('Conflict%', 'Replication%',
                                    'Opposite%', 'CI%', 'NHits'))) %>% 
    group_by(method) %>%
    mutate(res = ifelse(str_sub(var, 1, 1) == 'p', 
                        paste0(round2(res, 1), '%'),
                        as.character(res)),
           res = ifelse(str_detect(res, 'NA'), 'N/A', res),
           alpha = sgn_level1s)
  
  ress[[i]] <- res_concl
}


#Create Figure A2.3-------------------------------------------------------------
res_final <- ress %>% bind_rows %>% 
  group_by(method) %>%
  arrange(method, alpha, var) %>% 
  mutate(total = weighted.mean(std, w = rep(c(1, 2, .5, .5, 1), 3),
                               na.rm = T)) %>% 
  ungroup() %>% 
  mutate(method = fct_reorder(method, total))

glims <- c(min(c(res_final$std, -2), na.rm = T),
           max(c(res_final$std, 2.1), na.rm = T))

alphas2 <- rep(c('.01', '.05', '.10'), 5)
xt <- expression(paste("Significance level in exploratory datasets (", alpha, ")"))
(pc <- ggplot(res_final, aes(interaction(alpha, name), method, fill = std)) +
    geom_tile() +
    geom_vline(xintercept = c(3.5, 6.5, 9.5, 12.5), color = 'white',
               linewidth = 1) +
    geom_text(aes(label = res), color = 'grey20', size = 2.0, hjust = 0.55) +
    annotate(geom = 'text', x = 2, y = 26.5, vjust = 1, hjust = .5, color = 'grey0',
             label = 'Percentage of\nconflicting results\n(Conflict%)', size = 2.5) +
    annotate(geom = 'text', x = 5, y = 26.5, vjust = 1, hjust = .5, color = 'grey0',
             label = 'Replication\npercentage\n(Replication%)', size = 2.5) +
    annotate(geom = 'text', x = 8, y = 26.5, vjust = 1, hjust = .5, color = 'grey0',
             label = 'Percentage of\nopposite estimates\n(Opposite%)', size = 2.5) +
    annotate(geom = 'text', x = 11, y = 26.5, vjust = 1, hjust = .5, color = 'grey0',
             label = 'Overlap percentage\nof conf. intervals\n(CI%)', size = 2.5) +
    annotate(geom = 'text', x = 14, y = 26.5, vjust = 1, hjust = .5, color = 'grey0',
             label = 'Number of\nsignificant taxa\n(NHits)', size = 2.5) +
    annotate(geom = 'text', x = 1:15, y = 0, vjust = 0, hjust = .5, color = 'grey20',
             label = alphas2, size = 2.5) +
    scale_fill_gradient2(low = "#f7AD19", mid = "white", high = "#70aed1",
                         midpoint = 0, na.value = 'white',
                         breaks = -2:2,
                         limits = glims,
                         labels = c('-2 (Worse)', '-1', '0', '1', '2 (Better)'),
                         name = 'Standardized\nvalue') +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(title = 'Figure A2.3: Results of the separate study analyses', x = xt) +
    theme_minimal() +
    theme(axis.title.x = element_text(size = 10, color = 'grey0'),
          axis.text.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 8),
          legend.title = element_text(size = 10, color = 'grey0'),
          legend.text = element_text(size = 8),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0, size = 12, face = 'bold',
                                    margin = margin(0, 0, 10, 0)),
          plot.title.position = 'plot')
)

ggsave(filename = paste0('fig_A2.3_', date,'.png'),
       plot = pc, width = 170, height = 170, dpi = 300, unit = 'mm',
       bg = 'white')

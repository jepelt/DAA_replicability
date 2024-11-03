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
  filter(site == 'stool')


#Summarize results--------------------------------------------------------------
condition <- 6000

lim_alpha <- res_all %>% 
  mutate(q = ifelse(method == 'ANCOM-BC2', q20, q)) %>%
  filter(analysis == 'half' & half == 1 & !is.na(q)) %>% 
  group_by(method) %>%
  arrange(q) %>%
  slice(condition + 1) %>%
  ungroup() %>% 
  select(method, alpha = q)

res_exploratory <- res_all %>% 
  filter(analysis == 'half' & half == 1) %>% 
  mutate(q = ifelse(method == 'ANCOM-BC2', q20, q)) %>%
  left_join(., lim_alpha, by = 'method') %>% 
  filter(q < alpha) %>% 
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
  
hits <- res_exploratory %>% group_by(method) %>% summarize(r_nhit = n())
  
res_summary <- res_combined %>% 
  group_by(method) %>% 
  summarize_at(vars(repl, same, opposite, anti, nores),
               ~ sum(., na.rm = T)) %>% 
  ungroup() %>% 
  left_join(., hits, 'method') %>% 
  mutate(n = repl + same + opposite + anti + nores,
         r_confl = anti / n * 100,
         r_repl = repl / n * 100,
         condition = condition) %>% 
  left_join(., lim_alpha, 'method') %>% 
  select(condition, method, r_confl, r_repl, r_nhit, r_alpha = alpha)


#Create Figure------------------------------------------------------------------
date <- format(Sys.Date(), "%d%m%y")
round2 <- function(x, decimals = 2){
  format(round(x, decimals), nsmall = decimals)
}

res_figure <- bind_rows(res_summary) %>%
  mutate(r_nhit = condition,
         s_confl = -as.numeric(scale(sqrt(r_confl))),
         s_repl = as.numeric(scale(r_repl)),
         s_alpha = 0) %>% 
  ungroup() %>% 
  select(method, condition,
         r_confl, r_repl, r_alpha, s_confl, s_repl, s_alpha) %>% 
  pivot_longer(cols = -c(method, condition)) %>% 
  separate(name, into = c('type', 'var')) %>% 
  pivot_wider(names_from = type, values_from = value) %>% 
  mutate(name = factor(var,
                       levels = c('confl','repl', 'alpha'),
                       labels = c('Conflict%', 'Replication%', 'Nominal\nFDR level'))) %>% 
  group_by(method) %>% 
  mutate(total = mean(s)) %>% 
  ungroup() %>% 
  mutate(method = fct_reorder(method, total),
         res = case_when(var == 'alpha' ~ round2(r, 2),
                         var == 'confl' ~ paste0(round2(r, 1), '%'),
                         T ~ paste0(round2(r, 0), '%')))


(p1 <- ggplot(res_figure, aes(name, method, fill = s)) +
    geom_tile() + #linewidth = .1, color = 'grey40') +
    geom_text(aes(label = res), color = 'grey20', size = 2.5, hjust = 0.55) +
    geom_vline(xintercept = c(4.5, 8.5), color = 'white',
               linewidth = 1) +
    scale_fill_gradient2(low = "#f7AD19", mid = "white", high = "#70aed1",
                         midpoint = 0, na.value = 'white',
                         breaks = -2:1,
                         labels = c('-2 (Worse)', '-1', '0', '1 (Better)'),
                         name = 'Standardized\nvalue') +
    scale_x_discrete(expand = c(0, 0), position = 'top') +
    scale_y_discrete(expand = c(0, 0)) +
    #labs(title = 'Figure A3: Consistency on constant sensitivity\nin the split-data analyses (NHits = 6000)') +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x =element_text(size = 7),
          axis.title.y = element_blank(),
          axis.text.y = element_text(size = 7),
          legend.title = element_text(size = 8, color = 'grey0'),
          legend.text = element_text(size = 7),
          panel.spacing = unit(2, "lines"),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0, size = 12, face = 'bold',
                                    margin = margin(0, 0, 10, 0)),
          plot.title.position = 'plot')
)


ggsave(filename = paste0('fig_A3_nhits_6000_', date,'.png'),
       plot = p1, width = 100, height = 120, dpi = 300, unit = 'mm',
       bg = 'white')

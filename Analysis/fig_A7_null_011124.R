library(tidyverse)
library(data.table)


#Load files and prepare tibble of results---------------------------------------
load('res_permuted_301024.rds')
load('data_meta_041023.rds')
load('data_prevl_061023.rds')

res_all <- bind_rows(res_list) %>% 
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
  .[, method := tstrsplit(method, " ", fixed = T, keep = 1)]


#Create figures ... ------------------------------------------------------------
res_hist <- res_all %>% 
  mutate(type = ifelse(str_detect(data_id, '16S'), '16S', 'SG'),
         significance = ifelse(p < .05, 'p < .05', 'p > .05')) %>% 
  filter(!(prevl_control == 1 & prevl_case == 1 & method == 'LogR')) %>% 
  filter(!(!ss_pass_p & method == 'ANCOM-BC2')) %>% 
  filter(!is.na(p))

d_lines <- res_hist %>% 
  group_by(method) %>% 
  summarize(y = n() / 100)

cut_text <- d_lines %>%
  filter(method == 'edgeR') %>% 
  mutate(x = .15, y = 4900, significance = 'p > .05')

(p <- ggplot(res_hist, aes(p, fill = significance)) + 
    geom_histogram(breaks = seq(0, 1, .01)) +
    geom_hline(data = d_lines, aes(yintercept = y), color = '#E03E18',
               linewidth = .7, linetype = '33') +
    geom_text(data = cut_text, aes(x, y), label = "y-axis cut", hjust = 0) +
    scale_fill_manual(values = c('blue', 'grey60')) +
    facet_wrap(~ method, ncol = 5) +
    coord_cartesian(ylim = c(0, 5000)) +
    labs(#title = 'Figure A7.1: Distribution of p-values on datasets with randomly permuted group labels',
      y = 'Number of taxa', x = 'p-value', fill = 'Significance') +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, vjust = .40, size = 7),
          axis.title.x = element_text(size = 10),
          axis.text.y = element_text(size = 8),
          axis.title.y = element_text(size = 10),
          strip.text = element_text(size = 8),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8),
          legend.position = 'bottom',
          plot.title = element_text(face = 'bold', size = 11),
          plot.title.position = 'plot')
)

date <- format(Sys.Date(), "%d%m%y")
ggsave(filename = paste0('fig_A7_null_', date,'.png'),
       plot = p, width = 170, height = 200, dpi = 300, unit = 'mm',
       bg = 'white')

res_hist %>% group_by(method) %>% summarize(s = sum(p < .01, na.rm = T))
res_hist %>% group_by(method) %>% summarize(m = 100 * mean(p < .05, na.rm = T))


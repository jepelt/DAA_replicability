library(tidyverse)

load('data_abs_200324.rds')
load('res_absolute_011124.rds')

meta <- data_abs[[1]]$meta
counts_abs <- data_abs[[1]]$counts_abs

true_means <- cbind(meta, counts_abs) %>% 
  group_by(group) %>% 
  summarize_if(is.numeric, mean) %>% 
  column_to_rownames('group') %>% 
  t()

true_sign <- sign(true_means[, 'case'] - true_means[, 'control'])

res_true <- tibble(taxon = names(true_sign),
                   sign_true = true_sign)

res_abs <- bind_rows(res_list) %>% 
  left_join(., res_true, by = 'taxon') %>% 
  mutate(sign_est = sign(est)) %>% 
  select(method, taxon, sign_est, sign_true) %>% 
  group_by(method) %>% 
  summarize(acc = mean(sign_est == sign_true, na.rm = T)) %>% 
  filter(!str_detect(method, 'LRT|Score|ev = TRUE|Wilcoxon0')) %>% 
  mutate(norm = case_when(str_detect(method, 'TSS|corncob|NB|LDM') ~ 'TSS',
                          # str_detect(method, 'ALDEx|ANCOM|LinDA|CLR') ~ 'CLR',
                          T ~ 'Other'),
         norm = factor(norm, levels = c('TSS', 'Other')),
         method2 = factor(tolower(method),
                          labels = 
                            c('ALDEx2 (CLR)', 'ALDEx2 (scaled CLR)',
                              'ANCOM-BC2 (Bias corr)',
                              'corncob (TSS)',
                              'DESeq2 (RLE)', 'DESeq2 (GMPR)', 'DESeq2 (TSS)', 'DESeq2 (Wrench)',
                              'edgeR (RLE)', 'edgeR (TMM)', 'edgeR (TMMwsp)',
                              'fastANCOM (Ref. taxa)',
                              'LDM (TSS)',
                              'limma-voom (RLE)', 'limma-voom (TMM)',  'limma-voom (TMMwsp)',
                              'LinDA (Bias corr.)',
                              'LogR',
                              'MaAsLin2-AST (TSS)', 'MaAsLin2 (CSS)', 'MaAsLin2 (TMM)',
                              'MaAsLin2 (TSS)', 'MaAsLin2 (CLR)',
                              'metagenomeSeq (CSS)', 'metagenomeSeq (GMPR)',
                              'metagenomeSeq (TSS)', 'metagenomeSeq (Wrench)',
                              'NegBin (TSS)',
                              'ORM/Wilcox (GMPR)', 'ORM/Wilcox (TSS)', 'ORM/Wilcox (Wrench)',
                              'radEmu (Median DA)',
                              'ZicoSeq (Ref. taxa)')),
         acc_main = ifelse(method2 %in% c('ALDEx2 (CLR)',
                                          'ANCOM-BC2 (Bias corr)',
                                          'corncob (TSS)',
                                          'DESeq2 (RLE)',
                                          'edgeR (TMM)',
                                          'fastANCOM (Ref. taxa)',
                                          'LDM (TSS)',
                                          'limma-voom (TMM)',
                                          'LinDA (Bias corr.)',
                                          'LogR',
                                          'MaAsLin2 (TSS)',
                                          'metagenomeSeq (CSS)',
                                          'ORM/Wilcox (TSS)',
                                          'ZicoSeq (Ref. taxa)'
         ), acc, NA))

ff <- ifelse(1:33 %in% c(1, 3, 4, 5, 10, 12, 13, 15, 17, 18, 22, 24, 30, 33),
            'bold', 'plain')

(pabs <- ggplot(res_abs, aes(method2, acc, color = norm)) +
  geom_vline(xintercept = 
               c(2, 3, 4, 8, 11, 12, 13, 16, 17, 18, 23, 27, 28, 31, 32) + .5,
             color = 'grey40') +
  geom_hline(yintercept = 0.5, linetype = '44', color = 'grey40') +
  geom_point(size = 2) +
  geom_point(aes(y = acc_main), size = 3) +
  scale_y_continuous(breaks = seq(.1, .9, .1),
                     labels = scales::percent_format()) +
  labs(y = 'Accuracy', color = 'Norm.\ntype',
       # title = 'Figure A15: Accuracy of estimating the sign of absolute DA',
       #caption = '(Absolute DA = DA estimate based on measured absolute abundances)'
       ) +
  ggsci::scale_color_jco() +
  theme_light() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = -90, vjust = 0.2, hjust = 0, size = 6,
                                   face = ff),
        axis.title.y = element_text(size = 10),
        axis.text.y = element_text(size = 8),
        legend.title = element_text(),
        legend.position = 'right',
        plot.title.position = "plot",
        plot.title = element_text(face = 'bold', size = 12),
        plot.caption.position = "plot")
)

date <- format(Sys.Date(), "%d%m%y")
ggsave(filename = paste0('figure_A15_abs_', date,'.png'),
       plot = pabs, width = 170, height = 100, dpi = 300, units = 'mm',
       bg = 'white')

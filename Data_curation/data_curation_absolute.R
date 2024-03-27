library(tidyverse)
library(data.table)

data_abs <- list()

#Barlow data--------------------------------------------------------------------
rel_raw <- read.csv('Relative_Abundance_Table.csv') %>% 
  rename_at(vars(matches('__')), ~ paste0('genus_', 1:142)) %>% 
  filter(Site == 'Stool') %>% 
  mutate(Sample = str_replace_all(Sample, ' ', '_'))

counts0 <- rel_raw %>% 
  select(Sample, contains('genus')) %>% 
  column_to_rownames('Sample') %>% 
  mutate_all(~ . * 19500 / 100) %>% 
  mutate_all(as.integer) %>% 
  as.matrix()

counts <- counts0[, colMeans(counts0 > .5) >= .1] 


abs_raw <- read.csv('Absolute_Abundance_Table.csv') %>% 
  rename_at(vars(matches('__')), ~ paste0('genus_', 1:142)) %>% 
  filter(Site == 'Stool') %>% 
  mutate(Sample = str_replace_all(X, ' ', '_'))

counts_abs0 <- abs_raw %>% 
  select(Sample, contains('genus')) %>% 
  column_to_rownames('Sample') %>% 
  as.matrix()

counts_abs <- counts_abs0[, colMeans(counts_abs0 > .5) >= .1] 


meta <- rel_raw %>% 
  transmute(data_id = 'Barlow',
            group = factor(Diet, levels = c('Control', 'Keto'),
                           labels = c('control', 'case')),
            id = Sample) %>% 
  column_to_rownames('id') %>% 
  as.data.frame()

data_abs[[1]] <- list(meta = meta, counts = counts, counts_abs = counts_abs)


#Vieira-Silva data--------------------------------------------------------------
abs_raw <- read_tsv('QMP.matrix.tsv')
meta_raw <- readxl::read_xlsx('viera_meta.xlsx', sheet = 2)

meta_ids <- meta_raw %>% 
  select(id = 1, disease = 2) %>% 
  filter(disease %in% c('CD', 'mHC')) %>% 
  select(-disease)
abundance_ids <- abs_raw %>% select(id = 1)
ids <- inner_join(meta_ids, abundance_ids, 'id')


counts_abs0 <- abs_raw %>% 
  rename(id = '...1') %>% 
  left_join(ids, ., 'id') %>% 
  column_to_rownames('id') %>% 
  rename_all(~ paste0('genus_', 1:303)) %>% 
  as.matrix()

counts_abs <- counts_abs0[, colMeans(counts_abs0 > .5) >= .1]


librs <- rep(0, nrow(counts_abs0))
for(i in 1:nrow(counts_abs0)){
  x <- counts_abs0[i, ]
  librs[i] <- min(x[x > 0]) 
}

counts0 <- counts_abs / librs

counts <- counts0[, colMeans(counts0 > .5) >= .1]
mode(counts) <- 'integer'


meta <- meta_raw %>% 
  select(id = 1, group = 2) %>% 
  mutate(data_id = 'Vieira') %>% 
  left_join(ids, ., 'id') %>% 
  mutate(group = factor(group, levels = c('mHC', 'CD'),
                        labels = c('control', 'case'))) %>% 
  column_to_rownames('id') %>% 
  as.data.frame()

data_abs[[2]] <- list(meta = meta, counts = counts, counts_abs = counts_abs)


#Save data----------------------------------------------------------------------
save(data_abs, file = 'data_abs_200324.rds')

library(tidyverse)

load('data_171023.rds')

studies_sep <- map_dfr(d_all, ~ .$meta %>% slice(1)) %>% 
  remove_rownames() %>% 
  filter(site == 'stool') %>% 
  filter(str_detect(data_id, 'whole')) %>% 
  filter(!disease %in% c('HIV', 'schizophrenia', 'NASH')) %>% 
  filter(n_control >= 10 & n_case >= 10) %>% 
  pull(data_id)

d_sep <- keep(d_all, ~ .$meta$data_id[1] %in% studies_sep) 


d_permuted <- list()
l <- 0

set.seed(1)
for(i in 1:length(d_sep)){
  meta <- d_sep[[i]]$meta
  counts <- d_sep[[i]]$counts
  
  for(j in 1:10){
    permuted_ind <- sample(1:nrow(meta))
    meta_perm <- meta %>% 
      mutate_at(vars(group, bmi, age, sex), ~ .[permuted_ind]) %>%
      as.data.frame
    
    l <- l + 1
    d_permuted[[l]] <- list(meta = meta_perm, counts = counts)
  }
}

#save(d_permuted, file = 'data_permuted_171024.rds')

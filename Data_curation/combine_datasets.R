library(tidyverse)

load('data_sg_half_270923.rds')
d_sg_half <- keep(data_half, ~ .$meta$iter[1] < 5.5) 

load('data_sg_whole_041023.rds')
d_sg_whole <- data_whole

load('data_16s_half_031023.rds')
d_16s_half <- keep(data_half, ~ .$meta$iter[1] < 5.5)

load('data_16s_whole_031023.rds')
d_16s_whole <- data_whole


#Select and combine datasets used in split-data analyses and separate study
#analyses as one file-----------------------------------------------------------

diseases_sg <- d_sg_whole %>% 
  map(~ .$meta) %>% 
  bind_rows %>% 
  filter(!str_detect(study, 'ChngKR_2016')) %>% 
  distinct(study, .keep_all = T) %>% 
  group_by(disease) %>% 
  summarize(n = n()) %>% 
  filter(n >= 2) %>% 
  pull(disease)

d_sg_whole2 <- keep(d_sg_whole, ~ .$meta$disease[1] %in% diseases_sg)


diseases_16s <- d_16s_whole %>% 
  map(~ .$meta) %>% 
  bind_rows %>% 
  distinct(study, .keep_all = T) %>% 
  group_by(disease) %>% 
  summarize(n = n()) %>% 
  filter(n >= 2 | disease == 'EDD') %>% 
  pull(disease)

d_16s_whole2 <- keep(d_16s_whole, ~ .$meta$disease[1] %in% diseases_16s)


d_all <- c(d_16s_half, d_16s_whole2, d_sg_half, d_sg_whole2)

#save(d_all, file = 'data_171023.rds')



#Create files including all metadata--------------------------------------------
datasets <- c(sg_half, sg_whole, s16_half, s16_whole)

calc_prevls <- function(l){
  counts <- l$counts
  meta <- l$meta
  
  prevl <- colMeans(counts > .5)
  prevl_control <- colMeans(counts[meta$group == 'control', ] > .5)
  prevl_case <- colMeans(counts[meta$group == 'case', ] > .5)
  mean_count <- colMeans(counts)
  
  prevls <- cbind(prevl, prevl_control, prevl_case, mean_count) %>% 
    as.data.frame() %>% 
    rownames_to_column('taxon') %>% 
    mutate(data_id = meta$data_id[1])
  
  return(prevls)
}


metas <- map(datasets, ~ .$meta %>% dplyr::slice(1))

d_meta <- bind_rows(metas) %>% 
  remove_rownames() %>% 
  select(-c(bmi, age, sex, group)) 

#save(d_meta, file = 'data_meta_041023.rds')


prevalences <- map(datasets, calc_prevls)

d_prevl <- bind_rows(prevalences)

#save(d_prevl, file = 'data_prevl_061023.rds')

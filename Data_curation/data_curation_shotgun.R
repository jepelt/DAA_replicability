library(tidyverse)

smd <- curatedMetagenomicData::sampleMetadata

two_conds <- smd %>% 
  group_by(study_name) %>%
  distinct(study_condition) %>% 
  summarize(groups = n()) %>% 
  filter(groups > 1.5)

sites <- smd %>% 
  group_by(study_name) %>%
  distinct(body_site) %>% 
  summarize(sites = n())

nsa <- smd %>% 
  group_by(study_name) %>%
  summarize(samples = n())

nsu <- smd %>% 
  group_by(study_name) %>%
  distinct(subject_id) %>% 
  summarize(subjects = n())

nc <- smd %>% 
  filter(study_condition == 'control') %>% 
  group_by(study_name) %>%
  distinct(subject_id) %>% 
  summarize(control = n())

nc2 <- smd %>% 
  filter(study_condition != 'control') %>% 
  group_by(study_name) %>%
  distinct(subject_id, .keep_all = T) %>% 
  summarize(non_control = n(),
            site = first(body_site))

studies <- inner_join(two_conds, nsa, 'study_name') %>% 
  left_join(., nsu, 'study_name') %>% 
  left_join(., nc, 'study_name') %>% 
  left_join(., nc2, 'study_name') %>% 
  left_join(., sites, 'study_name') %>% 
  filter(control >= 10 & non_control >= 10)

rm(nc, nc2, nsa, nsu, sites, two_conds)


#Curate data--------------------------------------------------------------------
d1 <- smd %>% 
  filter(study_name == "BedarfJR_2017") %>% 
  mutate(dis = 'PD')

d2 <- smd %>% 
  filter(study_name == "Castro-NallarE_2015") %>% 
  mutate(dis = 'schizophrenia')

d3_0 <- smd %>% 
  filter(study_name == "ChngKR_2016") %>% 
  arrange(sample_id) %>% 
  group_by(subject_id) %>% 
  mutate(i = 1:2,
         dis = 'AD') %>% 
  ungroup()

d3_1 <- d3_0 %>% filter(i == 1) %>% mutate(study_name = "ChngKR_2016_1")
d3_2 <- d3_0 %>% filter(i == 2) %>% mutate(study_name = "ChngKR_2016_2")
rm(d3_0)

d4_1 <- smd %>% 
  filter(study_name == 'FengQ_2015' 
         & study_condition %in% c('control', 'CRC')) %>% 
  mutate(study_name = 'FengQ_2015_CRC',
         dis = 'CRC')

d4_2 <- smd %>% 
  filter(study_name == 'FengQ_2015' 
         & study_condition %in% c('control', 'adenoma')) %>% 
  mutate(study_name = 'FengQ_2015_adenoma',
         dis = 'adenoma')

d5_1 <- smd %>% 
  filter(study_name == 'GhensiP_2019' 
         & study_condition %in% c('control', 'mucositis')) %>% 
  group_by(subject_id) %>% 
  arrange(desc(study_condition)) %>%
  filter(row_number() == 1) %>% 
  ungroup() %>% 
  mutate(study_name = 'GhensiP_2019_mucositis',
         dis = 'mucositis')

d5_2 <- smd %>% 
  filter(study_name == 'GhensiP_2019' 
         & study_condition %in% c('control', 'peri-implantitis')) %>% 
  group_by(subject_id) %>% 
  arrange(desc(study_condition)) %>%
  filter(row_number() == 1) %>% 
  ungroup() %>% 
  mutate(study_name = 'GhensiP_2019_peri-implatitis',
         dis = 'peri-implantitis')

d6 <- smd %>% 
  filter(study_name == "GuptaA_2019") %>% 
  mutate(dis = 'CRC')

d7 <- smd %>% 
  filter(study_name == "HMP_2019_ibdmdb" & days_from_first_collection == 0) %>% 
  distinct(subject_id, .keep_all = T) %>% 
  mutate(dis = 'IBD')

d8 <- smd %>% 
  filter(study_name == "HallAB_2017") %>%
  group_by(subject_id) %>% 
  arrange(days_from_first_collection) %>% 
  filter(row_number() == 1) %>% 
  ungroup() %>% 
  mutate(dis = 'IBD')

d9_1 <- smd %>% 
  filter(study_name == 'HanniganGD_2017' 
         & study_condition %in% c('control', 'CRC')) %>% 
  mutate(study_name = 'HanniganGD_2017_CRC',
         dis = 'CRC')

d9_2 <- smd %>% 
  filter(study_name == 'HanniganGD_2017' 
         & study_condition %in% c('control', 'adenoma')) %>% 
  mutate(study_name = 'HanniganGD_2017_adenoma',
         dis = 'adenoma')

d10 <- smd %>% 
  filter(study_name == "Heitz-BuschartA_2016") %>%
  group_by(subject_id) %>% 
  arrange(days_from_first_collection) %>% 
  filter(row_number() == 1) %>% 
  ungroup() %>% 
  mutate(dis = 'T1D')

d11 <- smd %>% 
  filter(study_name == "IjazUZ_2017") %>%
  group_by(subject_id) %>% 
  arrange(days_from_first_collection) %>% 
  filter(row_number() == 1) %>% 
  ungroup() %>% 
  mutate(dis = 'IBD')

d12 <- smd %>% 
  filter(study_name == "JieZ_2017") %>% 
  mutate(dis = 'ACVD')

d13_1 <- smd %>% 
  filter(study_name == 'KarlssonFH_2013' 
         & study_condition %in% c('control', 'IGT')) %>% 
  mutate(study_name = 'KarlssonFH_2013_IGT',
         dis = 'IGT')

d13_2 <- smd %>% 
  filter(study_name == 'KarlssonFH_2013' 
         & study_condition %in% c('control', 'T2D')) %>% 
  mutate(study_name = 'KarlssonFH_2013_T2D',
         dis = 'T2D')

d14_1 <- smd %>% 
  filter(study_name == 'LiJ_2014' 
         & study_condition %in% c('control', 'IBD')) %>% 
  distinct(subject_id, .keep_all = T) %>% 
  mutate(study_name = 'LiJ_2014_IBD',
         dis = 'IBD')

d14_2 <- smd %>% 
  filter(study_name == 'LiJ_2014' 
         & study_condition %in% c('control', 'T1D')) %>% 
  distinct(subject_id, .keep_all = T) %>% 
  mutate(study_name = 'LiJ_2014_T1D',
         dis = 'T1D')

d14_3 <- smd %>% 
  filter(study_name == 'LiJ_2014' 
         & study_condition %in% c('control', 'T2D')) %>% 
  distinct(subject_id, .keep_all = T) %>% 
  mutate(study_name = 'LiJ_2014_T2D',
         dis = 'T2D')

d15_1 <- smd %>% 
  filter(study_name == 'LiJ_2017' 
         & study_condition %in% c('control', 'hypertension')) %>% 
  distinct(subject_id, .keep_all = T) %>% 
  mutate(study_name = 'LiJ_2017_hypertension',
         dis = 'hypertension')

d15_2 <- smd %>% 
  filter(study_name == 'LiJ_2017' 
         & study_condition %in% c('control', 'pre-hypertension')) %>% 
  distinct(subject_id, .keep_all = T) %>% 
  mutate(study_name = 'LiJ_2017_pre-hypertension',
         dis = 'pre-hypertension')

d16 <- smd %>% 
  filter(study_name == 'LoombaR_2017') %>% #Different groups!
  mutate(study_condition = ifelse(disease_stage %in% c('0', 'I', 'II'),
                                  'control', 'advanced fibrosis'),
         dis = 'advanced fibrosis')

#MS = metabolic disorder
#Control group includes also MS controls etc
#No counts! (as no number of reads)
# d17_1 <- smd %>% 
#   filter(study_name == 'MetaCardis_2020_a' 
#          & study_condition %in% c('control', 'CAD')) %>% 
#   distinct(subject_id, .keep_all = T) %>% 
#   mutate(study_name = 'MetaCardis_2020_a_CAD',
#          dis = 'CAD')
# 
# d17_2 <- smd %>% 
#   filter(study_name == 'MetaCardis_2020_a' 
#          & study_condition %in% c('control', 'HF')) %>% 
#   distinct(subject_id, .keep_all = T) %>% 
#   mutate(study_name = 'MetaCardis_2020_a_HF',
#          dis = 'HF')
# 
# d17_3 <- smd %>% 
#   filter(study_name == 'MetaCardis_2020_a' 
#          & study_condition %in% c('control', 'IGT')) %>% 
#   distinct(subject_id, .keep_all = T) %>% 
#   mutate(study_name = 'MetaCardis_2020_a_IGT',
#          dis = 'IGT')
# 
# d17_4 <- smd %>% 
#   filter(study_name == 'MetaCardis_2020_a' 
#          & study_condition %in% c('control', 'T2D')) %>% 
#   distinct(subject_id, .keep_all = T) %>% 
#   mutate(study_name = 'MetaCardis_2020_a_T2D',
#          dis = 'T2D')

d18 <- smd %>% 
  filter(study_name == "NagySzakalD_2017") %>% 
  mutate(dis = 'ME/CFS')

d19 <- smd %>% 
  filter(study_name == "NielsenHB_2014" & days_from_first_collection == 0) %>% 
  mutate(dis = 'IBD')
  
d20 <- smd %>% 
  filter(study_name == "QinJ_2012" & !is.na(study_condition)) %>% 
  mutate(dis = 'T2D')

d21 <- smd %>% 
  filter(study_name == "QinN_2014") %>% 
  mutate(dis = 'cirrhosis')

d22 <- smd %>% 
  filter(study_name == "RaymondF_2016") %>% 
  mutate(dis = 'cephalosporins')

d23 <- smd %>% 
  filter(study_name == "RubelMA_2020") %>% 
  mutate(dis = 'STH')

d24 <- smd %>% 
  filter(study_name == "SankaranarayananK_2015") %>% 
  mutate(dis = 'T2D')

#The same subjects as controls and cases
# d25 <- smd %>% 
#   filter(study_name == "TettAJ_2016") %>% 
#   mutate(dis = 'T2D')

d26_1 <- smd %>% 
  filter(study_name == 'ThomasAM_2018a' 
         & study_condition %in% c('control', 'CRC')) %>% 
  mutate(study_name = 'ThomasAM_2018a_CRC',
         dis = 'CRC')

d26_2 <- smd %>% 
  filter(study_name == 'ThomasAM_2018a' 
         & study_condition %in% c('control', 'adenoma')) %>% 
  mutate(study_name = 'ThomasAM_2018a_adenoma',
         dis = 'adenoma')

d27 <- smd %>% 
  filter(study_name == "ThomasAM_2018b") %>% 
  mutate(dis = 'CRC')

d28 <- smd %>% 
  filter(study_name == "ThomasAM_2019_c") %>% 
  mutate(dis = 'CRC')

#The same subjects (infants) have several study conditions. Very different ages
#if filtered.
# d29 <- smd %>% 
#   filter(study_name == 'VatanenT_2016' 
#          & study_condition %in% c('control', 'otitis')) %>%
#   group_by(subject_id) %>% 
#   arrange(desc(study_condition), days_from_first_collection) %>% 
#   filter(row_number() == 1) %>% 
#   ungroup() %>% 
#   mutate(dis = 'otitis')

d30 <- smd %>% 
  filter(study_name == "VogtmannE_2016" & !is.na(study_condition)) %>% 
  mutate(dis = 'CRC')

d31 <- smd %>% 
  filter(study_name == "WirbelJ_2018") %>% 
  mutate(dis = 'CRC')

d32_1 <- smd %>% 
  filter(study_name == 'XieH_2016' 
         & study_condition %in% c('control', 'asthma')) %>% 
  mutate(study_name = 'XieH_2016_asthma',
         dis = 'asthma')

d32_2 <- smd %>% 
  filter(study_name == 'XieH_2016' 
         & study_condition %in% c('control', 'migraine')) %>% 
  mutate(study_name = 'XieH_2016_migraine',
         dis = 'migraine')

d33_1 <- smd %>% 
  filter(study_name == 'YachidaS_2019' 
         & study_condition %in% c('control', 'CRC')) %>% 
  mutate(study_name = 'YachidaS_2019_CRC',
         dis = 'CRC')

d33_2 <- smd %>% 
  filter(study_name == 'YachidaS_2019' 
         & study_condition %in% c('control', 'adenoma')) %>% 
  mutate(study_name = 'YachidaS_2019_adenoma',
         dis = 'adenoma')

d34 <- smd %>% 
  filter(study_name == "YeZ_2018") %>% 
  mutate(dis = 'BD')

#CRC + T2D
d35 <- smd %>% 
  filter(study_name == "YuJ_2015") %>% 
  mutate(dis = 'CRC')

d36_1 <- smd %>% 
  filter(study_name == 'ZellerG_2014' 
         & study_condition %in% c('control', 'CRC')) %>% 
  mutate(study_name = 'ZellerG_2014_CRC',
         dis = 'CRC')

d36_2 <- smd %>% 
  filter(study_name == 'ZellerG_2014' 
         & study_condition %in% c('control', 'adenoma')) %>% 
  mutate(study_name = 'ZellerG_2014_adenoma',
         dis = 'adenoma')

d37 <- smd %>% 
  filter(study_name == "ZhuF_2020") %>% 
  mutate(dis = 'schizophrenia')


#Combine datasets---------------------------------------------------------------
#dl <- mget(ls()[1:46])

mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

rs <- curatedMetagenomicData::returnSamples
assays <- SummarizedExperiment::assay

meta_all <- bind_rows(dl) %>% 
  select(study = study_name, sample_id, study_condition, site = body_site, dis,
         reads = number_reads,
         bmi = BMI, age, sex = gender) %>% 
  mutate(type = 'Shotgun',
         group = ifelse(study_condition == 'control', 'control', 'case'),
         group = factor(group, levels = c('control', 'case'))) %>% 
  select(-study_condition) %>% 
  group_by(study) %>%
  mutate(n_control = sum(group == 'control'),
         n_case = sum(group == 'case'),
         n = n_control + n_case,
         sd_bmi = sd(bmi, na.rm = T),
         sd_bmi = ifelse(is.na(sd_bmi), 0, sd_bmi),
         sd_age = sd(age, na.rm = T),
         sd_age = ifelse(is.na(sd_age), 0, sd_age),
         n_sex = sum(!is.na(sex)),
         same_sex = sum(sex == 'male', na.rm = T) == n_sex | 
           sum(sex == 'female', na.rm = T) == n_sex,
         is_bmi = sum(!is.na(bmi)) / n >= .90 & sd_bmi > 0,
         is_age = sum(!is.na(age)) / n >= .90 & sd_age > 0,
         is_sex = sum(!is.na(sex)) / n >= .90 & !same_sex,
         is_sex = ifelse(str_detect(study, 'LiJ_2014'), F, is_sex), #Upd 260823
         half_elig = n_control >= 20 & n_case >= 20,
         avg_ls = exp(mean(log(reads)))) %>% 
  ungroup() %>% 
  group_by(study, group) %>% 
  mutate(bmi = ifelse(is.na(bmi) & is_bmi, median(bmi, na.rm = T), bmi),
         age = ifelse(is.na(age) & is_age, median(age, na.rm = T), age),
         sex = ifelse(is.na(sex) & is_sex, mode(sex), sex),
         sex = as.factor(sex)) %>% 
  ungroup() %>% 
  group_by(study) %>% 
  mutate_at(vars(bmi, age), ~ as.numeric(scale(.))) %>% 
  ungroup() %>% 
  select(type, site, disease = dis, avg_ls, study, half_elig,
         is_bmi, is_age, is_sex,
         n_control, n_case, n, id = sample_id, bmi, age, sex, group)

summary(meta_all)

tse <- smd %>%
  filter(study_name %in% studies$study_name) %>%
  filter(!study_name %in% c('VatanenT_2016', 'TettAJ_2016',
                            'MetaCardis_2020_a')) %>%
  rs("relative_abundance", counts = T, rownames = "short")

counts_all <- t(assays(tse)) %>% 
  as.data.frame() %>% 
  rename_all(~ paste0('genus_', 1:1363)) %>% 
  rownames_to_column('id')

d <- left_join(meta_all, counts_all, 'id')


#Half data----------------------------------------------------------------------
data_half <- list()
half_studies <- d %>% filter(half_elig) %>% distinct(study) %>% pull()
j <- 0

set.seed(1)
for(i in half_studies){
  data0 <- d %>% filter(study == i) %>% distinct(id, .keep_all = T)
  
  for(iter in 1:100){
    data1 <- data0 %>% 
      mutate(iter = iter) %>% 
      group_by(group) %>% 
      mutate(half = sample(rep(1:2, 10000)[1:n()])) %>% 
      ungroup()
    
    for(h in 1:2){
      data <- data1 %>% filter(half == h)
      
      meta <- data %>% 
        select(iter, half, type:group) %>% 
        column_to_rownames('id') %>% 
        mutate(data_id = paste('half', 'sg', study, iter, half)) %>% 
        as.data.frame()
      
      counts0 <- data %>% 
        select(id, matches('genus')) %>% 
        column_to_rownames('id') %>% 
        as.data.frame() %>% 
        as.matrix()
      
      counts1 <- counts0[, colMeans(!is.na(counts0)) == 1]
      
      counts <- counts1[, colMeans(counts1 > 0.5) >= .1]
      
      j <- j + 1
      data_half[[j]] <- list(meta = meta, counts = counts)
    }
  }
}

#save(data_half, file = 'data_sg_half_270923.rds')


#Whole_data---------------------------------------------------------------------
data_whole <- list()
whole_studies <- d %>% distinct(study) %>% pull()

for(i in 1:length(whole_studies)){
  data <- d %>% filter(study == whole_studies[i]) %>% 
    distinct(id, .keep_all = T)
  
  meta <- data %>% 
    select(type:group) %>% 
    column_to_rownames('id') %>% 
    mutate(data_id = paste('whole', 'sg', study)) %>%
    as.data.frame()
  
  counts0 <- data %>% 
    select(id, matches('genus')) %>% 
    column_to_rownames('id') %>% 
    as.data.frame() %>% 
    as.matrix()
  
  counts1 <- counts0[, colMeans(!is.na(counts0)) == 1]
  
  counts <- counts1[, colMeans(counts1 > 0.5) >= .1]
  
  data_whole[[i]] <- list(meta = meta, counts = counts)
}

#save(data_whole, file = 'data_sg_whole_041023.rds')

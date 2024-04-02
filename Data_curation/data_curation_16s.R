library(tidyverse)

#Read in and filter datasets----------------------------------------------------
studies <- tibble(folder = list.files(path = 'Data/Duvallet')) %>% 
  filter(str_detect(folder, 'results') & !str_detect(folder, 'tar\\.gz')) %>% 
  filter(!str_detect(folder, 'crc_zhu|ob_escobar')) %>% #No metadata
  mutate(folder = str_replace(folder, '_results', '')) %>% 
  pull(folder)

datasets <- list()
for(i in 1:length(studies)){
  counts_file <- paste0('Data/Duvallet/', studies[i], '_results/RDP/',
                        studies[i], '.otu_table.100.denovo.rdp_assigned')
  
  meta_file <- paste0('Data/Duvallet/', studies[i], '_results/', studies[i],
                      '.metadata.txt')
  
  counts0 <- read.table(file = counts_file, sep = '\t', header = T) %>% 
    tidytable::separate(X, into = c('kingdom', 'phylum', 'class', 'order',
                                    'family', 'genus', 'species', 'otu'),
                        sep = ';') %>% 
    group_by(kingdom, phylum, class, order, family, genus) %>% 
    summarize_if(is.numeric, ~ sum(.)) %>% 
    ungroup() %>% 
    mutate(taxon = paste(kingdom, phylum, class, order, family, genus,
                         sep = '-')) %>% 
    select(-c(kingdom:genus)) %>% 
    column_to_rownames('taxon') %>% 
    t() %>% 
    as.matrix()
  
  n_row_same <- F
  while(!n_row_same){
    counts1 <- counts0[, colMeans(counts0 > .5) >= .10]
    print(nrow(counts1))
    counts0 <- counts1[rowSums(counts1) > 200, ]
    print(nrow(counts0))
    n_row_same <- nrow(counts1) == nrow(counts0) 
  }
  
  counts <- counts0 %>% 
    as.data.frame() %>% 
    rownames_to_column('id') %>% 
    mutate(id = str_replace_all(id, '\\.', '-'),
           id = case_when(id == "Ademona1-2065" ~ "Adenoma1-2065",
                          id == "ADenoma12-2799" ~ "Adenoma12-2799",
                          id == "HG113" ~ "XHG113",
                          id == "HG237" ~ "XHG237",
                          id == "HG216" ~ "XHG216",
                          TRUE ~ id), 
           id = paste0(id, '_', studies[i]))
  
  rownames(counts) <- NULL
  
  header <- ifelse(studies[i] %in% c("autism_kb",
                                     "cdi_youngster",
                                     "crc_xiang",
                                     "crc_zhao",
                                     "hiv_noguerajulian",
                                     "ibd_engstrand_maxee",
                                     "ibd_huttenhower",
                                     "ra_littman",
                                     "t1d_mejialeon"), F, T)
  
  meta <- read.table(file = meta_file, sep = '\t', header = header,
                     fileEncoding = "latin1") %>% 
    rename(id = 1) %>% 
    mutate(study = studies[i],
           id = str_replace_all(as.character(id), '\\.', '-'),
           id = ifelse(study %in% c("crc_baxter", "hiv_dinh", "ibd_alm",
                                    "ibd_huttenhower", "ra_littman"),
                       paste0('X', id), id),
           id = ifelse(study == "ob_zupancic" & str_sub(id, 1, 1) != 'p',
                       paste0('X', id), id),
           id = paste0(id, '_', studies[i]))
  
  datasets[[i]] <- list(counts = counts, meta = meta)
  
  print(studies[i])
  print(c(dim(counts)[1], dim(meta)[1], dim(inner_join(meta, counts, 'id'))[1]))
}

#save(datasets, file = 'duv_uncurated_031023.rds')
load('duv_uncurated_031023.rds')


#Curate metadata----------------------------------------------------------------
m1 <- datasets[[1]]$meta %>% 
  mutate(disease = 'ASD',
         study = study,
         bmi = NA,
         age = NA,
         sex = factor(host_sex_s, levels = c('female', 'male')),
         group = factor(DiseaseState, levels = c('H', 'ASD'),
                        labels = c('control', 'case'))) %>% 
  select(id, disease, study, bmi, age, sex, group)

m2 <- datasets[[2]]$meta %>% 
  mutate(disease = 'ASD',
         study = study,
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(V6, levels = c('H', 'ASD'),
                        labels = c('control', 'case'))) %>% 
  select(id, disease, study, bmi, age, sex, group)

m3_1 <- datasets[[3]]$meta %>% 
  mutate(disease = 'CDI',
         study = paste0(study, '_cdi'),
         bmi = NA,
         age = age,
         sex = factor(gender, levels = c('F', 'M'),
                      labels = c('female', 'male')),
         group = factor(DiseaseState, levels = c('H', 'CDI'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m3_2 <- datasets[[3]]$meta %>% 
  mutate(disease = 'diarrhea',
         study = paste0(study, '_diarrhea'),
         bmi = NA,
         age = age,
         sex = factor(gender, levels = c('F', 'M'),
                      labels = c('female', 'male')),
         group = factor(DiseaseState, levels = c('H', 'nonCDI'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m4 <- datasets[[4]]$meta %>% 
  mutate(disease = 'CDI',
         study = study,
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(DiseaseState, levels = c('H', 'CDI'),
                        labels = c('control', 'case'))) %>% 
  select(id, disease, study, bmi, age, sex, group)

m5 <- datasets[[5]]$meta %>% 
  filter(V25 == 0 & V23 %in% c('H', 'CDI')) %>% 
  mutate(disease = 'CDI',
         study = study,
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(V23, levels = c('H', 'CDI'),
                        labels = c('control', 'case'))) %>% 
  select(id, disease, study, bmi, age, sex, group)

m6_1 <- datasets[[6]]$meta %>% 
  mutate(disease = 'adenoma',
         study = paste0(study, '_adenoma'),
         bmi = as.numeric(BMI_s) %>% if_else(. < 10, NA, ., NA),
         age = Age_s,
         sex = as.factor(Gender_s),
         group = factor(DiseaseState, levels = c('H', 'nonCRC'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m6_2 <- datasets[[6]]$meta %>% 
  mutate(disease = 'CRC',
         study = paste0(study, '_crc'),
         bmi = as.numeric(BMI_s) %>% if_else(. < 10, NA, ., NA),
         age = Age_s,
         sex = as.factor(Gender_s),
         group = factor(DiseaseState, levels = c('H', 'CRC'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m7 <- datasets[[7]]$meta %>% 
  mutate(disease = 'CRC',
         study = study,
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(V3, levels = c('H', 'CRC'),
                        labels = c('control', 'case'))) %>% 
  select(id, disease, study, bmi, age, sex, group)

m8_1 <- datasets[[8]]$meta %>% 
  mutate(disease = 'adenoma',
         study = paste0(study, '_adenoma'),
         bmi = NA,
         age = as.numeric(str_sub(age, 1, 2)),
         sex = factor(str_sub(sex, 1, 1), levels = c('f', 'm'),
                      labels = c('female', 'male')),
         group = factor(DiseaseState, levels = c('H', 'nonCRC'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m8_2 <- datasets[[8]]$meta %>% 
  mutate(disease = 'CRC',
         study = paste0(study, '_CRC'),
         bmi = NA,
         age = as.numeric(str_sub(age, 1, 2)),
         sex = factor(str_sub(sex, 1, 1), levels = c('f', 'm'),
                      labels = c('female', 'male')),
         group = factor(DiseaseState, levels = c('H', 'CRC'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m9 <- datasets[[9]]$meta %>% 
  mutate(disease = 'CRC',
         study = study,
         bmi = BMI..kg.m..,
         age = Age..years.,
         sex = factor(Gender, levels = c('F', 'M'),
                      labels = c('female', 'male')),
         group = factor(DiseaseState, levels = c('H', 'CRC'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

#Very low library sizes!
m10 <- datasets[[10]]$meta %>% 
  mutate(disease = 'CRC',
         study = study,
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(V2, levels = c('H', 'CRC'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m11 <- datasets[[11]]$meta %>% 
  mutate(disease = 'EDD',
         study = study,
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(DiseaseState, levels = c('H', 'EDD'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m12 <- datasets[[12]]$meta %>% 
  mutate(disease = 'HIV',
         study = study,
         bmi = as.numeric(BMI_s),
         age = as.numeric(age_s),
         sex = as.factor(tolower(sex_s)),
         group = factor(DiseaseState, levels = c('H', 'HIV'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

#For samples with ls < 100 (after taxa filtering)
m13 <- datasets[[13]]$meta %>% 
  filter(time_point == 1) %>% 
  mutate(disease = 'HIV',
         study = study,
         bmi = body_mass_index,
         age = age,
         sex = NA,
         group = factor(DiseaseState, levels = c('H', 'HIV'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

#A few low readcounts (60, 102, 150, 242...)
m14 <- datasets[[14]]$meta %>% 
  filter(V66 %in% c('BCN0', 'STK')) %>% 
  mutate(disease = 'HIV',
         study = study,
         bmi = NA,
         age = NA,
         sex = as.factor(V34),
         group = factor(V64, levels = c('H', 'HIV'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m15 <- datasets[[15]]$meta %>% 
  mutate(disease = 'IBD',
         study = study,
         bmi = NA,
         age = age,
         sex = as.factor(gender),
         group = case_when(DiseaseState == 'nonIBD' ~ 'control',
                           DiseaseState %in% c('UC', 'CD') ~ 'case',
                           T ~ NA),
         group = factor(group, levels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

#Smallest ls: 68  90 130 143 145 191
m16 <- datasets[[16]]$meta %>% 
  filter(V4 == 'stool') %>% 
  mutate(disease = 'IBD',
         study = study,
         bmi = NA,
         age = V12,
         sex = factor(V10, levels = c('F', 'M'),
                      labels = c('female', 'male')),
         group = case_when(V6 == 'H' ~ 'control',
                           V6 %in% c('UC', 'CD') ~ 'case',
                           T ~ NA),
         group = factor(group, levels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m17 <- datasets[[17]]$meta %>% 
  mutate(disease = 'IBD',
         study = study,
         bmi = NA,
         age = NA,
         sex = NA,
         group = case_when(DiseaseState == 'nonIBD' ~ 'control',
                           DiseaseState %in% c('UC', 'CD') ~ 'case',
                           T ~ NA),
         group = factor(group, levels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m18 <- datasets[[18]]$meta %>% 
  filter(V3 == 'stool') %>% 
  mutate(disease = 'IBD',
         study = study,
         bmi = NA,
         age = NA,
         sex = factor(V15, levels = c('F', 'M'),
                      labels = c('female', 'male')),
         group = case_when(V5 == 'H' ~ 'control',
                           V5 %in% c('UC', 'CD') ~ 'case',
                           T ~ NA),
         group = factor(group, levels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m19_1 <- datasets[[19]]$meta %>% 
  mutate(disease = 'CIRR',
         study = paste0(study, '_CIRR'),
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(DiseaseState, levels = c('H', 'CIRR'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m19_2 <- datasets[[19]]$meta %>% 
  mutate(disease = 'MHE',
         study = paste0(study, '_MHE'),
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(DiseaseState, levels = c('H', 'MHE'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m20 <- datasets[[20]]$meta %>% 
  filter(Status_s == 'Baseline') %>% 
  mutate(disease = 'NASH',
         study = study,
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(DiseaseState, levels = c('H', 'NASH'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m21_1 <- datasets[[21]]$meta %>% 
  mutate(disease = 'NASH',
         study = paste0(study, '_NASH'),
         bmi = NA,
         age = age,
         sex = factor(sex, levels = c('F', 'M'),
                      labels = c('female', 'male')),
         group = factor(DiseaseState, levels = c('H', 'NASH'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m21_2 <- datasets[[21]]$meta %>% 
  mutate(disease = 'OB',
         study = paste0(study, '_OB'),
         bmi = NA,
         age = age,
         sex = factor(sex, levels = c('F', 'M'),
                      labels = c('female', 'male')),
         group = factor(DiseaseState, levels = c('H', 'nonNASH-OB'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m22_1 <- datasets[[22]]$meta %>% 
  filter(n_sample == 0) %>% 
  mutate(disease = 'OB',
         study = paste0(study, '_OB'),
         bmi = body_mass_index,
         age = age,
         sex = NA,
         group = factor(DiseaseState, levels = c('H', 'OB'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m22_2 <- datasets[[22]]$meta %>% 
  filter(n_sample == 0) %>% 
  mutate(disease = 'OW',
         study = paste0(study, '_OW'),
         bmi = body_mass_index,
         age = age,
         sex = NA,
         group = factor(DiseaseState, levels = c('H', 'OW'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m23_1 <- datasets[[23]]$meta %>% 
  mutate(disease = 'OB',
         study = paste0(study, '_OB'),
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(DiseaseState, levels = c('H', 'OB'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m23_2 <- datasets[[23]]$meta %>% 
  mutate(disease = 'OW',
         study = paste0(study, '_OW'),
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(DiseaseState, levels = c('H', 'OW'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m24_1 <- datasets[[24]]$meta %>% 
  mutate(disease = 'OB',
         study = paste0(study, '_OB'),
         bmi = BMI,
         age = age_at_visit,
         sex = factor(sex, levels = c('F', 'M'),
                      labels = c('female', 'male')),
         group = factor(DetailedDiseaseState, levels = c('H', 'OB'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m24_2 <- datasets[[24]]$meta %>% 
  mutate(disease = 'OW',
         study = paste0(study, '_OW'),
         bmi = BMI,
         age = age_at_visit,
         sex = factor(sex, levels = c('F', 'M'),
                      labels = c('female', 'male')),
         group = factor(DetailedDiseaseState, levels = c('H', 'OW'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m25 <- datasets[[25]]$meta %>% 
  group_by(submitted_subject_id_s) %>% 
  arrange(visit_number) %>% 
  filter(row_number() == 1) %>% 
  ungroup() %>% 
  mutate(disease = 'OB',
         study = study,
         bmi = NA,
         age = NA,
         sex = as.factor(sex_s),
         group = factor(DiseaseState, levels = c('H', 'OB'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m26 <- datasets[[26]]$meta %>% 
  mutate(disease = 'parkinson',
         study = study,
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(DiseaseState, levels = c('H', 'PAR'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m27_1 <- datasets[[27]]$meta %>% 
  mutate(disease = 'NORA',
         study = paste0(study, '_NORA'),
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(V6,
                        levels = c('Healthy', 'New onset rheumatoid arthritis'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m27_2 <- datasets[[27]]$meta %>% 
  mutate(disease = 'CRA',
         study = paste0(study, '_CRA'),
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(V6,
                        levels = c('Healthy',
                                   'Treated chronic rheumatoid arthritis'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m27_3 <- datasets[[27]]$meta %>% 
  mutate(disease = 'PSA',
         study = paste0(study, '_PSA'),
         bmi = NA,
         age = NA,
         sex = NA,
         group = factor(V6,
                        levels = c('Healthy',
                                   'Psoriatic arthritis'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m28 <- datasets[[28]]$meta %>% 
  mutate(disease = 'T1D',
         study = study,
         bmi = NA,
         age = Age,
         sex = factor(Gender, levels = c('F', 'M'),
                      labels = c('female', 'male')),
         group = factor(DiseaseState, levels = c('H', 'T1D'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)

m29 <- datasets[[29]]$meta %>% 
  mutate(disease = 'T1D',
         study = study,
         bmi = NA,
         age = NA,
         sex = factor(V19, levels = c('f', 'm'),
                      labels = c('female', 'male')),
         group = factor(V52, levels = c('H', 'T1D'),
                        labels = c('control', 'case'))) %>% 
  filter(!is.na(group)) %>% 
  select(id, disease, study, bmi, age, sex, group)


#Combine datasets---------------------------------------------------------------
mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

metas <- bind_rows(m1, m10, m11, m12, m13, m14, m15, m16, m17, m18, m19_1,
                   m19_2, m2, m20, m21_1, m21_2, m22_1, m22_2, m23_1, m23_2,
                   m24_1, m24_2, m25, m26, m27_1, m27_2, m27_3, m28, m29, m3_1,
                   m3_2, m4, m5, m6_1, m6_2, m7, m8_1, m8_2, m9)

counts_all <- bind_rows(map(datasets, ~ .$counts))
colnames(counts_all)[2:295] <- paste0('genus_', 1:294)

meta_all <- inner_join(metas, counts_all, 'id') %>%
  mutate(reads = rowSums(select(., matches('genus')), na.rm = T),
         type = '16S',
         site = 'stool',
         sex = as.character(sex)) %>% 
  select(study, type, site, id, group, disease, reads, bmi, age, sex) %>% 
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
         half_elig = n_control >= 20 & n_case >= 20,
         avg_ls = exp(mean(log(reads)))) %>% 
  ungroup() %>% 
  group_by(study, group) %>% 
  mutate(bmi = ifelse(is.na(bmi) & is_bmi, median(bmi, na.rm = T), bmi),
         age = ifelse(is.na(age) & is_age, median(age, na.rm = T), age),
         sex = ifelse(is.na(sex) & is_sex, mode(sex), sex),
         sex = as.factor(sex)
         ) %>% 
  ungroup() %>% 
  group_by(study) %>% 
  mutate_at(vars(bmi, age), ~ as.numeric(scale(.))) %>% 
  ungroup() %>% 
  select(type, site, disease, avg_ls, study, half_elig,
         is_bmi, is_age, is_sex,
         n_control, n_case, n, id, bmi, age, sex, group)

summary(meta_all)

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
        mutate(data_id = paste('half', type, study, iter, half)) %>% 
        as.data.frame()
      
      counts0 <- data %>% 
        select(id, matches('genus')) %>% 
        column_to_rownames('id') %>% 
        as.data.frame() %>% 
        as.matrix()
      
      counts1 <- counts0[, colSums(!is.na(counts0)) > 1.5]
      
      counts <- counts1[, colMeans(counts1 > 0.5) >= .10]
      
      j <- j + 1
      data_half[[j]] <- list(meta = meta, counts = counts)
    }
  }
}

save(data_half, file = 'data_16s_half_031023.rds')


#Whole_data---------------------------------------------------------------------
data_whole <- list()
whole_studies <- d %>% distinct(study) %>% pull()

for(i in 1:length(whole_studies)){
  data <- d %>% filter(study == whole_studies[i]) %>% 
    distinct(id, .keep_all = T)
  
  meta <- data %>% 
    select(type:group) %>% 
    column_to_rownames('id') %>% 
    mutate(data_id = paste('whole', type, study)) %>%
    as.data.frame()
  
  counts0 <- data %>% 
    select(id, matches('genus')) %>% 
    column_to_rownames('id') %>% 
    as.data.frame() %>% 
    as.matrix()
  
  counts1 <- counts0[, colSums(!is.na(counts0)) > 1.5]
  
  counts <- counts1[, colMeans(counts1 > 0.5) >= .10]
  
  data_whole[[i]] <- list(meta = meta, counts = counts)
}

save(data_whole, file = 'data_16s_whole_031023.rds')
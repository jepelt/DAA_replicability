library(tidyverse)

#Helper functions---------------------------------------------------------------
make_physeq <- function(counts, meta){
  otu <- phyloseq::otu_table(t(counts), taxa_are_rows = TRUE)
  meta_data <- phyloseq::sample_data(meta)
  phyloseq::phyloseq(otu, meta_data)
}


#ALDEx2-------------------------------------------------------------------------
time_aldex <- function(counts, meta, fm = ~ group){
  mm <- model.matrix(fm, meta)
  clrs <- ALDEx2::aldex.clr(t(counts), mm)
  obj_aldex <- ALDEx2::aldex.glm(clrs, mm)
}


#ANCOM-BC2----------------------------------------------------------------------
time_ancombc2 <- function(counts, meta, fm = ~ group){
  physeq <- make_physeq(counts, meta)
  
  obj <- ANCOMBC::ancombc2(data = physeq, fix_formula = 'group',
                           verbose = F)
}


#corncob------------------------------------------------------------------------
time_corncob <- function(counts, meta, fm = ~ group, ev = T){
  nvar <- length(all.vars(fm))
  f_phi <- ifelse(ev, '~ 1', '~ group')
  fm0 <- if(nvar > 1){update.formula(fm, ~ . - group)}else{ ~ 1}
  
  obj <- corncob::differentialTest(data = t(counts),
                                   sample_data = meta,
                                   formula = fm,
                                   formula_null = fm0,
                                   phi.formula = formula(f_phi),
                                   phi.formula_null = formula(f_phi),
                                   test = "LRT")
}


#DESeq2-------------------------------------------------------------------------
time_deseq <- function(counts, meta, fm = ~ group, norm = 'Default'){
  m <- t(counts)
  mode(m) <- "integer"
  
  deseq_data <- DESeq2::DESeqDataSetFromMatrix(countData = m,
                                               colData = meta,
                                               design = fm)
  
  obj <- DESeq2::DESeq(deseq_data, sfType = 'poscounts', quiet = T)
}


#edgeR--------------------------------------------------------------------------
time_edger <- function(counts, meta, fm = ~ group, norm = 'TMM', test = 'QL'){
  mm <- model.matrix(fm, meta)
  
  dge <- edgeR::DGEList(t(counts))
  normalized <- edgeR::calcNormFactors(dge, method = norm)
  dispersions <- edgeR::estimateDisp(normalized, mm)
  
  fit_ql <- edgeR::glmQLFit(dispersions, mm)
  obj <- edgeR::glmQLFTest(fit_ql, coef = 2)
}


#fastANCOM----------------------------------------------------------------------
time_fastancom <- function(counts, meta, fm = ~ group, pseudo = 0.5, sig = 0.05,
                           ref.rate = 0.05){
  
  obj <- fastANCOM::fastANCOM(Y = counts, x = meta$group, Z = NULL,
                              zero_cut = 1,
                              pseudo = pseudo,
                              sig = sig,
                              ref.rate = ref.rate)
}


#LDM----------------------------------------------------------------------------
time_ldm <- function(counts, meta, fm = ~ group, clr = F){
  
  obj <- LDM::ldm(counts ~ group, data = meta, comp.anal = clr,
                  verbose = F, n.cores = 1)
}


#limma-voom---------------------------------------------------------------------
time_limma <- function(counts, meta, fm = ~ group, norm = 'TMM'){
  mm <- model.matrix(fm, meta)
  
  dge <- edgeR::DGEList(t(counts))
  normalized <- edgeR::calcNormFactors(dge, method = norm)
  obj_voom <- limma::voom(normalized, mm, plot = F)
  lm_fit <- limma::lmFit(obj_voom, mm)
  obj <- limma::eBayes(lm_fit)
}


#LinDA--------------------------------------------------------------------------
time_linda <- function(counts, meta, fm = ~ group, zeroh = 'Adaptive', pc = 0.5){
  
  obj <- LinDA::linda(otu.tab = t(counts), meta = meta,
                      formula = '~ group',
                      adaptive = T, imputation = F,
                      pseudo.cnt = pc)
}


#Logistic regression (Firth)----------------------------------------------------
time_logr_firth <- function(counts, meta, fm = ~ group){
  presences <- 1 * (counts > 0.5)
  meta_glm <- meta %>% select(all.vars(fm))
  
  obj <- apply(presences, 2, da_logr_firth, data = meta_glm)
}


da_logr_firth <- function(y, data){
  data$y <- y
  
  #Uncommenting the line below may cause paste() function work improperly!
  m <- logistf::logistf(y ~ ., data = data, plconf = 2, alpha = .166)
  
  if(sum(y) < length(y)){
    c(est = m$coefficients['groupcase'],
      p = m$prob['groupcase'],
      lwr = m$ci.lower['groupcase'],
      upr = m$ci.upper['groupcase'])
  }else{
    c(0, 1, NA, NA)
  }
}


#MaAsLin2-----------------------------------------------------------------------
time_maaslin2 <- function(counts, meta, fm = ~ group, tr = 'LOG', norm = 'TSS'){
  
  obj <- Maaslin2::Maaslin2(input_data = t(counts), input_metadata = meta,
                            output = 'output', min_prevalence = 0,
                            normalization = norm, transform = tr,
                            fixed_effects = all.vars(fm), standardize = F,
                            plot_heatmap = F, plot_scatter = F)
}


#metagenomeSeq------------------------------------------------------------------
time_mgs <- function(counts, meta, norm = 'CSS'){
  mgs_data <- metagenomeSeq::newMRexperiment(counts = t(counts),
                                             phenoData = afr(meta))
  normalized <- metagenomeSeq::cumNorm(mgs_data, p = 0.5)
  
  obj_mgs <- metagenomeSeq::fitFeatureModel(normalized,
                                            model.matrix(~ 1 + group,
                                                         data = meta))
}

afr <- Biobase::AnnotatedDataFrame


#Ordinal regression model-------------------------------------------------------
time_orm <- function(counts, meta, fm = ~ group, norm = 'TSS'){
  m <- t(counts)
  normalized <- counts / colSums(m)
  
  mf <- meta %>% transmute(group = ifelse(meta$group == 'control', 0, 1)) %>% 
    as.matrix
  
  res <- apply(normalized, 2, da_orm, mf = mf) %>%
    t() %>%
    as.data.frame()
}

da_orm <- function(y, mf){
  m1 <- rms::orm(as.numeric(y) ~ mf)
  
  if(ncol(mf) > 1){
    dv1 <- m1$deviance[2]
    score1 <- m1$stats['Score']
    
    m0 <- rms::orm(as.numeric(y) ~ mf[, -1])
    dv0 <- m0$deviance[2]
    score0 <- m0$stats['Score']
    
    p <- ifelse(!m0$fail, 1 - pchisq((dv0 - dv1), df = 1), NA)
    p_score <- ifelse(!m0$fail, 1 - pchisq((score1 - score0), df = 1), NA)
  }else{
    p <- ifelse(!m1$fail, m1$stats['P'], NA)
    p_score <- ifelse(!m1$fail, m1$stats['Score P'], NA)
  }
  
  if(!m1$fail){
    values <- c(m1$coefficients['group'], sqrt(m1$var[2, 2]), p, p_score)
  }else{
    values <- rep(NA, 4)
  }
  
  return(as.numeric(values))
}


#radEmu-------------------------------------------------------------------------
time_rademu <- function(counts, meta, fm = ~ group){
  
  obj <- radEmu::emuFit(formula = fm, data = meta, Y = counts,
                        return_wald_p = T, run_score_tests = F)
}


#ZicoSeq------------------------------------------------------------------------
time_zicoseq <- function(counts, meta, fm = ~ group, version = 'Default'){
  
  obj <- GUniFrac::ZicoSeq(
    meta.dat = meta, feature.dat = t(counts),
    grp.name = 'group', adj.name = NULL,
    feature.dat.type = "count", perm.no = 999,
    verbose = F, return.feature.dat = T)
}


#Benchmark for runnig time------------------------------------------------------
load('data_time_300924.rds')

benchmarks <- benchmarks_logr <- list()
set.seed(1)

for(b in 1:length(d_time)){
  counts <- d_time[[b]]$counts
  meta <- d_time[[b]]$meta
  
  benchmark <- microbenchmark::microbenchmark(
    time_aldex(counts, meta),
    time_ancombc2(counts, meta),
    time_corncob(counts, meta),
    time_deseq(counts, meta),
    time_edger(counts, meta),
    time_fastancom(counts, meta),
    time_ldm(counts, meta),
    time_limma(counts, meta),
    time_linda(counts, meta),
    time_maaslin2(counts, meta),
    time_mgs(counts, meta),
    time_orm(counts, meta),
    time_zicoseq(counts, meta),
    times = 1
  )
  
  benchmarks[[b]] <- benchmark %>% 
    as.data.frame() %>% 
    mutate(type = meta$type[1],
           data_id = meta$data_id[1],
           n_control = meta$n_control[1],
           n_case = meta$n_case[1],
           n = n_control + n_case,
           avg_ls = meta$avg_ls[1],
           n_taxa = ncol(counts))
}

#This is run lastly as logistf() causes problems with paste()
for(b in 1:length(d_time)){
  counts <- d_time[[b]]$counts
  meta <- d_time[[b]]$meta
  
  benchmark <- microbenchmark::microbenchmark(
    time_logr_firth(counts, meta),
    times = 1
  )
  
  benchmarks_logr[[b]] <- benchmark %>% 
    as.data.frame() %>% 
    mutate(type = meta$type[1],
           data_id = meta$data_id[1],
           n_control = meta$n_control[1],
           n_case = meta$n_case[1],
           n = n_control + n_case,
           avg_ls = meta$avg_ls[1],
           n_taxa = ncol(counts))
}

benchmark_res <- bind_rows(benchmarks, benchmarks_logr)
#save(benchmark_res, file = 'res_benchmark_241024.rds')

#Calculate typical runnig times and create Figure A..---------------------------
#load('res_benchmark_241024.rds')

run_times <- benchmark_res %>% 
  mutate(m = gsub(".*_(.*?)\\(.*", "\\1", expr),
         method = recode(m, aldex = 'ALDEx2',
                         ancombc2 = 'ANCOM-BC2',
                         corncob = 'corncob',
                         deseq = 'DESeq2',
                         edger = 'edgeR',
                         fastancom = 'fastANCOM',
                         ldm = 'LDM',
                         limma = 'limma-voom',
                         linda = 'LinDA',
                         firth = 'LogR',
                         maaslin2 = 'MaAsLin2/t-test',
                         mgs = 'metagenomeSeq',
                         orm = 'ORM/Wilcoxon',
                         zicoseq = 'ZicoSeq'),
         time = time / 10 ^ 9,
         dataset = as.factor(paste0('N = ', n, '; N_taxa = ', n_taxa)),
         dataset = fct_reorder(dataset, n * n_taxa)) %>% 
  group_by(method) %>% 
  mutate(m = exp(mean(log(time)))) %>% 
  ungroup() %>% 
  mutate(method = fct_reorder(method, m))

run_times_gm <- run_times %>% 
  group_by(method) %>% 
  summarize(m = exp(mean(log(time))))

#writexl::write_xlsx(run_times_gm, 'table_times_241024.xlsx')

(p1 <- ggplot(run_times, aes(method, time, color = dataset)) + 
    geom_point(size = 2, alpha = .5) +
    geom_point(data = run_times_gm, aes(y = m), color = 'black', size = 7,
               shape = 124) +
    scale_y_log10(breaks = c(.01, .1, 1, 10, 100),
                  labels = c('0.01', '0.1', '1', '10', '100'),
                  name = 'Run-time (s)') +
    scale_x_discrete(limits = rev) +
    #labs(title = 'Figure A2: Run-times of 14 DAA methods on six datasets') +
    ggsci::scale_color_jco(name = 'Sample size and\nnumber of taxa\n in the dataset') +
    coord_flip() +
    theme_light() +
    theme(axis.title.y = element_blank(),
          axis.text.x = element_text(size = 8),
          axis.title.x = element_text(size = 10),
          legend.text = element_text(size = 8),
          plot.title = element_blank(),
          legend.position = 'bottom')
)

date <- format(Sys.Date(), "%d%m%y")
ggsave(filename = paste0('fig_A2_times_', date,'.png'),
       plot = p1, width = 170, height = 150, dpi = 300, unit = 'mm',
       bg = 'white')

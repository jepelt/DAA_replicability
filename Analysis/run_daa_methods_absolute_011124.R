library(tidyverse)

#Helper functions---------------------------------------------------------------
make_physeq <- function(counts, meta){
  otu <- phyloseq::otu_table(t(counts), taxa_are_rows = TRUE)
  meta_data <- phyloseq::sample_data(meta)
  phyloseq::phyloseq(otu, meta_data)
}


#ALDEx2-------------------------------------------------------------------------
run_aldex <- function(counts, meta, fm = ~ group, glm = T, gamma = NULL){
  if(!glm){
    g <- ifelse(meta$group == 'case', 1, 0)
    obj <- ALDEx2::aldex(reads = t(counts), conditions = g,
                         verbose = F, gamma = gamma, effect = T)
    
    res <- obj %>%
      rownames_to_column('taxon') %>%
      pivot_longer(cols = c(we.ep, wi.ep)) %>%
      mutate(method = ifelse(name == 'we.ep',
                             paste0('ALDEx2 (t-test', gamma, ')'),
                             paste0('ALDEx2 (Wilcoxon', gamma, ')'))
      ) %>% 
      select(taxon, est = effect, p = value, method)
    
  }else{
    mm <- model.matrix(fm, meta)
    clrs <- ALDEx2::aldex.clr(t(counts), mm)
    obj_aldex <- ALDEx2::aldex.glm(clrs, mm)
    
    res <- obj_aldex %>%
      rownames_to_column('taxon') %>%
      select(taxon,
             est = 'groupcase:Est',
             se = 'groupcase:SE',
             p = 'groupcase:pval') %>%
      mutate(df = nrow(meta) - length(all.vars(fm)) - 1,
             method = 'ALDEx2 (glm)')
  }
  
  return(res %>% mutate(data_id = meta$data_id[1]))
}


#ANCOM-BC2----------------------------------------------------------------------
run_ancombc2 <- function(counts, meta, fm = ~ group, p_adj = 'BH'){
  physeq <- make_physeq(counts, meta)
  fmc <- paste(as.character(fm)[-1], collapse = ' + ')
  
  obj_p <- ANCOMBC::ancombc2(data = physeq, fix_formula = fmc,
                             prv_cut = 0, p_adj_method = 'none',
                             verbose = F, alpha = .05)
  
  r <- obj_p$res
  res <- tibble(taxon = r$taxon,
                est = r$lfc_groupcase,
                se = r$se_groupcase,
                df = nrow(counts) - length(all.vars(fm)) - 1,
                p0 = r$p_groupcase,
                ss_pass_p = r$passed_ss_groupcase,
                method = paste0('ANCOM-BC2 (', p_adj, ')')) %>% 
    mutate(p = ifelse(p0 < 0.05 & !ss_pass_p, runif(length(p0), 0.05, 1), p0)) %>% 
    select(taxon, est, se, df, p, method, ss_pass_p)
  
  return(res %>% mutate(data_id = meta$data_id[1]))
}


#corncob------------------------------------------------------------------------
run_corncob <- function(counts, meta, fm = ~ group, ev = T){
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
  
  res_list <- list()
  for(i in 1:nrow(t(counts))){
    if(!is.logical(obj$all_models[[i]])){
      res_list[[i]] <- obj$all_models[[i]]$coefficients[2, ]
      df <- obj$all_models[[i]]$df.residual
    }else{
      res_list[[i]] <- rep(NA, 4)
      names(res_list[[i]]) <- c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)')
    }
  }
  
  res <- bind_rows(res_list) %>%
    dplyr::rename(est = Estimate, se = 'Std. Error') %>%
    mutate(taxon = rownames(t(counts)),
           df = Inf,
           p = 2 * (1 - pnorm(abs(est / se))),
           p_lrt = obj$p) %>%
    pivot_longer(cols = c(p, p_lrt)) %>%
    mutate(method = paste0('corncob (',
                           ifelse(name == 'p', 'Wald', 'LRT'), ', ev = ',
                           ev, ')')) %>%
    select(taxon, est, se, df, p = value, method)
  
  return(res %>% mutate(data_id = meta$data_id[1]))
}


#DESeq2-------------------------------------------------------------------------
run_deseq <- function(counts, meta, fm = ~ group, norm = 'Default'){
  m <- t(counts)
  mode(m) <- "integer"
  nvar <- length(all.vars(fm))
  fm0 <- if(nvar > 1){update.formula(fm, ~ . - group)}else{ ~ 1}
  
  deseq_data <- DESeq2::DESeqDataSetFromMatrix(countData = m,
                                               colData = meta,
                                               design = fm)
  
  if(norm == 'Default'){
    sf <- DESeq2::estimateSizeFactorsForMatrix(m, type = 'poscounts')}
  if(norm == 'TSS'){
    reads <- colSums(m)
    sf <- reads / exp(mean(log(reads)))}
  if(norm == 'Wrench'){
    cond <- meta$group
    sf <- Wrench::wrench(mat = m, condition = cond)$nf}
  if(norm == 'GMPR'){
    sf <- GUniFrac::GMPR(m)
    sf <- ifelse(is.na(sf), 1, sf)}
  
  DESeq2::sizeFactors(deseq_data) <- sf
  
  obj <- DESeq2::DESeq(deseq_data, sfType = 'poscounts', quiet = T)
  
  obj_lrt <- DESeq2::DESeq(obj, test = "LRT", reduced = fm0, quiet = T)
  
  res_wald <- DESeq2::results(obj, tidy = T, format = "DataFrame") %>%
    mutate(method = paste0('DESeq2 (Wald, ', norm, ')'),
           df = Inf) %>%
    select(taxon = row, est = log2FoldChange, se = lfcSE, df,
           p = pvalue, p_fdr = padj, method)
  
  res_lrt <- DESeq2::results(obj_lrt, tidy = T, format = "DataFrame") %>%
    mutate(method = paste0('DESeq2 (LRT, ', norm, ')'),
           df = Inf) %>%
    select(taxon = row, est = log2FoldChange, p = pvalue, p_fdr = padj, method)
  
  res <- bind_rows(res_wald, res_lrt)
  return(res %>% mutate(data_id = meta$data_id[1]))
}


#edgeR--------------------------------------------------------------------------
run_edger <- function(counts, meta, fm = ~ group, norm = 'TMM', test = 'QL'){
  mm <- model.matrix(fm, meta)
  
  dge <- edgeR::DGEList(t(counts))
  normalized <- edgeR::calcNormFactors(dge, method = norm)
  dispersions <- edgeR::estimateDisp(normalized, mm)
  
  if(test == 'QL'){
    fit_ql <- edgeR::glmQLFit(dispersions, mm)
    obj <- edgeR::glmQLFTest(fit_ql, coef = 2)
  }else if(test == 'LRT'){
    fit_glm <- edgeR::glmFit(dispersions, mm)
    obj <- edgeR::glmLRT(fit_glm, coef = 2)
  }
  
  res <- edgeR::topTags(obj, n = ncol(counts))$table %>%
    as.data.frame() %>%
    rownames_to_column('taxon') %>%
    mutate(method = paste0('edgeR (', test, ', ', norm,')')) %>%
    select(taxon, est = logFC, p = PValue, method)
  
  return(res %>% mutate(data_id = meta$data_id[1]))
}


#fastANCOM----------------------------------------------------------------------
run_fastancom <- function(counts, meta, fm = ~ group, pseudo = 0.5, sig = 0.05,
                          ref.rate = 0.05){
  covs <- NULL
  nvar <- length(all.vars(fm))
  if(nvar > 1){
    meta$sex <- as.numeric(meta$sex)
    covs <- meta %>% select(all_of(all.vars(fm)[2:nvar])) %>% as.matrix
  }
  
  obj <- fastANCOM::fastANCOM(Y = counts, x = meta$group, Z = covs,
                              zero_cut = 1,
                              pseudo = pseudo,
                              sig = sig,
                              ref.rate = ref.rate)
  
  res <- obj$results$final %>%
    rownames_to_column('taxon') %>%
    mutate(df = nrow(counts) - nvar - 1,
           method = paste0('fastANCOM (',
                           paste(pseudo, sig, ref.rate, sep = ', '), ')')) %>%
    select(taxon, est = log2FC, se = log2FC.SD, df,
           p = log2FC.pval, reject = REJECT, method)
  
  return(res %>% mutate(data_id = meta$data_id[1]))
}


#LDM----------------------------------------------------------------------------
run_ldm <- function(counts, meta, fm = ~ group, clr = F){
  covs <- paste(all.vars(fm)[-1], collapse = ' + ')
  
  if(covs == ""){
    obj <- LDM::ldm(counts ~ group, data = meta, comp.anal = clr,
                    verbose = F, n.cores = 1)
  }else{
    f <- paste0('counts | (', covs, ') ~ group') %>% as.formula
    obj <- LDM::ldm(f, data = meta, comp.anal = clr,
                    verbose = F, n.cores = 1)
  }
  
  stopifnot(rownames(obj$beta) %in% c('groupcase', 'groupcontrol'))
  res <- tibble(taxon = names(obj$q.otu.omni),
                name = rownames(obj$beta),
                dir = ifelse(name == 'groupcontrol', -1, 1),
                est_abs = obj$beta[1, ] %>% as.numeric,
                est = dir * est_abs,
                p = obj$p.otu.omni %>% as.numeric,
                p_fdr = obj$q.otu.omni %>% as.numeric,
                method = paste0('LDM (CLR = ', clr, ')')) %>%
    select(taxon, est, p, p_fdr, method)
  
  return(res %>% mutate(data_id = meta$data_id[1]))
}


#limma-voom---------------------------------------------------------------------
run_limma <- function(counts, meta, fm = ~ group, norm = 'TMM'){
  mm <- model.matrix(fm, meta)
  
  dge <- edgeR::DGEList(t(counts))
  normalized <- edgeR::calcNormFactors(dge, method = norm)
  obj_voom <- limma::voom(normalized, mm, plot = F)
  lm_fit <- limma::lmFit(obj_voom, mm)
  obj <- limma::eBayes(lm_fit)
  
  res <- limma::topTable(obj, coef = 2, n = ncol(counts), confint = .834) %>%
    rownames_to_column('taxon') %>%
    mutate(method = paste0('limma-voom (', norm, ')')) %>%
    select(taxon, est = logFC, p = P.Value, lwr = CI.L, upr = CI.R, method)
  
  return(res %>% mutate(data_id = meta$data_id[1]))
}


#LinDA--------------------------------------------------------------------------
run_linda <- function(counts, meta, fm = ~ group, zeroh = 'Adaptive', pc = 0.5){
  fmc <- paste('~', paste(as.character(fm)[-1], collapse = ' + '))
  stopifnot(zeroh %in% c('Adaptive', 'Imputation', 'Pseudo'))
  adapt <- zeroh == 'Adaptive'
  imp <- zeroh == 'Imputation'
  
  obj <- LinDA::linda(otu.tab = t(counts), meta = meta,
                      formula = fmc,
                      adaptive = adapt, imputation = imp,
                      pseudo.cnt = pc)
  
  res <- obj$output$groupcase %>%
    rownames_to_column('taxon') %>%
    mutate(method = paste0('LinDA (', zeroh, ', ', pc, ')')) %>%
    select(taxon, est = log2FoldChange, se = lfcSE, df, p = pvalue, method)
  
  return(res %>% mutate(data_id = meta$data_id[1]))
}


#Logistic regression (Firth)----------------------------------------------------
run_logr_firth <- function(counts, meta, fm = ~ group){
  presences <- 1 * (counts > 0.5)
  meta_glm <- meta %>% select(all.vars(fm))
  
  obj <- apply(presences, 2, da_logr_firth, data = meta_glm)
  
  res <- obj %>%
    t() %>%
    as.data.frame() %>%
    rename(est = 1, p = 2, lwr = 3, upr = 4) %>%
    rownames_to_column('taxon') %>%
    mutate(method = 'LogR (Firth)')
  
  return(res %>% mutate(data_id = meta$data_id[1]))
}


da_logr_firth <- function(y, data){
  data$y <- y
  
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
run_maaslin2 <- function(counts, meta, fm = ~ group, tr = 'LOG', norm = 'TSS'){
  
  obj <- Maaslin2::Maaslin2(input_data = t(counts), input_metadata = meta,
                            output = 'output', min_prevalence = 0,
                            normalization = norm, transform = tr,
                            fixed_effects = all.vars(fm), standardize = F,
                            plot_heatmap = F, plot_scatter = F)
  
  res <- obj$results %>%
    filter(metadata == 'group' & value == 'case') %>%
    mutate(df = N - length(all.vars(fm)) - 1,
           method = paste0('MaAsLin2 (', tr, ', ', norm, ')')) %>%
    select(taxon = feature, est = coef, se = stderr, df, p = pval, method)
  
  return(res %>% mutate(data_id = meta$data_id[1]))
}


#metagenomeSeq------------------------------------------------------------------
run_mgs <- function(counts, meta, norm = 'CSS'){
  mgs_data <- metagenomeSeq::newMRexperiment(counts = t(counts),
                                             phenoData = afr(meta))
  
  if(norm == 'CSS'){
    normalized <- metagenomeSeq::cumNorm(mgs_data, p = 0.5)}
  if(norm == 'Wrench'){
    normalized <- metagenomeSeq::wrenchNorm(mgs_data, meta$group)}
  if(norm == 'TSS'){
    normalized <- mgs_data
    metagenomeSeq::normFactors(normalized) <- rowSums(counts)}
  if(norm == 'GMPR'){
    normalized <- mgs_data
    sf <- GUniFrac::GMPR(t(counts))
    sf <- ifelse(is.na(sf), 1, sf)
    sf <- sf * rowSums(counts)
    metagenomeSeq::normFactors(normalized) <- sf
  }
  
  obj_mgs <- metagenomeSeq::fitFeatureModel(normalized,
                                            model.matrix(~ 1 + group,
                                                         data = meta))
  
  res <- metagenomeSeq::MRfulltable(obj_mgs, number = ncol(counts)) %>%
    rownames_to_column('taxon') %>%
    mutate(df = Inf,
           method = paste0('mSeq (', norm, ')')) %>%
    select(taxon, est = logFC, se, df, p = pvalues, method)
  
  return(res %>% mutate(data_id = meta$data_id[1]))
}

afr <- Biobase::AnnotatedDataFrame


#Negativebinomial regression----------------------------------------------------
run_nb <- function(counts, meta, fm = ~ group){
  d <- tibble(x = meta$group, libr = log(rowSums(counts)))
  
  res <- apply(counts, 2, da_nb, d = d) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column('taxon') %>%
    mutate(method = paste0('NB')) %>%
    select(taxon, est, se, df, p, method)
  
  return(res %>% mutate(data_id = meta$data_id[1]))
}

da_nb <- function(y, d){
  d$y <- y
  m <- glmmTMB::glmmTMB(y ~ x + offset(libr),
                        family = glmmTMB::nbinom2, data = d)
  
  cf <- coef(summary(m))$cond
  c(est = cf[2, 1], se = cf[2, 2], df = Inf, p = cf[2, 4])
}


#Ordinal regression model-------------------------------------------------------
run_orm <- function(counts, meta, fm = ~ group, norm = 'TSS'){
  m <- t(counts)
  meta2 <- meta %>%
    mutate(#sex = as.numeric(meta$sex),
      group = ifelse(meta$group == 'control', 0, 1))
  mf <- meta2 %>% select(all.vars(fm)) %>% as.matrix
  
  if(norm == 'TSS'){
    sf <- colSums(m)}
  if(norm == 'Wrench'){
    cond <- meta$group
    sf <- Wrench::wrench(mat = m, condition = cond)$nf}
  if(norm == 'GMPR'){
    sf <- GUniFrac::GMPR(m)
    sf <- ifelse(is.na(sf), 1, sf)}
  
  normalized <- counts / sf
  
  res <- apply(normalized, 2, da_orm, mf = mf) %>%
    t() %>%
    as.data.frame() %>%
    rename(est = 1, se = 2, p_lrt = 3, p_score = 4) %>%
    mutate(taxon = colnames(counts),
           df = Inf,
           p = 2 * (1 - pnorm(abs(est / se)))) %>%
    pivot_longer(cols = c(p, p_lrt, p_score)) %>%
    mutate(method = paste0('ORMW (',
                           case_when(name == 'p' ~ 'Wald, ',
                                     name == 'p_lrt' ~ 'LRT, ',
                                     name == 'p_score' ~ 'Score, '),
                           norm, ')')) %>%
    select(taxon, est, se, df, p = value, method)
  
  return(res %>% mutate(data_id = meta$data_id[1]))
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
run_rademu <- function(counts, meta, fm = ~ group){
  
  obj <- radEmu::emuFit(formula = fm, data = meta, Y = counts,
                        alpha = 1 - .834,
                        return_wald_p = T, run_score_tests = F)
  
  res <- obj$coef %>% 
    transmute(taxon = category,
              est = estimate,
              p = wald_p,
              lwr = lower,
              upr = upper,
              method = 'radEmu (Wald)')
  
  return(res %>% mutate(data_id = meta$data_id[1]))
}


#ZicoSeq------------------------------------------------------------------------
run_zicoseq <- function(counts, meta, fm = ~ group, version = 'Default'){
  nvar <- length(all.vars(fm))
  fm0 <- if(nvar > 1){update.formula(fm, ~ . - group)}else{NULL}
  covs <- if(nvar > 1){all.vars(fm0)}else{NULL}
  
  if(version == 'Default'){
    obj <- GUniFrac::ZicoSeq(
      meta.dat = meta, feature.dat = t(counts),
      grp.name = 'group', adj.name = covs,
      feature.dat.type = "count", perm.no = 999,
      verbose = F, return.feature.dat = T)
  }else if(version == 'Links'){
    obj <- GUniFrac::ZicoSeq(
      meta.dat = meta, feature.dat = t(counts), 
      grp.name = 'group', adj.name = covs,
      feature.dat.type = "count", perm.no = 999,
      link.func = list(function (x) x ^ 0.25,
                       function (x) x ^ 0.5,
                       function (x) x ^ 0.75), 
      verbose = F, return.feature.dat = T)
  }
  
  res <- data.frame(taxon = names(obj$p.raw),
                    est = obj$coef.list[[1]]['groupcase', ] %>% as.vector,
                    p = as.numeric(obj$p.raw),
                    p_fdr = as.numeric(obj$p.adj.fdr),
                    method = paste0('ZicoSeq (', version, ')'))
  
  return(res %>% mutate(data_id = meta$data_id[1]))
}


#Run methods--------------------------------------------------------------------
load('data_abs_200324.rds')

set.seed(1)
datasets <- data_abs[1]
lcounts <- map(datasets, ~ .$counts)
lmeta <- map(datasets, ~ .$meta)

funs <- list(
   ~ run_aldex(.x, .y, glm = T),
   ~ run_aldex(.x, .y, glm = F, gamma = 0.5),
   ~ run_ancombc2(.x, .y),
   ~ run_corncob(.x, .y, ev = F),
   ~ run_corncob(.x, .y, ev = T),
   ~ run_deseq(.x, .y, norm = 'Default'),
   ~ run_deseq(.x, .y, norm = 'TSS'),
   ~ run_deseq(.x, .y, norm = 'GMPR'),
   ~ run_deseq(.x, .y, norm = 'Wrench'),
   ~ run_edger(.x, .y, norm = 'TMM'),
   ~ run_edger(.x, .y, norm = 'TMMwsp'),
   ~ run_edger(.x, .y, norm = 'RLE'),
   ~ run_fastancom(.x, .y),
   ~ run_limma(.x, .y, norm = 'TMM'),
   ~ run_limma(.x, .y, norm = 'TMMwsp'),
   ~ run_limma(.x, .y, norm = 'RLE'),
   ~ run_linda(.x, .y),
   ~ run_mgs(.x, .y, norm = 'CSS'),
   ~ run_mgs(.x, .y, norm = 'TSS'),
   ~ run_mgs(.x, .y, norm = 'GMPR'),
   ~ run_mgs(.x, .y, norm = 'Wrench'),
   ~ run_nb(.x, .y),
   ~ run_orm(.x, .y, norm = 'TSS'),
   ~ run_orm(.x, .y, norm = 'GMPR'),
   ~ run_orm(.x, .y, norm = 'Wrench'),
   ~ run_rademu(.x, .y),
   ~ run_zicoseq(.x, .y)
)

resl <- list()
Sys.time()
for(k in 1:length(funs)){
  rs <- NULL
  try(rs <- map2(lcounts, lmeta, funs[[k]]))
  
  resl[[k]] <- rs 
  print(k)
  print(Sys.time())
}


res_ldm1 <- res_ldm2 <- res_maaslin1 <- res_maaslin2 <- res_maaslin3 <- 
  res_maaslin4 <- res_maaslin5 <- res_logr <- list()

set.seed(1)
for(i in 1:length(datasets)){
  counts <- lcounts[[i]]
  meta <- lmeta[[i]]
  
  res_ldm1[[i]] <- run_ldm(counts, meta)
  res_maaslin1[[i]] <- run_maaslin2(counts, meta)
  res_maaslin2[[i]] <- run_maaslin2(counts, meta, tr = 'AST')
  res_maaslin3[[i]] <- run_maaslin2(counts, meta, tr = 'NONE', norm = 'CLR')
  res_maaslin4[[i]] <- run_maaslin2(counts, meta, norm = 'CSS')
  res_maaslin5[[i]] <- run_maaslin2(counts, meta, norm = 'TMM')
  res_logr[[i]] <- run_logr_firth(counts, meta)
}

res_list <- c(resl, res_ldm1, res_ldm2, res_maaslin1, res_maaslin2, 
              res_maaslin3, res_maaslin4, res_maaslin5, res_logr)

#save(res_list, file = 'res_absolute_011124.rds')
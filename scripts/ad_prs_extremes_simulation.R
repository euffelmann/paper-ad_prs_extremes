
                    ##### step 0: load packages & define functions #####

rm(list = ls())
                    
# theory from Dudbridge see https://github.com/DudbridgeLab/avengeme, install with "devtools::install_github("DudbridgeLab/avengeme")"
suppressPackageStartupMessages(library(avengeme)) 
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse)) 

## transforms explained variance from observed to liability scale
prs_r2obs_to_r2liab <- function(K,P,prs_r2obs){
  ## Lee et al. 2012 Genet Epidemiology
  t = -qnorm(K,mean=0,sd=1) # disease threshold
  z<-dnorm(t)               # height of the normal distribution at T
  i1<-z/K                   # mean liability of A1 (eg Falconer and Mackay) 
  k1=i1*(i1-t)              # reduction in variance in A1
  i0<--z/(1-K)              # mean liability of A0
  k0=i0*(i0-t)              # reduction in variance in A0
  
  theta = i1*(P-K)/(1-K)*(i1*(P-K)/(1-K)-t) # theta in equation (15) Lee et al. 2012 Genet Epidemiology
  cv = K*(1-K)/z^2*K*(1-K)/(P*(1-P))        # C in equation (15)
  R2 = prs_r2obs*cv/(1+prs_r2obs*theta*cv)
  
  return(R2)
}

## transforms heritability from observed to liability scale
h2o_to_h2l<-function(K,P,h2o){
  ## Based on Eq 23 of Lee et al 2011 AJHG
  t = -qnorm(K,0,1) ; z = dnorm(t)
  return(h2o*K*K*(1-K)*(1-K)/{z*z*P*(1-P)})
}

                        ##### step 1: set parameters #####

## run interactively on local computer with parameters that speed up computation or (else) on LISA for full analysis
interactive <- F 
if (interactive == T) {
  
  plink <- "/usr/local/bin/plink"
  setwd("/Users/schneemil/Documents/CTG/projects/alzheimers/temp")
  run <- 1; r2l_prs <- 0.05; h2l <- 0.3; K <- 0.05; Pd <- Pt <- K 
  disorder_name <- "scz"
  nsnp <- 1200; pi0 <- (nsnp - 100) / nsnp; nsnp_causal <- round(nsnp*(1-pi0)) 
  threshold <- round(-qnorm(p = K, mean = 0, sd = 1), digits = 1)
  total_testdata_n <- 2000
  n_ind_per_run <- 1000
  min_maf <- 0.2
  prs_thresholds <- c(0.05, 0.5, 1)
  
} else {
  
  plink <- "/home/emilu/programs/plink1.9/plink1.9"
  # prs_r2obs_to_r2liab(K = 0.05, P = (379 / 761), prs_r2obs = 0.039) based on Jansen et al. (2019) excluding APOE
  r2l_prs <- 0.05 # variance explained my PRS on liability scale 
  h2l <- 0.1 # h2 of disease on liability scale
  K <- 0.05 # prevalence
  Pd <- K # proportion of cases in discovery cohort ~ as in Jansen et al. 
  Pt <- K # proportion of cases in target cohort 
  disorder_name <- "scz"
  nsnp <- 170000 ## number of SNPs to be simulated
  pi0 <- (nsnp - 1200) / nsnp # proportion of null SNPs with 1200 causal SNPs
  nsnp_causal <- round(nsnp*(1-pi0)) # number of causal SNPs: 1200 as estimated by Holland et al. (2020)
  threshold <- round(-qnorm(p = K, mean = 0, sd = 1), digits = 1)
  total_testdata_n <- 300000
  n_ind_per_run <- 1000
  min_maf <- 0.001
  prs_thresholds <- c(5e-8, 0.05, 0.5, 1)
}


                    ##### step 2: simulate individuals, makes binary plink files #####

# See Dudbridge 2013 Plos Genetics
# Calculates the size of discovery sample sample (GWAS sample) to achieve a given AUC, R2 or power with a PRS in the target sample.
temp.total_n <- sampleSizeForGeneScore( targetQuantity = "R2", ## see: https://github.com/DudbridgeLab/avengeme/blob/master/avengeme-manual.pdf
                                        targetValue = r2l_prs,
                                        nsnp = nsnp,
                                        vg1 = h2l,
                                        pi0 = pi0,
                                        binary = TRUE,
                                        prevalence = K,
                                        sampling = Pd)$n

dudbridge_p <- sampleSizeForGeneScore( targetQuantity = "R2", ## see: https://github.com/DudbridgeLab/avengeme/blob/master/avengeme-manual.pdf
                                        targetValue = r2l_prs,
                                        nsnp = nsnp,
                                        vg1 = h2l,
                                        pi0 = pi0,
                                        binary = TRUE,
                                        prevalence = K,
                                        sampling = Pd)$p

## I will simulate individuals in batches to be more memory-efficient
batches <- ceiling(temp.total_n / n_ind_per_run)

# simulate MAFs and snp_effects
MAFs <- runif(n = nsnp, min = min_maf, max = 0.5) ## small sample-sizes, have large MAF to prevent 0 allele counts
snp_effects <- rnorm(n = nsnp_causal, mean = 0, sd = sqrt(h2l / nsnp_causal)) # see derivation below

snp_effects <- c(snp_effects, rep(0, nsnp-nsnp_causal))
snp_names <- c(sprintf("snp_causal%05d",1:nsnp_causal),sprintf("snp_null%05d",(nsnp_causal+1):nsnp))
snp_names_allels <- paste( rep(snp_names,each=2),c("_A1","_A2"),sep="" )

simulate_individual <- function(){
  allele1 <- as.numeric(MAFs - runif(n=nsnp,min=0,max=1) > 0) # eu: simulating alleles corresponding to the frequencies.
  # the smaller the allele frequency, the less likely this will be > 1
  allele2 <- as.numeric(MAFs - runif(n=nsnp,min=0,max=1) > 0)
  
  # standardize all genotypes and multiply them by their SNP_effects
  # explanation: allele1+allele2 -> individual dosage; 2*MAFs -> mean; sqrt(2*MAFs*(1-MAFs)) -> standard deviation
  g_liab <- as.matrix(rbind({(allele1+allele2-2*MAFs)/sqrt(2*MAFs*(1-MAFs))})) %*% as.matrix(snp_effects)
  e_liab <- rnorm(n = 1, mean = 0, sd = sqrt(1 - h2l))
  liab <- g_liab + e_liab
  dis <- as.numeric(liab > threshold) + 1 # 1 = control, 2 = case
  snps <- c(rbind( (allele1+1),(allele2+1) ))
  
  names(g_liab) <- "g_liab"
  names(e_liab) <- "e_liab"
  names(liab) <- "liab"
  names(dis) <- disorder_name
  names(snps) <- snp_names_allels
  
  return(c(g_liab,e_liab,liab,dis,snps))
}


for (i in 1:batches) {
  print(paste0("batch: ",i, " out of ", batches))
  
  ## simulate discovery data
  # the last round of the loop needs fewer individuals, so loop here makes sure I simulate fewer
  if (!is.integer(temp.total_n / n_ind_per_run) & i == batches) {
    last_N <- round((temp.total_n - {batches - 1} * n_ind_per_run))
    N1 <- round(Pd * last_N) # cases
    N0 <- round((1-Pd) * last_N) # controls
  } else {
    N1 <- round(Pd * n_ind_per_run) # round(temp.total_n * Pd); N1 # cases
    N0 <- round((1-Pd) * n_ind_per_run) # round(temp.total_n * (1 - Pd)); N0 # controls
  }
  
  print(paste0("Simulate discovery data"))
  temp_discdata1 <- list()
  temp_discdata2 <- list()
  count <- 0
  while( length(temp_discdata1) + length(temp_discdata2) <  N1 + N0 ){ # eu: only store as many individuals as I need
    individual <- simulate_individual()
    # eu: here I save only those individuals I need. This allows me to save time, as I don't need so simulate a whole population first
    if( as.numeric(individual[disorder_name])==1 &&  length(temp_discdata1) < N0  ){ temp_discdata1[[length(temp_discdata1)+1]] <- individual }
    if( as.numeric(individual[disorder_name])==2 &&  length(temp_discdata2) < N1  ){ temp_discdata2[[length(temp_discdata2)+1]] <- individual }
  }
  
  ## collect discovery data
  temp_discdata1 <- cbind(phenotype = 1, do.call("rbind",temp_discdata1))
  temp_discdata2 <- cbind(phenotype = 2, do.call("rbind",temp_discdata2))
  discdata <-  as.data.frame(rbind(temp_discdata2,temp_discdata1))
  
  ## save discovery data as .ped file
  ##    - Family ID
  ##    - Individual ID
  ##    - Paternal ID
  ##    - Maternal ID
  ##    - Sex (1=male; 2=female; other=unknown)
  ##    - Phenotype (  -9 missing ; 0 missing ; 1 unaffected ; 2 affected)
  
  discdata$FID <- sprintf("FAM_DISC%06d", {1 + ({i - 1} * n_ind_per_run)}:{({i - 1} * n_ind_per_run) + {(N1 + N0)}})
  discdata$IID <- sprintf("IID_DISC%06d", {1 + ({i - 1} * n_ind_per_run)}:{({i - 1} * n_ind_per_run) + {(N1 + N0)}})
  discdata$pat <- discdata$mat <- 0; discdata$sex <- 2
  disc_pedfile <- discdata[, c("FID", "IID", "pat", "mat", "sex", "phenotype", snp_names_allels)]

  fwrite(as.data.frame(disc_pedfile),file = paste0("discdata_", i,".ped"), quote = F, sep= "\t" , col.names = F, row.names = F)
  
  ## save discovery data as .map file
  ##    - chromosome (1-22, X, Y or 0 if unplaced)
  ##    - rs# or snp identifier
  ##    - Genetic distance (morgans)
  ##    - Base-pair position (bp units)
  
  disc_mapfile <- as.data.frame(cbind(1,snp_names, 0, ((1:length(snp_names))) ))
  fwrite(disc_mapfile, file = paste0("discdata_", i,".map"), quote = F, sep= "\t", col.names = F, row.names = F)
  
  ## convert to binary files
  system(paste0(plink, " --file discdata_", i, " --make-bed  --out discdata_", i))
  system(paste0("rm discdata_", i,".ped discdata_", i,".map"))
  
  ## simulate test data
  print(paste0("Simulate test data"))
  if (n_ind_per_run * (i - 1) < total_testdata_n) {

    temp_testdata1 <- list()
    temp_testdata2 <- list()
    temp_N <- n_ind_per_run
    temp_N0 <- temp_N * (1 - Pt)
    temp_N1 <- temp_N * Pt
    while( length(temp_testdata1) + length(temp_testdata2) <  temp_N0 + temp_N1){ # eu: only store as many individuals as I need
      individual <- simulate_individual()
      # eu: here I save only those individuals I need. This allows me to save time, as I don't need so simulate a whole population first
      if( as.numeric(individual[disorder_name])==1 &&  length(temp_testdata1) < N0  ){ temp_testdata1[[length(temp_testdata1)+1]] <- individual }
      if( as.numeric(individual[disorder_name])==2 &&  length(temp_testdata2) < N1  ){ temp_testdata2[[length(temp_testdata2)+1]] <- individual }
    }
    
    ## collect testdata
    temp_testdata1 <- cbind(phenotype = 1, do.call("rbind",temp_testdata1))
    temp_testdata2 <- cbind(phenotype = 2, do.call("rbind",temp_testdata2))
    testdata <-  as.data.frame(rbind(temp_testdata2,temp_testdata1))
    
    ## save testdata as .ped file
    testdata$FID <- sprintf("FAM_TEST%06d", {1 + ({i - 1} * n_ind_per_run)}:{n_ind_per_run * i})
    testdata$IID <- sprintf("IID_TEST%06d", {1 + ({i - 1} * n_ind_per_run)}:{n_ind_per_run * i})
    testdata$pat <- testdata$mat <- 0; testdata$sex <- 2
    test_pedfile <- testdata[, c("FID", "IID", "pat", "mat", "sex", "phenotype", snp_names_allels)]
    fwrite(as.data.frame(test_pedfile),file=paste0("testdata_", i,".ped"), quote = F, sep= "\t" , col.names = F, row.names = F)
    
    ## save discovery data as .map file
    test_mapfile <- as.data.frame(cbind(1,snp_names, 0, ((1:length(snp_names))) ))
    fwrite(test_mapfile, file = paste0("testdata_", i,".map"), quote = F, sep= "\t", col.names = F, row.names = F)
    
    ## convert to binary files
    system(paste0(plink, " --file testdata_", i, " --make-bed  --out testdata_", i))
    system(paste0("rm testdata_", i,".ped testdata_", i,".map"))
  }
}


                      ##### step 3: merge individual plink files #####

for (sample in c("disc", "test")) {
  files <- unlist(fread(cmd = paste0("ls ", sample, "data_*.bim"), header=FALSE))
  names(files) <- NULL
  files <- gsub(files, pattern = ".bim", replacement = "")
  fileConn <- file(paste0(sample, "_files.txt")); writeLines(files, fileConn); close(fileConn)
  system(paste0(plink, " --merge-list ", sample, "_files.txt --make-bed --out ", sample, "data"))
  system(paste0("rm ", sample, "data_* ", sample, "_files.txt"))
}

if (as.numeric(str_extract(system("wc -l discdata.fam", intern = T), "\\-*\\d+\\.*\\d*")) != round(temp.total_n)) {print("Error: Missing individuals in disocvery cohort")}
if (as.numeric(str_extract(system("wc -l testdata.fam", intern = T), "\\-*\\d+\\.*\\d*")) != total_testdata_n) {print("Error: Missing individuals in testdata cohort")}
if (as.numeric(str_extract(system("wc -l discdata.bim", intern = T), "\\-*\\d+\\.*\\d*")) != nsnp) {print("Error: Missing SNPs in discovery cohort")}
if (as.numeric(str_extract(system("wc -l testdata.bim", intern = T), "\\-*\\d+\\.*\\d*")) != nsnp) {print("Error: Missing SNPs in testdata cohort")}

                      ##### step 4: run discovery GWAS #####

## list object to store all results
run_summary <- NULL
for (prs_threshold in prs_thresholds) {
  
  system( paste0(plink, " --bfile discdata --logistic --freq --out gwas_alz_sumstats") )
  temp_gwas_alz_sumstats <- fread("gwas_alz_sumstats.assoc.logistic")
  temp_gwas_alz_maf <- fread("gwas_alz_sumstats.frq")[, c("SNP", "MAF")]
  gwas_alz_sumstats <- merge(x = temp_gwas_alz_sumstats, y = temp_gwas_alz_maf, by = ("SNP"))
  gwas_alz_sumstats$LOG_OR <- log(gwas_alz_sumstats$OR)
  gwas_alz_sumstats$LOG_OR_SE <- abs(log(gwas_alz_sumstats$OR)/qnorm(gwas_alz_sumstats$P/2))
  gwas_alz_sumstats_prs_threshold <- gwas_alz_sumstats[gwas_alz_sumstats$P <= prs_threshold, c("SNP", "A1", "LOG_OR")]
  fwrite(x = gwas_alz_sumstats, file = "gwas_alz_sumstats.assoc.logistic", sep = "\t", row.names = F, col.names = T)
  fwrite(x = gwas_alz_sumstats_prs_threshold, file = "gwas_alz_sumstats_prs_threshold.txt", sep = "\t", row.names = F, col.names = F)
  
                      ##### step 5: create testdata files with varying degree of overlap with discovery data #####
  
  rm(list=ls()[grep("temp",ls())])
  
  ## create list of individuals to be extracted
  discdata.fam <- fread("discdata.fam")
  testdata.fam <- fread("testdata.fam")
  
  testdata_index50 <- sample(x = seq(1:dim(testdata.fam)[1]), replace = F)[1:(total_testdata_n/2)] # index for 50% overlap
  discdata_index50 <- sample(x = seq(1:dim(discdata.fam)[1]), replace = F)[1:(total_testdata_n/2)] # index for 50% overlap
  discdata_index100 <- sample(x = seq(1:dim(discdata.fam)[1]), replace = F)[1:(total_testdata_n)] # index for 100% overlap
  
  testdata50.fam <- testdata.fam[testdata_index50,]
  discdata50.fam <- discdata.fam[discdata_index50,]
  discdata100.fam <- discdata.fam[discdata_index100,]
  
  fwrite(x = testdata50.fam[,1:2], file = "testdata50.txt", quote = F, sep = "\t", row.names = F, col.names = F)
  fwrite(x = discdata50.fam[,1:2], file = "discdata50.txt", quote = F, sep = "\t", row.names = F, col.names = F)
  fwrite(x = discdata100.fam[,1:2], file = "discdata100.txt", quote = F, sep = "\t", row.names = F, col.names = F)
  
  ## extract individuals
  system(paste0(plink, " --bfile testdata --keep testdata50.txt --make-bed --out temp_testdata50"))
  system(paste0(plink, " --bfile discdata --keep discdata50.txt --make-bed --out temp_discdata50"))
  system(paste0(plink, " --bfile discdata --keep discdata100.txt --make-bed --out testdata100"))
  
  
  system(paste0(plink, " --bfile temp_testdata50 --bmerge temp_discdata50 --make-bed --out testdata50"))
  system("cp testdata.bim testdata0.bim"); system("cp testdata.fam testdata0.fam"); system("cp testdata.bed testdata0.bed") # testdata with 0% overlap
  system("rm temp*")
  
  
                      ##### step 6: ALZ GWAS model checks #####
  
  # saving explained variance (as model check)
  gwas_alz_sumstats_dudbridge_p <- gwas_alz_sumstats[gwas_alz_sumstats$P < dudbridge_p,]
  fwrite(x = gwas_alz_sumstats_dudbridge_p, file = "gwas_alz_sumstats_dudbridge_p.txt", sep = "\t", row.names = F, col.names = T)
  system(paste0(plink, " --bfile testdata0 --score gwas_alz_sumstats_dudbridge_p.txt 1 4 11 header sum --out r2l_prs"))
  prs <- fread("r2l_prs.profile")
  prs_r2obs <- summary(lm(data = prs, formula = PHENO ~ SCORESUM))$r.sq
  prs_r2liab <- prs_r2obs_to_r2liab(K = K, P = Pt, prs_r2obs = prs_r2obs)
  disc_fam <- fread("discdata.fam")
  K_disc <- table(disc_fam$V6)[2] / sum(table(disc_fam$V6))
  ncase <- sum(disc_fam$V6 == 2)
  ncontrol <- sum(disc_fam$V6 == 1)
  neff <- 4 / ((1 / ncase) + (1 / ncontrol))
  ldsc_intercept <- var(gwas_alz_sumstats$STAT[(nsnp_causal+1):nsnp])
  ldsc_h2o <- nsnp * (mean(gwas_alz_sumstats$STAT^2) - 1) / neff# eu. suppl. table 57 "LD Score regression distinguishes confounding from polygenicity in genome-wide association studies"
  ldsc_h2l <- h2o_to_h2l(K=K,P=0.5,h2o=ldsc_h2o) 
  lambda <- median(gwas_alz_sumstats$STAT^2) / qchisq(0.5, 1)
  
  gwas_alz_sumstats_summary <- data.frame(prs_r2liab = round(prs_r2liab, 3),
                                          ldsc_intercept = round(ldsc_intercept, 2),
                                          ldsc_h2l = round(ldsc_h2l, 2),
                                          K = round(K_disc, 2),
                                          lambda = round(lambda, 2),
                                          neff = round(neff, 2),
                                          ncase = round(ncase, 2),
                                          ncontrol = round(ncontrol, 2))
  
  
                      ##### step 7: Polygenic Risk Scoring #####
  
  for (overlap in c("0", "50", "100")) { # percentage overlap
    system(paste0(plink, " --bfile testdata", overlap, " --score gwas_alz_sumstats_prs_threshold.txt 1 2 3 sum --out testdata_prs", overlap))
    
  }
  
  
                      ##### step 8: run GWAS on PRS extremes #####
  
  for (overlap in c("0", "50", "100")) { 
    # load alzheimers PRS
    testdata_prs <- fread(paste0("testdata_prs", overlap, ".profile"))
    # define PRS extremes
    quants <- quantile(testdata_prs$SCORESUM,probs = seq(0, 1, 0.05))
    prs_extremes <- testdata_prs[testdata_prs$SCORESUM < quants[2] | testdata_prs$SCORESUM > quants[20],]
    prs_extremes$phenotype_new[prs_extremes$SCORESUM > quants[20]] <- 2 # top of PRS
    prs_extremes$phenotype_new[prs_extremes$SCORESUM < quants[2]] <- 1 # bottom of PRS
    # write list of individuals that fall within the top and bottom 5%
    fwrite(x = prs_extremes[,1:2], file = paste0("prs_extremes", overlap, ".txt"), quote = F, sep = "\t", row.names = F, col.names = F)
    # write new plink files with selected individuals
    system(paste0(plink, " --bfile testdata", overlap, " --keep prs_extremes", overlap, ".txt --make-bed --out prs_extremes", overlap))
    # create new .fam file with dichotomous PRS extremes phenotype
    fam <- fread(paste0("prs_extremes", overlap, ".fam"))
    prs_extremes_fam <- merge(x = fam, y = prs_extremes[, c("FID", "IID", "phenotype_new")], by.x = c("V1", "V2") , by.y = c("FID", "IID"))
    prs_extremes_fam$V6 <- NULL
    fwrite(x = prs_extremes_fam, file = paste0("prs_extremes", overlap, ".fam"), quote = F, sep = "\t", row.names = F, col.names = F)
    # run GWAS on new PRS extremes phenotype
    system(paste0(plink, " --bfile prs_extremes", overlap," --logistic --freq --out gwas_prs_extremes", overlap, "_sumstats"))
  }
  
  
                      ##### step 9: testdata key statistics #####
  
  sumstats0 <- fread("gwas_prs_extremes0_sumstats.assoc.logistic")
  null_snps_in_prs <- gwas_alz_sumstats_prs_threshold %>% filter(str_detect(SNP, "null")) %>% .$SNP # all null SNPs that are in the PRS after applying the threshold
  null_snps_not_in_prs <- sumstats0 %>% filter(str_detect(SNP, "null")) %>% filter(!SNP %in% null_snps_in_prs) %>% .$SNP # all null SNPs that are NOT in the PRS
  all_null_snps <- sumstats0[(nsnp_causal+1):nsnp]$SNP # all null SNPs used in the GWAS
  
  fpr <- NULL
  prs_extremes_summary <- NULL
  z_stat <- NULL
  for (overlap in c("0", "50", "100")) { 
    
    sumstats <- fread(paste0("gwas_prs_extremes", overlap, "_sumstats.assoc.logistic"))
    
    # cases in PRS risk extremes
    cases_controls <- fread(paste0("testdata", overlap, ".fam"))
    prs_extremes <- fread(paste0("prs_extremes", overlap, ".fam"))
    n_cases_in_prs_top <- sum(prs_extremes$V2[prs_extremes$V6 == 2] %in% cases_controls$V2[cases_controls$V6 == 2])
    n_cases_in_prs_bottom <- sum(prs_extremes$V2[prs_extremes$V6 == 1] %in% cases_controls$V2[cases_controls$V6 == 2])  
    
    # N
    n_prs_top <- sum(prs_extremes$V6 == 2)
    n_prs_bottom <- sum(prs_extremes$V6 == 1)
    neff <- 4 / ((1 / n_prs_top) + (1 / n_prs_bottom))
    
    # ldsc h2 and genomic inflation factor
    ldsc_h2o <- nsnp * (mean(sumstats$STAT^2) - 1) / neff # eu. suppl. table 57 "LD Score regression distinguishes confounding from polygenicity in genome-wide association studies"
    ldsc_h2l <- h2o_to_h2l(K = K, P = 0.5, h2o = ldsc_h2o) 
    lambda <- median(sumstats$STAT^2) / qchisq(0.5, 1)
    
    # summary table
    prs_extremes_summary[[length(prs_extremes_summary) + 1]] <- c(overlap, n_cases_in_prs_top, n_cases_in_prs_bottom, n_prs_top, 
                                                                  n_prs_bottom, neff, ldsc_h2o, ldsc_h2l, lambda)
    
    for (snp_set_name in c("all_null_snps", "null_snps_in_prs", "null_snps_not_in_prs")) {
      
      if (snp_set_name == "all_null_snps") {snp_set <- all_null_snps}
      if (snp_set_name == "null_snps_in_prs") {snp_set <- null_snps_in_prs}
      if (snp_set_name == "null_snps_not_in_prs") {snp_set <- null_snps_not_in_prs}
      
      ## distribution of Z-stat of null snps (should be mean = 0 and var = 1)
      z_stat_var <- var(sumstats$STAT[sumstats$SNP %in% snp_set], na.rm = T)
      z_stat_mean <- mean(sumstats$STAT[sumstats$SNP %in% snp_set], na.rm = T)
      z_stat[[length(z_stat) + 1]] <- c(overlap, snp_set_name, z_stat_var, z_stat_mean)
      
      for (sig_threshold in c(0.05, 5e-8)) {
        
        ## false-positive rate
        fpr_temp <- sum(sumstats$P[sumstats$SNP %in% snp_set] < sig_threshold, na.rm = T) / length(snp_set)
        fpr[[length(fpr) + 1]] <- c(overlap, snp_set_name, sig_threshold, fpr_temp)
        
      }
    }
  }
  
  prs_extremes_summary_dfr <- as.data.frame(do.call("rbind", prs_extremes_summary))
  names(prs_extremes_summary_dfr) <- c("overlap", "n_cases_in_prs_top", "n_cases_in_prs_bottom", "n_prs_top", 
                                       "n_prs_bottom", "neff", "ldsc_h2o", "ldsc_h2l", "lambda")
  
  fpr_dfr <- as.data.frame(do.call("rbind", fpr))
  names(fpr_dfr) <- c("overlap", "snp_set_name", "sig_threshold", "fpr")
  
  z_stat_dfr <- as.data.frame(do.call("rbind", z_stat))
  names(z_stat_dfr) <- c("overlap", "snp_set_name", "z_stat_var", "z_stat_mean")  
  
                        ##### step 10: save all data #####
  
  rm(list=ls()[grep("temp",ls())])
  
  temp_prs_extremes_sumstats0 <- fread("gwas_prs_extremes0_sumstats.assoc.logistic")
  temp_prs_extremes_sumstats0$LOG_OR <- log(temp_prs_extremes_sumstats0$OR)
  temp_prs_extremes_sumstats0$LOG_OR_SE <- abs(log(temp_prs_extremes_sumstats0$OR)/qnorm(temp_prs_extremes_sumstats0$P/2))
  temp_prs_extremes_sumstats0_maf <- fread("gwas_prs_extremes0_sumstats.frq")[, c("SNP", "MAF")]
  prs_extremes_sumstats0 <- merge(x = temp_prs_extremes_sumstats0, y = temp_prs_extremes_sumstats0_maf, by = ("SNP"))
  
  
  temp_prs_extremes_sumstats50 <- fread("gwas_prs_extremes0_sumstats.assoc.logistic")
  temp_prs_extremes_sumstats50$LOG_OR <- log(temp_prs_extremes_sumstats50$OR)
  temp_prs_extremes_sumstats50$LOG_OR_SE <- abs(log(temp_prs_extremes_sumstats50$OR)/qnorm(temp_prs_extremes_sumstats50$P/2))
  temp_prs_extremes_sumstats50_maf <- fread("gwas_prs_extremes0_sumstats.frq")[, c("SNP", "MAF")]
  prs_extremes_sumstats50 <- merge(x = temp_prs_extremes_sumstats50, y = temp_prs_extremes_sumstats50_maf, by = ("SNP"))
  
  temp_prs_extremes_sumstats100 <- fread("gwas_prs_extremes0_sumstats.assoc.logistic")
  temp_prs_extremes_sumstats100$LOG_OR <- log(temp_prs_extremes_sumstats100$OR)
  temp_prs_extremes_sumstats100$LOG_OR_SE <- abs(log(temp_prs_extremes_sumstats100$OR)/qnorm(temp_prs_extremes_sumstats100$P/2))
  temp_prs_extremes_sumstats100_maf <- fread("gwas_prs_extremes0_sumstats.frq")[, c("SNP", "MAF")]
  prs_extremes_sumstats100 <- merge(x = temp_prs_extremes_sumstats100, y = temp_prs_extremes_sumstats100_maf, by = ("SNP"))
  
  rm(list=ls()[grep("temp",ls())])
  
  run_summary[[as.character(prs_threshold)]] <- list(gwas_alz_sumstats = gwas_alz_sumstats,
                      gwas_alz_summary = gwas_alz_sumstats_summary,
                      prs_extremes0_sumstats = prs_extremes_sumstats0,
                      prs_extremes50_sumstats = prs_extremes_sumstats50,
                      prs_extremes100_sumstats = prs_extremes_sumstats100,
                      prs_extremes_summary = prs_extremes_summary_dfr,
                      prs_extremes_fpr = fpr_dfr,
                      prs_extremes_z_stat = z_stat_dfr
  )
  
}

saveRDS(object = run_summary, file = "run_summary.rds")



                      ##### step 11: run 100 times on LISA #####

## paste this to start interactive node on LISA:
# module load 2020
# module load R/4.0.2-intel-2020a
# R
# 
# # load libraries
# library(data.table)
# setwd("/home/emilu/projects/paper-ad_prs_extremes/simulation_results/")

if (FALSE) { # to prevent it running as part of the simulation
  
  source("/home/emilu/projects/alzheimers/scripts/custom_functions/save_job.R")
  
  for (run in 1:100) {
    
    save_job(paste0("#!/bin/bash
                 #SBATCH --nodes=1
                 #SBATCH -t 30:00:00
                 #SBATCH --job-name=", run, "_ad_prs_extremes
                 #SBATCH --output=ad_prs_extremes_", run, ".out
                 #SBATCH --error=ad_prs_extremes_", run, ".error
                 
                 cd $TMPDIR
                 cp /home/emilu/projects/paper-ad_prs_extremes/ad_prs_extremes_simulation.R ./
                 
                 Rscript ad_prs_extremes_simulation.R
                 
                 cp run_summary.rds /home/emilu/projects/paper-ad_prs_extremes/simulation_results/run_summary_", run, ".rds"))
    
    ## submit job script and then delete it
    system("sbatch job.sh")
    system("rm job.sh")
    
  }
  
  ## rerun for failed cores:
  files <- list.files(path = "/home/emilu/projects/paper-ad_prs_extremes/simulation_results", pattern = glob2rx("run_summary_*.rds"), full.names = F)
  files <- gsub(pattern = "run_summary_", replacement = "", x = files)
  files <- as.numeric(gsub(pattern = ".rds", replacement = "", x = files))
  one_hundred <- 1:100
  missing <- one_hundred[!one_hundred %in% files]
  
  
  for (run in missing) {
    
    save_job(paste0("#!/bin/bash
                 #SBATCH --nodes=1
                 #SBATCH -t 30:00:00
                 #SBATCH --job-name=", run, "_ad_prs_extremes
                 #SBATCH --output=ad_prs_extremes_", run, ".out
                 #SBATCH --error=ad_prs_extremes_", run, ".error
                 
                 cd $TMPDIR
                 cp /home/emilu/projects/paper-ad_prs_extremes/ad_prs_extremes_simulation.R ./
                 
                 Rscript ad_prs_extremes_simulation.R
                 
                 cp run_summary.rds /home/emilu/projects/paper-ad_prs_extremes/simulation_results/run_summary_", run, ".rds"))
    
    ## submit job script and then delete it
    system("sbatch job.sh")
    system("rm job.sh")
    
  }
  
}

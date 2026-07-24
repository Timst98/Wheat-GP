
rm(list = ls())

library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scoringRules)
library(MASS)
library(doParallel)
library(doRNG)
library(purrr)
library(pracma)

`%||%` <- function(a,b) if (!is.null(a) && !is.na(a) && a != "") a else b 


id=as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))%||%501
keri=as.integer(Sys.getenv("KERNEL"))%||%5     # which kernel combination
additive_only=as.integer(Sys.getenv("ADDONLY"))%||%0
product_only=as.integer(Sys.getenv("PRODONLY"))%||%0
long_Adam=as.integer(Sys.getenv("LONG_ADAM"))%||%1
GA_kernel=as.integer(Sys.getenv("GA_KERNEL"))%||%0
NLOPTR=as.integer(Sys.getenv("NLOPTR"))%||%0 #set to 0 for ADAM
grid_search_params=as.integer(Sys.getenv("GS_PARAMS"))%||%0
train_prop=as.numeric(Sys.getenv("train_prop"))%||%0.1
leak_spec='random'
noisy_obs = TRUE
noisy_sd= TRUE
model=c("f","y")[as.numeric(Sys.getenv("model"))%||%1]

get_int <- function(name, default) {
  x <- Sys.getenv(name)
  if (!nzchar(x)) return(default)
  y <- suppressWarnings(as.integer(x))
  if (is.na(y)) default else y
}

tar <- c("yield","prot")[ get_int("TARGET", 1) ]
print(c(tar,train_prop,model))
nr_outer_split = 1


############################################################################################
# Optimization tuning parameters
batchsize= 1 # Batch the training data for 8 time speedup
NLOPT_alg='NLOPT_LD_LBFGS' 
if(long_Adam){
  lr=Lr=c(0.001)
  max_iter=1000
}else{
  lr=Lr=.1
  max_iter=100}
if(NLOPTR){Lr=1}

num_cores = length(Lr)*8
learning_rate_decay = 1


set = "1"
leakage='no'
setupp=1


nr_add=50


##############################################

source('GP_modeling/Functions/Start.R')
source('GP_modeling/Functions/GP_functions.R')



data <-
  read.table(
    "Data/DATA.txt",
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
  )
meteo  <-
  read.table(
    "Data/METEO.txt",
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
  )
geno   <-
  read.table(
    "Data/GENO.txt",
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
  )

normalize_zscore <- function(x) {
  y <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  if (sd(x, na.rm = TRUE) == 0) {
    y <- 0
  }
  return(y)
}

c2 <- ncol(meteo)
meteo[, -c2] <- lapply(meteo[, -c2], normalize_zscore)
varNames <-
  sub("_.*", "", colnames(meteo)[3:ncol(meteo)]) %>% unique()
meteo1 <- meteo %>%
  dplyr::select(matches(paste(varNames[1:7], collapse = "|")))
meteo <- cbind(meteo$Env, meteo1)
colnames(meteo)[1] <- 'Env'

GENO_SPEC = list()
for (i in 1:10) {
  GENO_SPEC[[i]] = read.table(
    paste0("Kernels/GENO_SPEC", i, ".txt"),
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
  )
}

hamming_distance_matrix <- function(maat) {
  n <- nrow(maat)
  hamming_matrix <- outer(1:n, 1:n, Vectorize(function(i, j) {
    v <- (maat[i,] == maat[j,]) | is.na(maat[i,]) | is.na(maat[j,])
    return((ncol(maat) - sum(v)) / ncol(maat))
  }))
  return(hamming_matrix)
}
geno_dist_hamming   <-
  read.table(
    "Kernels/GENO_HAMMING.txt",
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
  )
colnames(geno_dist_hamming)[colnames(geno_dist_hamming) == "SCARO_.ex_SARO."] <-
  "SCARO_(ex_SARO)"
colnames(geno_dist_hamming)[colnames(geno_dist_hamming) == "EMBLEM_.NEW."] <-
  "EMBLEM_(NEW)"
rownames(geno_dist_hamming) <- colnames(geno_dist_hamming)

#### GBLUP 
iupac_to_dosage <- function(genotypes, ref1, ref2) {
  # Map each IUPAC code to its possible alleles
  iupac_map <- list(
    A = c("A", "A"), T = c("T", "T"), C = c("C", "C"), G = c("G", "G"),
    R = c("A", "G"), Y = c("C", "T"), S = c("G", "C"), W = c("A", "T"),
    K = c("G", "T"), M = c("A", "C"))
  
  
  # Apply function to each genotype
  sapply(genotypes, function(gt) {
    alleles <- iupac_map[[gt]]
    if (is.null(alleles)) return(NA)  # Handle missing/invalid codes
    if(alleles[1] == alleles[2]){return(sum(alleles == ref1))}
    if(alleles[1] != alleles[2]){return(1)}
    
  })
}

ref <- NA
geno_t <- geno
for(i in 2:ncol(geno)){
  snp_vector <- geno[,i] %>% unname() %>% t()
  # ref
  ss <- c(sum(snp_vector=='A', na.rm = TRUE),
          sum(snp_vector=='T', na.rm = TRUE),
          sum(snp_vector=='C', na.rm = TRUE),
          sum(snp_vector=='G', na.rm = TRUE))
  ss1 <- max(ss[which(ss/nrow(geno)>0.05)])
  ss2 <- min(ss[which(ss/nrow(geno)>0.05)])
  ref1 <- c('A', 'T', 'C', 'G')[ss==ss1]
  ref2 <- c('A', 'T', 'C', 'G')[ss==ss2]
  geno_t[,i]<-iupac_to_dosage(snp_vector, ref1, ref2) %>% unname()
}

geno_numeric <- geno_t[, -1]
geno_numeric <- apply(geno_numeric, 2, function(x) as.numeric(as.character(x)))
geno_dosage <- t(geno_numeric)
rownames(geno_dosage) <- colnames(geno_t)[-1] 
colnames(geno_dosage) <- geno_t[[1]]          
D_gblup <- as.matrix(dist(t(geno_dosage)))
rownames(D_gblup) <- geno$variety_name
rownames(D_gblup) <- geno$variety_name

meteo_list <- lapply(1:nrow(meteo), function(j1) {
  meteoM <- matrix(NA, nrow = 6, ncol = 7)
  for (jj in 1:7) {
    meteoM[1:length(meteo[j1, grep(varNames[jj], colnames(meteo))]), jj] <-
      as.numeric(meteo[j1, grep(varNames[jj], colnames(meteo))])
  }
  return(meteoM)
})

data <- merge(data, meteo, by='Env')

##############################################

sr_fct_cov <- function(ind_train, ind_test,params_opt=0){
  train <- data[ind_train, , drop = FALSE]
  vv <- ind_train
  if(params_opt){
    if(!exists('opt_par')){
      source('GP_modeling/Functions/Initial_params_BO.R', local = TRUE)
      cat('initial params: ', initial_par, '\n')
      initial_params_old <<- initial_par
    }else{initial_par<<-opt_par}
  }
  if(batchsize!=1){
    cl <- makeCluster(num_cores)  
    registerDoParallel(cl)}
  source('GP_modeling/Functions/GP_training_BO2.R', local = TRUE)
  if(batchsize!=1){stopCluster(cl)}
  Kobs <- alpha * Kmat1[ind_train, ind_train, drop = FALSE] +
    (1 - alpha) * diag(nrow(train))
  Kobs_inv <- tryCatch(
    chol2inv(chol(Kobs)),
    error = function(e) ginv(as.matrix(Kobs)))
  Kobs_inv <- matrix(as.numeric(Kobs_inv), nrow = nrow(Kobs_inv), ncol = ncol(Kobs_inv))
  betahat <- sum(Kobs_inv %*% train[[tar]]) / sum(Kobs_inv)
  v <- as.numeric((t(train[[tar]] - betahat) %*% Kobs_inv %*% (train[[tar]] - betahat)) / nrow(train))
  Kobs_inv <- (1 / v) * Kobs_inv
  Kmat <- v * alpha * Kmat1
  
  tau2=v*(1-alpha)
  
  n <- nrow(train); m <- length(ind_test)
  ones <- matrix(1, n, 1)
  y <- matrix(train[[tar]], ncol = 1)
  y_centered <- y - betahat * ones
  vv <- ind_train
  
  K_starX <- Kmat[ind_test, vv, drop = FALSE]
  K_Xstar <- Kmat[vv, ind_test, drop = FALSE]
  K_starstar <- Kmat[ind_test, ind_test, drop = FALSE]
  alphaa <- Kobs_inv %*% y_centered
  pred_mean <- as.numeric(betahat + K_starX %*% alphaa)
  core <- K_starstar - K_starX %*% Kobs_inv %*% K_Xstar
  u <- matrix(1, m, 1) - K_starX %*% (Kobs_inv %*% ones)
  denom <- as.numeric(t(ones) %*% Kobs_inv %*% ones)
  cov_extra <- (u %*% t(u)) / max(denom, 1e-12)
  covGP <- core + cov_extra
  covGP <- (covGP + t(covGP)) / 2
  pred_sd <- sqrt(pmax(diag(covGP), 0))


  if(model=='y'){pred_sd=sqrt(pmax(diag(covGP), 0)+tau2)}  
  return(list(pred = pred_mean, pred_sd = pred_sd, covGP = covGP))  }

select_train_ind <- function(setupp, leakage, leak_spec = "random") {
  if (setupp == '1') {
    unique_envs <- unique(data$Env)
    split_envs <- sample(unique_envs)
    train_envs <-
      split_envs[1:floor(length(split_envs) * train_prop)]
    test_envs <-
      split_envs[(floor(length(split_envs) * train_prop) + 1):length(split_envs)]
  }
  
  if (setupp == '2') {
    unique_vars <- unique(data$variety)
    split_vars <- sample(unique_vars)
    train_vars <-
      split_vars[1:floor(length(split_vars) * train_prop)]
    test_vars <-
      split_vars[(floor(length(split_vars) * train_prop) + 1):length(split_vars)]
  }
  
  if (setupp == '1' & leakage == 'no') {
    ind_train <- which(data$Env %in% train_envs)
  }
  
  if (setupp == '1' & leakage == 'yes' & leak_spec == 'standard') {
    ind_train <-
      which(data$Env %in% train_envs |
              data$variety_name %in% leak_name)
  }
  
  if (setupp == '1' & leakage == 'yes' & leak_spec == 'random') {
    ind_train <- which(data$Env %in% train_envs)
    # one random leak per env
    ind_leak <- NaN
    for (ttt in 1:length(test_envs)) {
      vvs <- data %>% filter(Env == test_envs[ttt]) %>% pull(variety_name)
      v_rand <- sample(vvs)[1]
      ind_leak[ttt] <-
        which(data$Env == test_envs[ttt] & data$variety_name == v_rand)
    }
    ind_train <- c(ind_train, ind_leak) %>% sort()
  }
  
  if (setupp == '2' & leakage == 'no') {
    ind_train <- which(data$variety_name %in% train_vars)
  }
  
  if (setupp == '2' & leakage == 'yes' & leak_spec == 'standard') {
    ind_train <-
      which(data$variety_name %in% train_vars |
              data$Env %in% leak_name)
  }
  
  if (setupp == '2' & leakage == 'yes' & leak_spec == 'random') {
    ind_train <- which(data$variety_name %in% train_vars)
    # one random leak per env
    ind_leak <- NaN
    for (ttt in 1:length(test_vars)) {
      vvs <- data %>% filter(variety_name == test_vars[ttt]) %>% pull(Env)
      v_rand <- sample(vvs)[1]
      ind_leak[ttt] <-
        which(data$variety_name == test_vars[ttt] & data$Env == v_rand)
    }
    ind_train <- c(ind_train, ind_leak) %>% sort()
  }
  
  return(ind_train)
  
}



##############################################
## for evaluatiobnn
# optimizte hyperparameters and fit

filename=paste0(
       'Results/Results BO/Full_fit_',
       lr,
       '_',
       max_iter,
       '_',
       tar,
       keri,
       leakage,
       setupp,
       '_',
       id,
       '.rda'
     )


      batchsize=1
      xx_train <- sr_fct_cov(1:nrow(data), 1:nrow(data),params_opt=1)
      data$tar_clean <- xx_train$pred
      res=list(data=data,opt_par=opt_par,initial_par=initial_par)
      save(list = "res",
        file = filename
     )
cat("full data initial params: ",initial_par,"\n")
  cat("full data optimal params: ",opt_par,"\n")


cat('Saved results under: ', filename, '.\n')

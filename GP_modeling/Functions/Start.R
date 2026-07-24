
suppressPackageStartupMessages({
library(nloptr)
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(dtw)
library(dbscan)
library(FNN)
library(caTools)
library(Metrics)
library(rsample)
library(Rfast)
library(lubridate)
library(patchwork)
library(tibble)
library(kernlab)
library(rrBLUP)
library(Rcpp)
library(Metrics)
library(dplyr)
library(tidyverse)
library(scoringRules)
library(foreach)
library(doParallel)
library(kernlab)
library(Matrix)
library(MASS)
library(progress)
library(mvtnorm)
library(BGLR)
library(lme4)
})

UB=c(1,20,1,1,1)
LB=rep(1e-6,5)
`%||%` <- function(a,b) if (!is.null(a) && !is.na(a) && a != "") a else b
mm2 = inds = ntr = #train = 
  c()
Export=c(
  "Adam",
  "additive_only",
  "batchsize",
  "create_kernelmat_derivative",
  "crps_norm",
  "D_gblup",
  "data",
  "GA_kernel",
  "GENO_SPEC",
  "geno",
  "geno_dist_hamming",
  "ginv",
  "GP_test",
  "grid",
  "id",
  "inds",
  "kernel_geno",
  "kernel_meteo",
  "keri",
  "learning_rate_decay",
  "likelihood_gradient",
  "log_likelihood",
  "logGAK",
  "logGAK_derivative",
  "LOG0",
  "LOGP",
  'LogS',
  "logs_norm",
  "long_Adam",
  "lr",
  "lr_tar",
  "Lr",
  "max_iter",
  "meteo",
  "meteo_list",
  "mm2",
  "mse",
  "NLOPT_alg",
  "NLOPTR",
  "noisy_obs",
  "notcluster",
  "nr_outer_split",
  "Packages",
  "product_only",
  "tar",
  "train",
  "LB",
  "UB",
  'grid_search_params'
)


if (keri == 1) {
  kernel_geno='dos'
  kernel_meteo='rbf'
}
if (keri == 2) {
  kernel_geno='ham'
  kernel_meteo='rbf'
}

if (keri == 3) {
  kernel_geno='spe'
  kernel_meteo='rbf'
}

if (keri == 4) {
  kernel_geno='dos'
  kernel_meteo='exp'
}
if (keri == 5) {
  kernel_geno='ham'
  kernel_meteo='exp'
}

if (keri == 6) {
  kernel_geno='spe'
  kernel_meteo='exp'
}

if (keri == 7) {
  kernel_geno='dos'
  GA_kernel=1
  kernel_meteo='ali'
}
if (keri == 8) {
  kernel_geno='ham'
  GA_kernel=1
  kernel_meteo='ali'
}

if (keri == 9) {
  kernel_geno='spe'
  GA_kernel=1
  kernel_meteo='ali'
}

Packages=c(
  'BGLR',
  'nloptr',
  'parallel',
  'foreach',
  'doParallel',
  'dplyr',
  'scoringRules',
  'Metrics',
  'progress',
  'rrBLUP',
  'pracma'
)

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




# normalize the different weather variables
normalize_zscore <- function(x) {
  y <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  if (sd(x, na.rm = TRUE) == 0) {
    y <- 0
  }
  return(y)
}
normalize_min_max <- function(x) {
  y <-
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  if (max(x, na.rm = TRUE) == 0 & min(x, na.rm = TRUE)) {
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




#################################################################################
# HAMMING for GENO

# count the number of equal or NA entries
hamming_distance_matrix <- function(maat) {
  n <- nrow(maat)
  hamming_matrix <- outer(1:n, 1:n, Vectorize(function(i, j) {
    v <- (maat[i,] == maat[j,]) | is.na(maat[i,]) | is.na(maat[j,])
    return((ncol(maat) - sum(v)) / ncol(maat))
  }))
  return(hamming_matrix)
}
# geno_dist_hamming <- hamming_distance_matrix(geno[,-1])
# rownames(geno_dist_hamming) <- geno$variety_name
# colnames(geno_dist_hamming) <- geno$variety_name
# write.table(geno_dist_hamming, file = "GENO_HAMMING.txt", sep = "\t", row.names = FALSE, quote = FALSE)
geno_dist_hamming   <-
  read.table(
    "Kernels/GENO_HAMMING.txt",
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE
  )







#################################################################################
# Alignment for meteo  (from Marco Cuturi's webpage)
# Useful constants
LOG0 <- -10000  # log(0)

# LOGP function: stable computation of log(exp(x) + exp(y))
LOGP <- function(x, y) {
  if (x > y) {
    return(x + log1p(exp(y - x)))
  } else {
    return(y + log1p(exp(x - y)))
  }
}

# Global Alignment Kernel function
logGAK <- function(seq1, seq2, sigma, triangular = 0) {
  nX <- nrow(seq1)  # Number of rows in seq1
  nY <- nrow(seq2)  # Number of rows in seq2
  dimvect <- ncol(seq1)  # Dimension of the time series vectors
  
  # Length of a column for dynamic programming
  cl <- nY + 1
  
  # Initialize logM to store two successive columns (dynamic programming table)
  logM <- matrix(LOG0, nrow = cl, ncol = 2)
  logM[1, 1] <- 0  # Initialize the lower-left cell with log(1) = 0
  
  # Maximum of abs(i - j) when 1 <= i <= nX and 1 <= j <= nY
  trimax <- max(nX - 1, nY - 1)
  
  # Triangular coefficients initialization
  logTriangularCoefficients <- rep(0, trimax + 1)
  if (triangular > 0) {
    for (i in seq_len(min(trimax + 1, triangular)) - 1) {
      logTriangularCoefficients[i + 1] <- log(1 - i / triangular)
    }
  }
  
  # Sigma factor for Gaussian kernel
  Sig <- -1 / (2 * sigma ^ 2)
  
  # Dynamic programming to compute the log of the Global Alignment Kernel
  cur <- 2
  old <- 1
  
  for (i in 1:nX) {
    logM[, cur] <- LOG0  # Reset the current column
    for (j in 1:nY) {
      if (logTriangularCoefficients[abs(i - j) + 1] > LOG0) {
        # Compute the Gaussian kernel value
        diff_sq <- sum((seq1[i,] - seq2[j,]) ^ 2)
        gram <-
          logTriangularCoefficients[abs(i - j) + 1] + diff_sq * Sig
        gram <- gram - log(2 - exp(gram))
        
        # Update logM for dynamic programming
        frompos1 <- logM[j + 1, old]
        frompos2 <- logM[j, cur]
        frompos3 <- logM[j, old]
        aux <- LOGP(frompos1, frompos2)
        logM[j + 1, cur] <- LOGP(aux, frompos3) + gram
      }
    }
    # Swap cur and old
    cur <- 3 - cur
    old <- 3 - old
  }
  
  # Return the final result
  return(logM[nY + 1, old])
}

logGAK_derivative <- function(seq1, seq2, sigma, triangular = 0) {
  nX <- nrow(seq1)  # Number of rows in seq1
  nY <- nrow(seq2)  # Number of rows in seq2
  dimvect <- ncol(seq1)  # Dimension of the time series vectors
  
  # Length of a column for dynamic programming
  cl <- nY + 1
  
  # Initialize logM to store two successive columns (dynamic programming table)
  logM <- matrix(LOG0, nrow = cl, ncol = 2)
  logM[1, 1] <- 0  # Initialize the lower-left cell with log(1) = 0
  
  # Maximum of abs(i - j) when 1 <= i <= nX and 1 <= j <= nY
  trimax <- max(nX - 1, nY - 1)
  
  # Triangular coefficients initialization
  logTriangularCoefficients <- rep(0, trimax + 1)
  if (triangular > 0) {
    for (i in seq_len(min(trimax + 1, triangular)) - 1) {
      logTriangularCoefficients[i + 1] <- log(1 - i / triangular)
    }
  }
  
  # Sigma factor for Gaussian kernel
  Sig <- -1 / (2 * sigma ^ 2)
  del_Sig = sigma ^ (-3)
  # Dynamic programming to compute the log of the Global Alignment Kernel
  cur <- 2
  old <- 1
  
  for (i in 1:nX) {
    logM[, cur] <- LOG0  # Reset the current column
    for (j in 1:nY) {
      if (logTriangularCoefficients[abs(i - j) + 1] > LOG0) {
        # Compute the Gaussian kernel value
        diff_sq <- sum((seq1[i,] - seq2[j,]) ^ 2)
        gram <-
          logTriangularCoefficients[abs(i - j) + 1] + diff_sq * Sig
        gram <- gram - log(2 - exp(gram))
        del_gram = diff_sq * del_Sig * (1 + 1 / (2 - exp(gram)) * exp(gram))
        
        # Update logM for dynamic programming
        frompos1 <- logM[j + 1, old]
        frompos2 <- logM[j, cur]
        frompos3 <- logM[j, old]
        aux <- LOGP(frompos1, frompos2)
        logM[j + 1, cur] <- LOGP(aux, frompos3) + del_gram
      }
    }
    # Swap cur and old
    cur <- 3 - cur
    old <- 3 - old
  }
  
  # Return the final result
  return(logM[nY + 1, old])
}

varNames <-
  sub("_.*", "", colnames(meteo)[3:ncol(meteo)]) %>% unique()


################################################################################
# FUNCTIONS
make_load_env <- function(parent_env,
                          tar, leakage, setupp,
                          data, meteo,
                          id, keri,
                          additive_only, product_only,
                          grid_search_params,
                          NLOPTR, NLOPT_alg, batchsize, lr, max_iter,
                          MODEL) {

  e <- new.env(parent = parent_env)

  # task-specific
  e$tar     <- tar
  e$leakage <- leakage
  e$setupp  <- setupp

  # big objects
  e$data  <- data
  e$meteo <- meteo

  # run config
  e$id <- id
  e$keri <- keri
  e$additive_only <- additive_only
  e$product_only  <- product_only
  e$grid_search_params <- grid_search_params
  e$NLOPTR <- NLOPTR
  e$NLOPT_alg <- NLOPT_alg
  e$batchsize <- batchsize
  e$lr <- lr
  e$max_iter <- max_iter
  e$MODEL <- MODEL

  e
}
select_train_ind <- function(setupp, leakage, leak_spec) {
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


########## ADAM optimizer ######################################################
notcluster=is.na(id)
Adam = function(f,
                grad_f,
                params_init,
                lr = 0.01,
                beta1 = 0.9,
                beta2 = 0.999,
                epsilon = 1e-8,
                max_iter = 100,
                tol = 1e-10,
                lb = rep(0, length(params_init)),
                ub = rep(Inf, length(params_init))) {
  # Initialize progress bar if not in cluster

    pb = progress_bar$new(
      format = "  Progress [:bar] :percent in :elapsed",
      total = max_iter,
      clear = FALSE,
      width = 60
    )
  
  
  params = params_init
  m = numeric(length(params))
  v = numeric(length(params))
  iter = 0
  
  while (iter < max_iter) {
    iter = iter + 1
    
    # Compute gradient
    grad = grad_f(params)
    
    # Update moments
    m = beta1 * m + (1 - beta1) * grad
    v = beta2 * v + (1 - beta2) * (grad ^ 2)
    
    # Compute bias-corrected moments
    m_hat = m / (1 - beta1 ^ iter)
    v_hat = v / (1 - beta2 ^ iter)
    
    # Update parameters
    params_new = params - lr * m_hat / (sqrt(v_hat) + epsilon)
    
    # Handle boundary constraints
    params_new[is.na(params_new)] = params[is.na(params_new)]
    params_new[params_new <= lb] = params[params_new <= lb]
    params_new[params_new > ub] = params[params_new > ub]
    if ((iter / max_iter * 100) %in% seq(0, 100, by = 5)) {
      print(c("iter=", iter, ", grad=", grad))
      if (learning_rate_decay) {
        lr = .9 * lr
      }
    }
    #if(iter%%50==0){print(grad)}
    #print(grad)
    # Check convergence
    if (max(abs(params_new - params)) < tol) {
      break
    }
    
    params = params_new
    
    
    
      pb$tick()
    
  }
  
  return(params)
}



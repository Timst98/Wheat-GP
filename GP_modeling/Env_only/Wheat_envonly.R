
print('Version 27.06 10:10')
id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
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
  #library(dtwclust)
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
############################################################################################
keri <- as.integer(Sys.getenv("KERNEL"))     # which kernel combination
setupp <- ifelse(as.integer(Sys.getenv("SETUP"))==1,'1','2') # new variety or new environment
leakage <- ifelse(as.integer(Sys.getenv("LEAKAGE"))==0,'no','yes')
leak_spec <- 'random'
noisy_obs = TRUE
train_prop = 0.8
nr_outer_split = 1
GA_kernel=as.integer(Sys.getenv("GA_KERNEL"))
notcluster=is.na(id)
noisy_sd=TRUE
############################################################################################
# Optimization tuning parameters
batchsize=.5 # Batch the training data for 8 time speedup

NLOPTR=as.integer(Sys.getenv("NLOPTR")) #set to 0 for ADAM
NLOPT_alg='NLOPT_LD_LBFGS' 
additive_only=1
product_only=0

#"NLOPT_LD_TNEWTON_PRECOND_RESTART"
long_Adam=as.integer(Sys.getenv("LONG_ADAM"))
if(long_Adam){
  lr=.01
  max_iter=1000
}else{
  lr=.1
  max_iter=100}

num_cores = 2 / batchsize*ifelse(keri%in%3:4&&GA_kernel==1,25,1)
learning_rate_decay = 1

#for conditional likelihood 
validation_gradient=0 
val_ratio=.2
val_sets=10

if (validation_gradient) {
  num_cores = val_sets*num_cores
}

####################### KERNELS #############################################################
if (keri == 1) {
  kernel_geno <- 'dos'
  kernel_meteo <- 'rbf'
}
if (keri == 2) {
  kernel_geno <- 'ham'
  kernel_meteo <- 'rbf'
}

if (keri == 3) {
  kernel_geno <- 'spe'
  kernel_meteo <- 'rbf'
}

if (keri == 4) {
  kernel_geno <- 'dos'
  kernel_meteo <- 'exp'
}
if (keri == 5) {
  kernel_geno <- 'ham'
  kernel_meteo <- 'exp'
}

if (keri == 6) {
  kernel_geno <- 'spe'
  kernel_meteo <- 'exp'
}

if (keri == 7) {
  kernel_geno <- 'dos'
  GA_kernel=1
  kernel_meteo <- 'ali'
}
if (keri == 8) {
  kernel_geno <- 'ham'
  GA_kernel=1
  kernel_meteo <- 'ali'
}

if (keri == 9) {
  kernel_geno <- 'spe'
  GA_kernel=1
  kernel_meteo <- 'ali'
}

grid= as.matrix(read.csv(paste0("GP_modeling/grids/grid_keri_ENV",keri,".csv")))
  
  LB=apply(grid,2,min)
  UB=apply(grid,2,max)

  


print(LB)
print(UB)
Export = c(
  'kernel_geno',
  'kernel_meteo',
  'D_gblup',
  'grid',
  'keri',
  'train',
  'ntr',
  'batchsize',
  'val_ratio',
  'create_kernelmat_derivative',
  'create_kernelmat_derivative',
  'Adam',
  'max_iter',
  'lr',
  'inds',
  'log_likelihood',
  'val_sets',
  'NLOPT_alg',
  'GENO_SPEC',
  'mm2',
  'likelihood_gradient',
  'ginv',
  'tar',
  'noisy_obs',
  'meteo',
  'data',
  'geno',
  'GA_kernel',
  'additive_only',
  'product_only',
  'geno_dist_hamming',
  'noisy_sd'
  
)
Packages=c(
  'nloptr',
  'parallel',
  'foreach',
  'doParallel',
  'dplyr',
  'scoringRules',
  'Metrics',
  'progress','pracma'
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

#data$yield <- (data$yield - mean(data$yield)) / sd(data$yield)
#data$prot <- (data$prot - mean(data$prot)) / sd(data$prot)

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
# create the kernel matrices


meteo_list <- lapply(1:nrow(meteo), function(j1) {
  meteoM <- matrix(NA, nrow = 6, ncol = 7)
  for (jj in 1:7) {
    meteoM[1:length(meteo[j1, grep(varNames[jj], colnames(meteo))]), jj] <-
      as.numeric(meteo[j1, grep(varNames[jj], colnames(meteo))])
  }
  return(meteoM)
})

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

create_kernelmat_derivative <-
    function(kernel_meteo,
             kernel_geno,
             theta1,
             mat) {
      if (kernel_meteo == 'exp') {
        meteo_dist_eucl <- as.matrix(dist(meteo[, -1], method = "euclidean"))
        rownames(meteo_dist_eucl) <- meteo$Env
        colnames(meteo_dist_eucl) <- meteo$Env
        int <- theta1 * max(meteo_dist_eucl)  
        Dmat <- meteo_dist_eucl[mat$Env, mat$Env]
        Kmat_meteo <- exp(-Dmat / (int))
        Kmat_meteo_deriv <-
          exp(-Dmat / (int)) * Dmat / (int ^ 2) * max(meteo_dist_eucl)
      }
      
      # rbf kernel 
      if (kernel_meteo == 'rbf') {
        meteo_dist_eucl <- as.matrix(dist(meteo[, -1], method = "euclidean"))
        rownames(meteo_dist_eucl) <- meteo$Env
        colnames(meteo_dist_eucl) <- meteo$Env
        int <-
          theta1 * max(meteo_dist_eucl)   # c(0.2,0.3,0.4,0.5,0.6)[i1]
        Dmat <- meteo_dist_eucl[mat$Env, mat$Env]
        Kmat_meteo <- exp(-Dmat^2 / (2*int^2))
        Kmat_meteo_deriv <-
          Kmat_meteo * Dmat^2 / (int ^ 3) * max(meteo_dist_eucl)
      }
      
      # alignments kernel
      if (kernel_meteo == 'ali') {
        Kernel_GA = Kernel_GA_deriv =
          matrix(NA, nrow = nrow(meteo), ncol = nrow(meteo))
        rownames(Kernel_GA) = meteo$Env
        colnames(Kernel_GA) = meteo$Env
        rownames(Kernel_GA_deriv) = meteo$Env
        colnames(Kernel_GA_deriv) = meteo$Env
        
        results =
          foreach(
            j1 = 1:nrow(meteo),
            .combine = 'c',
            .export = c(
              'meteo_list',
              'logGAK',
              'logGAK_derivative',
              'LOG0',
              'LOGP'
            )
          ) %dopar% {
            meteoM1 = meteo_list[[j1]]
            local_results =
              list() 
            
            for (j2 in 1:j1) {
              meteoM2 = meteo_list[[j2]]
              
              t1 = logGAK(meteoM1,
                          meteoM2,
                          sigma = theta1,
                          triangular = 0)
              t1_deriv =logGAK_derivative(meteoM1,
                                          meteoM2,
                                          sigma = theta1,
                                          triangular = 0)
              t2 =logGAK(meteoM2,
                         meteoM2,
                         sigma = theta1,
                         triangular = 0)
              t2_deriv =logGAK_derivative(meteoM2,
                                          meteoM2,
                                          sigma = theta1,
                                          triangular = 0)
              t3 =logGAK(meteoM1,
                         meteoM1,
                         sigma = theta1,
                         triangular = 0)
              t3_deriv = logGAK_derivative(meteoM1,
                                           meteoM1,
                                           sigma = theta1,
                                           triangular = 0)
              
              exp_val = exp(t1 - 0.5 * (t2 + t3))
              deriv_val =exp_val * (t1_deriv - 0.5 * (t2_deriv + t3_deriv))
              
              
              local_results[[length(local_results) + 1]] =
                list(
                  j1 = j1,
                  j2 = j2,
                  kernel_val = exp_val,
                  kernel_deriv = deriv_val
                )
            }
            
            return(local_results)  
          } 
        for (res in results) {
          j1 = res$j1
          j2 = res$j2
          Kernel_GA[j1, j2] = res$kernel_val
          Kernel_GA[j2, j1] = res$kernel_val
          Kernel_GA_deriv[j1, j2] = res$kernel_deriv
          Kernel_GA_deriv[j2, j1] = res$kernel_deriv
        }
        Kmat_meteo_deriv = Kernel_GA_deriv[mat$Env, mat$Env]
        Kmat_meteo = Kernel_GA[mat$Env, mat$Env]
      }
      
      
      Kmat=as.matrix(Kmat_meteo)
      return(
        list(
          Kmat = Kmat,
          del_theta1 = as.matrix( Kmat_meteo_deriv)
           
        )
      )
    }
  
################################################################################
# FUNCTIONS
GP_test <-
  function(data,
           train,
           test,
           betahat,
           Kobs_inv,
           Kmat,
           tt,
           vv) {
    ind_tt <-
      which(data$Env == test$Env[tt] &
              data$variety == test$variety[tt])[1]
    # where is the current test data in data
    Kcross <- Kmat[ind_tt, vv] %>% as.numeric()
    
    meanGP <-
      betahat + t(Kcross) %*% (Kobs_inv %*% (train[[tar]] - betahat * matrix(
        1, nrow = nrow(train), ncol = 1
      )))
    varGP <- Kmat[ind_tt, ind_tt] -
      t(Kcross) %*% Kobs_inv %*% Kcross +
      (1 - t(matrix(
        1, nrow = nrow(train), ncol = 1
      )) %*% Kobs_inv %*% Kcross) ^ 2 / (t(matrix(
        1, nrow = nrow(train), ncol = 1
      )) %*% Kobs_inv %*% matrix(1, nrow = nrow(train), ncol = 1))
    
    pred <- meanGP
    sd <- ifelse(noisy_sd,sqrt(varGP)+sqrt(v*(1-alpha)),sqrt(varGP))
    LogS=-2*dnorm(test[[tar]][tt],meanGP,sd,log = TRUE)
    return(c(pred, sd,LogS))
  }

LogS <-
  function(data,
           train,
           test,
           betahat,
           Kobs_inv,
           Kmat,
           tt,
           vv) {
    ind_tt <-
      sapply(tt, function(t) which(data$Env == test$Env[t] &
                                     data$variety == test$variety[t])[1])
    # where is the current test data in data
    Kcross <- Kmat[ind_tt, vv] %>% as.matrix()
    
    meanGP <-
      betahat + (Kcross) %*% (Kobs_inv %*% (train[[tar]] - betahat * matrix(
        1, nrow = nrow(train), ncol = 1
      )))
    
    Kobs_invsum=(t(matrix(
      1, nrow = nrow(train), ncol = 1
    )) %*% Kobs_inv %*% matrix(1, nrow = nrow(train), ncol = 1))%>%as.numeric
    outer_vec <- 1 - Kcross %*% (Kobs_inv %*% rep(1, nrow(train)))
    varGP <- Kmat[ind_tt, ind_tt] -
      (Kcross) %*% Kobs_inv %*% t(Kcross) +(outer_vec %*% t(outer_vec))%>%as.matrix / Kobs_invsum
    
    
    
    nte=length(test[[tar]])
    y=test[[tar]]
    prec=chol2inv(varGP)
    LogS=t(y-meanGP)%*%prec%*%(y-meanGP)+
      determinant(varGP)$modulus[1]+log(2*pi)*nte
    
    return(LogS)
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

################Continuous likelihood function##################################

log_likelihood = function(theta1,
                          alpha,
                          train) {
  
  
  Kmat1 <- create_kernelmat_derivative(kernel_meteo,
                                       kernel_geno,
                                         theta1,
                                         train)$Kmat

  
  if (noisy_obs == FALSE) {
    Kobs <- Kmat1
  }
  if (noisy_obs == TRUE) {
    Kobs <- alpha * Kmat1 + (1 - alpha) * diag(nrow(train))
  }
  Kobs_inv <- tryCatch(
    chol2inv(chol(Kobs)),
    # Attempt Cholesky-based inversion
    error = function(e) {
      # Handle errors
      message("Cholesky decomposition failed. Falling back to generalized inverse (ginv).")
      ginv(as.matrix(Kobs))          # Use ginv() from MASS as a fallback
    }
  )
  
  Kobs_inv <-
    matrix(as.numeric(Kobs_inv),
           nrow = nrow(Kobs_inv),
           ncol = ncol(Kobs_inv))
  betahat <- sum(Kobs_inv %*% train[[tar]]) / sum(Kobs_inv)
  
  v = as.numeric(1 / nrow(train) * t(train[[tar]] - betahat) %*% Kobs_inv %*%
                   (train[[tar]] - betahat))
  #Kmat=v*alpha*Kmat1
  Kobs_inv = 1 / v * Kobs_inv
  Kobs = v * Kobs
  zn = train[[tar]]
  
  return(
    nrow(train) * log(2 * pi) + determinant(as.matrix(Kobs))$modulus[1] +
      t(zn - betahat) %*% Kobs_inv %*% (zn - betahat)
  )
  
}



likelihood_gradient = function(theta1,
                               alpha,
                               train) {
  Trace = function(x) {
    sum(diag(x))
  }
  
  
    
    kernel_matrices  <- create_kernelmat_derivative(kernel_meteo,
                                                      kernel_geno,
                                                      theta1,
                                                      train)
      
    
    
    Kmat1 <- kernel_matrices$Kmat
    if (noisy_obs == FALSE) {
      Kobs <- Kmat1
    }
    if (noisy_obs == TRUE) {
      Kobs <- alpha * Kmat1 + (1 - alpha) * diag(nrow(train))
    }
    Kobs_inv <- tryCatch(
      chol2inv(chol(Kobs)),
      # Attempt Cholesky-based inversion
      error = function(e) {
        # Handle errors
        message("Cholesky decomposition failed. Falling back to generalized inverse (ginv).")
        ginv(as.matrix(Kobs))          # Use ginv() from MASS as a fallback
      }
    )
    
    Kobs_inv <-
      matrix(as.numeric(Kobs_inv),
             nrow = nrow(Kobs_inv),
             ncol = ncol(Kobs_inv))
    betahat <- sum(Kobs_inv %*% train[[tar]]) / sum(Kobs_inv)
    
    v = as.numeric(1 / nrow(train) * t(train[[tar]] - betahat) %*% Kobs_inv %*%
                     (train[[tar]] - betahat))
    
    Kmat = v * alpha * Kmat1
    Kobs_inv = 1 / v * Kobs_inv
    Kobs = v * Kobs
    zn = train[[tar]]
    
   
      ############################################################################
      
      Kobs_inv_deriv = -(v * Kobs_inv) %*% as.matrix(alpha * kernel_matrices$del_theta1) %*%
        (v * Kobs_inv)
      
      del_theta1_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
        betahat / (sum(v * Kobs_inv) ^ 2) %*% sum(Kobs_inv_deriv)
      
      del_v_theta1 = -2 * del_theta1_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                                   (zn - betahat)) -
        as.numeric(
          1 / nrow(train) * t(zn - betahat) %*% (v * Kobs_inv) %*% Kobs_inv_deriv %*%
            (v * Kobs_inv) %*% (zn - betahat)
        )
      del_v_theta1 = del_v_theta1 %>% as.numeric()
      
      
      Kobs_inv_deriv_theta1 = -Kobs_inv %*% as.matrix(v * alpha * kernel_matrices$del_theta1) %*%
        Kobs_inv +
        
        del_v_theta1 * v * Kobs_inv
      
      ############################################################################
      
      Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(alpha * kernel_matrices$del_theta1) %*%
        as.matrix(v * Kobs_inv)
      
      
      del_theta1_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
        sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
      
      #v=as.numeric(1/nrow(train)*t(zn-betahat)%*%Kobs_inv%*%(zn-betahat))
      
      del_v_theta1 = -2 * del_theta1_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                                   (zn - betahat)) +
        as.numeric(1 / nrow(train) * t(zn - betahat) %*% Kobs_inv_deriv %*%
                     (zn - betahat)) %>% as.numeric()
      
      
      Kobs_inv_deriv_theta1 = as.matrix(1 / v * Kobs_inv_deriv - 1 / (v ^
                                                                        2) * del_v_theta1 * (Kobs_inv * v))
      
      
      #########################################################################################################
      
      
      #########################################################################################################
     
      #########################################################################################################
      
      
      ###########################################################################################################
      
      Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(Kmat1 - diag(nrow(train))) %*%
        as.matrix(v * Kobs_inv)
      
      
      del_alpha_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
        sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
      
      #v=as.numeric(1/nrow(train)*t(zn-betahat)%*%Kobs_inv%*%(zn-betahat))
      
      del_v_alpha = -2 * del_alpha_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                                 (zn - betahat)) +
        as.numeric(1 / nrow(train) * t(zn - betahat) %*% Kobs_inv_deriv %*%
                     (zn - betahat))
      del_v_alpha = del_v_alpha %>% as.numeric()
      
      Kobs_inv_deriv_alpha = as.matrix(1 / v * Kobs_inv_deriv - 1 / (v ^ 2) *
                                         del_v_alpha * (Kobs_inv * v))
      
      #########################################################################################################
      
      grad = c(
        Trace(
          Kobs_inv %*% as.matrix(v * alpha * kernel_matrices$del_theta1 + del_v_theta1 *
                                   Kobs / v)
        ) -
          2 * del_theta1_betahat %*% rep(1, nrow(train)) %*% (Kobs_inv %*% (zn -
                                                                              betahat)) -
          t(zn - betahat) %*% Kobs_inv_deriv_theta1 %*% (zn - betahat),
        
        
        Trace(Kobs_inv %*% as.matrix(v * (
          Kmat1 - diag(nrow(train))
        ) + del_v_alpha * Kobs / v)) -
          2 * del_alpha_betahat %*% rep(1, nrow(train)) %*% Kobs_inv %*% (zn -
                                                                            betahat) +
          t(zn - betahat) %*% Kobs_inv_deriv_alpha %*% (zn - betahat)
        
        
        
      )
      
      
  
  
  return(grad)
  
  
}


########## ADAM optimizer ######################################################

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
  if (notcluster) {
    pb = progress_bar$new(
      format = "  Progress [:bar] :percent in :elapsed",
      total = max_iter,
      clear = FALSE,
      width = 60
    )
  }
  
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
        lr = .8 * lr
      }
    }
    #if(iter%%50==0){print(grad)}
    #print(grad)
    # Check convergence
    if (max(abs(params_new - params)) < tol) {
      break
    }
    
    params = params_new
    
    
    if (notcluster) {
      pb$tick()
    }
  }
  
  return(params)
}



# ######## UNCOMMENT FOR GRADIENT CHECK ##########################################
# 
if(id==1){
    ind_train=sample(nrow(data),10);ind_test=sample((1:nrow(data))[-ind_train],10)
    train <- data[ind_train, ]
    test <- data[-ind_train, ]
    vv <- ind_train
    ntr = nrow(train)
    cl <- makeCluster(num_cores)   # Adjust the number of cores as needed
    registerDoParallel(cl)
    tryCatch({
         nloptr(runif(2),function(x) {log_likelihood(x[1], x[2],,train)},
                 function(x) {likelihood_gradient(x[1], x[2],train)},
                 opts=list(algorithm="NLOPT_LD_LBFGS",check_derivatives=TRUE))
    }, error = function(e) { message("") 
      
    })
    stopCluster(cl)
    
}
################################################################################


mm2 = inds = ntr = train = c()

cl <- makeCluster(num_cores)  # Adjust the number of cores as needed
registerDoParallel(cl)
Res = foreach(
  tar = c('yield','prot'),
  .combine = 'list',
  .export = c(Export, 'tar', 'noisy_obs'),
  .packages = Packages
) %dopar% {
  hyperpp <- matrix(NA, 3, 10)
  
  #  tau <- mean(data[[tar]]*0.01)
  
  #res <- matrix(NA, nr_outer_split, ifelse(  kernel_geno == 'spe',24,25))
  res <- matrix(NA, nr_outer_split, 26)
  
  for (i in 1:nr_outer_split) {
    if(GA_kernel){
      
        load(paste0('Results/Results_env_only/Results grid search/Initial_par_',
                    tar, keri,'GA', leakage, setupp, '_', id,
                    '.rda'))
      
      
    }else{
      
        load(paste0('Results/Results_env_only/Results grid search/Initial_par_',
                    tar, keri, leakage, setupp, '_', id,
                    '.rda'))
      
      
    }
    
    Initial_params = Initial_res
    ind_train = Initial_params[[1]]
    train <- data[ind_train, ]
    test <- data[-ind_train, ]
    vv <- ind_train
    ntr = nrow(train)
    
    initial_par = Initial_params[[2]]
    
    if(sum(LB>initial_par)>0){
      LB[which(LB>initial_par)]=initial_par[which(LB>initial_par)]-1e-4
    }
    
    if(sum(UB<initial_par)>0){
      UB[which(UB<initial_par)]=initial_par[which(UB<initial_par)]+1e-4
    }
    
      if (NLOPTR) {
      
          opt_par = nloptr(initial_par, function(x) {
            set.seed(sum(x))
            inds = matrix(sample(1:nrow(train)), ncol = 1 / batchsize)
            apply(
              foreach(
                ii = 1:(1 / batchsize),
                .combine = 'cbind',
                .export = c(
                  Export,
                  'meteo_list',
                  'logGAK',
                  'LOG0',
                  'LOGP',
                  'logGAK_derivative'
                ),
                .packages =Packages
              ) %dopar% {
                log_likelihood(x[1], x[2],
                               train[inds[, ii], ])
              },
              1,
              mean
            )
          }, function(x) {
            set.seed(sum(x))
            inds = matrix(sample(1:nrow(train)), ncol = 1 / batchsize)
            apply(
              foreach(
                ii = 1:(1 / batchsize),
                .combine = 'cbind',
                .export = c(
                  Export,
                  'meteo_list',
                  'logGAK',
                  'LOG0',
                  'LOGP',
                  'logGAK_derivative'
                ),
                .packages = Packages
              ) %dopar% {
                likelihood_gradient(x[1], x[2],
                                    train[inds[, ii], ])
              },
              1,
              mean
            )
            
          }, opts = list(algorithm = NLOPT_alg, print_level = 1),
          lb = LB,
          ub =UB)$solution
      
        
      } else{
       tic()
          opt_par = Adam(function(x) {
            log_likelihood(x[1], x[2],  train)
          }, function(x) {
            inds = matrix(sample(1:nrow(train)), ncol = 1 / batchsize)
            apply(foreach(
              ii = 1:(1 / batchsize),
              .combine = 'cbind',
              .export = Export,
              .packages=Packages
            ) %dopar% {
              likelihood_gradient(x[1], x[2],
                                  train[inds[, ii], ])
            },
            1,
            mean)
          }, initial_par, lr = lr, max_iter = max_iter, beta1 = .6, beta2 = .95,
          lb = LB,
          ub =UB)
        time_optimization=toc()
        
      }
      
tic()
        Kmat1 <-create_kernelmat_derivative(kernel_meteo, 
                                            kernel_geno,
                                            opt_par[1],
                                         
                                            data)$Kmat
        alpha=opt_par[2]
        
     
      
      if (noisy_obs == FALSE) {
        Kobs <- Kmat1[vv, vv]
      }
      if (noisy_obs == TRUE) {
        Kobs <- alpha * Kmat1[vv, vv] + (1 - alpha) * diag(nrow(train))
      }
      
      
      Kobs_inv <- tryCatch(
        chol2inv(chol(Kobs)),
        # Attempt Cholesky-based inversion
        error = function(e) {
          # Handle errors
          message("Cholesky decomposition failed. Falling back to generalized inverse (ginv).")
          ginv(as.matrix(Kobs))          # Use ginv() from MASS as a fallback
        }
      )
      Kobs_inv <-
        matrix(as.numeric(Kobs_inv),
               nrow = nrow(Kobs_inv),
               ncol = ncol(Kobs_inv))
      betahat <- sum(Kobs_inv %*% train[[tar]]) / sum(Kobs_inv)
      
      v = as.numeric(1 / nrow(train) * t(train[[tar]] - betahat) %*% Kobs_inv %*%
                       (train[[tar]] - betahat))
      
      Kobs_inv = 1 / v * Kobs_inv
      Kmat = v * alpha * Kmat1
      
      results <-
        sapply(1:nrow(test), function(tt)
          GP_test(data, train, test, betahat, Kobs_inv, Kmat, tt, vv))
      
      test$pred <- results[1, ]
      test$sd <- (results[2, ])
      test$LogS=results[3, ]
      # test <- test[complete.cases(test),]
      
      res[i, 1] = mse(test$pred, test[[tar]])
      res[i, 2] = mean(crps_norm(test[[tar]], test$pred, test$sd)) # include tau2!
      time_inference=toc()
               time_total=time_inference+time_optimization
      # global average
      mean_train <- train[[tar]] %>% mean()
      sd_train <- train[[tar]] %>% sd()
      res[i, 3] = mse(test[[tar]], mean_train)
      res[i, 4] = mean(abs(test[[tar]] - mean_train))
      
      # variety average
      if (leakage == 'yes' || setupp == '1') {
        pred_var <- NaN
        for (j in 1:nrow(test)) {
          pred_var[j] <- train %>%
            filter(variety_name == test$variety_name[j]) %>%
            pull(tar) %>% mean() %>% as.numeric()
        }
        loca <- 1 - is.nan(pred_var)
        res[i, 5] = mse(test[[tar]][loca == 1], pred_var[loca == 1])
        res[i, 6] = mean(abs(test[[tar]][loca == 1] - pred_var[loca == 1]))
      }
      
      
      # environmental average
      if (leakage == 'yes' || setupp == '2') {
        pred_loc <- NaN
        for (j in 1:nrow(test)) {
          pred_loc[j] <- train %>% ungroup() %>%
            filter(Env == test$Env[j]) %>%
            pull(tar) %>% mean() %>% as.numeric()
        }
        loca <- 1 - is.nan(pred_loc)
        res[i, 7] = mse(test[[tar]][loca == 1], pred_loc[loca == 1])
        res[i, 8] = mean(abs(test[[tar]][loca == 1] - pred_loc[loca == 1]))
      }
      
        res[i, 9] = log_likelihood(opt_par[1],
                                   opt_par[2],
                                   
                                   train)
        #res[i,10]=mean(test$LogS)
        res[i,10]=mean(logs_norm(test[[tar]], test$pred, test$sd))
        res[i,11]=LogS(data, train, test, betahat, Kobs_inv,
                       Kmat, 1:nrow(test), vv)
        res[i, 13:14] = opt_par
        opt_par = initial_par
        res[i, 12] = log_likelihood(opt_par[1],
                                    opt_par[2],
                                    
                                    train)
        res[i, 19:20] = opt_par
        
      res[i, 24:26] = c(time_optimization,time_inference,time_total)
        
      }
      
  res
  
}


if(GA_kernel==0){
  for (i in 1:2) {
    tar = c('yield', 'prot')[i]
    res = as.numeric(Res[[i]])
    if(additive_only){
      if (NLOPTR) {
        if (validation_gradient) {
          
          save(
            list = "res",
            file = paste0(
              'Results/Results_env_only/Results LBFGS/Results_hyper_VAL_',
              NLOPT_alg,
              '_',
              batchsize,
              '_',
              tar,
              keri,
              leakage,
              setupp,
              '_',
              ifelse(noisy_sd,'noisy_sd_',''),id,
              '.rda'
            )
          )
        } else{
          
          save(
            list = "res",
            file = paste0(
              'Results/Results_env_only/Results LBFGS/Results_hyper_',
              NLOPT_alg,
              '_',
              batchsize,
              '_',
              tar,
              keri,
              leakage,
              setupp,
              '_',
              ifelse(noisy_sd,'noisy_sd_',''),id,
              '.rda'
            )
          )
        }
      } else{
        if (validation_gradient) {
          save(
            list = "res",
            file = paste0(
              'Results/Results_env_only/Results Adam/Results_hyper_VAL_ADAM_',
              lr,
              '_',
              max_iter,
              '_',
              tar,
              keri,
              leakage,
              setupp,
              '_',
              ifelse(noisy_sd,'noisy_sd_',''),id,
              '.rda'
            )
          )
        } else{
          save(
            list = "res",
            file = paste0(
              'Results/Results_env_only/Results Adam/Results_hyper_ADAM_',
              lr,
              '_',
              max_iter,
              '_',
              batchsize,
              '_',
              tar,
              keri,
              leakage,
              setupp,
              '_',
              ifelse(noisy_sd,'noisy_sd_',''),id,
              '.rda'
            )
          )
        }
      }
    }else{
      if(product_only){
        
        if (NLOPTR) {
          if (validation_gradient) {
            save(
              list = "res",
              file = paste0(
                'Results/Results_product_only/Results LBFGS/Results_hyper_VAL_',
                NLOPT_alg,
                '_',
                batchsize,
                '_',
                tar,
                keri,
                leakage,
                setupp,
                '_',
                ifelse(noisy_sd,'noisy_sd_',''),id,
                '.rda'
              )
            )
          } else{
            
            save(
              list = "res",
              file = paste0(
                'Results/Results_product_only/Results LBFGS/Results_hyper_',
                NLOPT_alg,
                '_',
                batchsize,
                '_',
                tar,
                keri,
                leakage,
                setupp,
                '_',
                ifelse(noisy_sd,'noisy_sd_',''),id,
                '.rda'
              )
            )
          }
        } else{
          if (validation_gradient) {
            save(
              list = "res",
              file = paste0(
                'Results/Results_product_only/Results Adam/Results_hyper_VAL_ADAM_',
                lr,
                '_',
                max_iter,
                '_',
                tar,
                keri,
                leakage,
                setupp,
                '_',
                ifelse(noisy_sd,'noisy_sd_',''),id,
                '.rda'
              )
            )
          } else{
            save(
              list = "res", 
              file = paste0(
                'Results/Results_product_only/Results Adam/Results_hyper_ADAM_',
                lr,
                '_',
                max_iter,
                '_',
                batchsize,
                '_',
                tar,
                keri,
                leakage,
                setupp,
                '_',
                ifelse(noisy_sd,'noisy_sd_',''),id,
                '.rda'
              )
            )
          }
        }
      }else{
        
        if (NLOPTR) {
          if (validation_gradient) {
            save(
              list = "res",
              file = paste0(
                'Results/Results LBFGS/Results_hyper_VAL_',
                NLOPT_alg,
                '_',
                batchsize,
                '_',
                tar,
                keri,
                leakage,
                setupp,
                '_',
                ifelse(noisy_sd,'noisy_sd_',''),id,
                '.rda'
              )
            )
          } else{
            
            save(
              list = "res",
              file = paste0(
                'Results/Results LBFGS/Results_hyper_',
                NLOPT_alg,
                '_',
                batchsize,
                '_',
                tar,
                keri,
                leakage,
                setupp,
                '_',
                ifelse(noisy_sd,'noisy_sd_',''),id,
                '.rda'
              )
            )
          }
        } else{
          if (validation_gradient) {
            save(
              list = "res",
              file = paste0(
                'Results/Results Adam/Results_hyper_VAL_ADAM_',
                lr,
                '_',
                max_iter,
                '_',
                tar,
                keri,
                leakage,
                setupp,
                '_',
                ifelse(noisy_sd,'noisy_sd_',''),id,
                '.rda'
              )
            )
          } else{
            save(
              list = "res", 
              file = paste0(
                'Results/Results Adam/Results_hyper_ADAM_',
                lr,
                '_',
                max_iter,
                '_',
                batchsize,
                '_',
                tar,
                keri,
                leakage,
                setupp,
                '_',
                ifelse(noisy_sd,'noisy_sd_',''),id,
                '.rda'
              )
            )
          }
        }
      }
    }
    
  }
}else{
  for (i in 1:2) {
    tar = c('yield', 'prot')[i]
    res = as.numeric(Res[[i]])
    if(additive_only){
      if (NLOPTR) {
        if (validation_gradient) {
          
          save(
            list = "res",
            file = paste0(
              'Results/Results_env_only/Results LBFGS/Results_hyper_VAL_',
              NLOPT_alg,
              '_',
              batchsize,
              '_',
              tar,
              keri,
              'GA',
              leakage,
              setupp,
              '_',
              ifelse(noisy_sd,'noisy_sd_',''),id,
              '.rda'
            )
          )
        } else{
          
          save(
            list = "res",
            file = paste0(
              'Results/Results_env_only/Results LBFGS/Results_hyper_',
              NLOPT_alg,
              '_',
              batchsize,
              '_',
              tar,
              keri,
              'GA',
              leakage,
              setupp,
              '_',
              ifelse(noisy_sd,'noisy_sd_',''),id,
              '.rda'
            )
          )
        }
      } else{
        if (validation_gradient) {
          save(
            list = "res",
            file = paste0(
              'Results/Results_env_only/Results Adam/Results_hyper_VAL_ADAM_',
              lr,
              '_',
              max_iter,
              '_',
              tar,
              keri,
              'GA',
              leakage,
              setupp,
              '_',
              ifelse(noisy_sd,'noisy_sd_',''),id,
              '.rda'
            )
          )
        } else{
          save(
            list = "res",
            file = paste0(
              'Results/Results_env_only/Results Adam/Results_hyper_ADAM_',
              lr,
              '_',
              max_iter,
              '_',
              batchsize,
              '_',
              tar,
              keri,
              'GA',
              leakage,
              setupp,
              '_',
              ifelse(noisy_sd,'noisy_sd_',''),id,
              '.rda'
            )
          )
        }
      }
    }else{
      if(product_only){
        
        if (NLOPTR) {
          if (validation_gradient) {
            save(
              list = "res",
              file = paste0(
                'Results/Results_product_only/Results LBFGS/Results_hyper_VAL_',
                NLOPT_alg,
                '_',
                batchsize,
                '_',
                tar,
                keri,
                'GA',
                leakage,
                setupp,
                '_',
                ifelse(noisy_sd,'noisy_sd_',''),id,
                '.rda'
              )
            )
          } else{
            
            save(
              list = "res",
              file = paste0(
                'Results/Results_product_only/Results LBFGS/Results_hyper_',
                NLOPT_alg,
                '_',
                batchsize,
                '_',
                tar,
                keri,
                'GA',
                leakage,
                setupp,
                '_',
                ifelse(noisy_sd,'noisy_sd_',''),id,
                '.rda'
              )
            )
          }
        } else{
          if (validation_gradient) {
            save(
              list = "res",
              file = paste0(
                'Results/Results_product_only/Results Adam/Results_hyper_VAL_ADAM_',
                lr,
                '_',
                max_iter,
                '_',
                tar,
                keri,
                'GA',
                leakage,
                setupp,
                '_',
                ifelse(noisy_sd,'noisy_sd_',''),id,
                '.rda'
              )
            )
          } else{
            save(
              list = "res", 
              file = paste0(
                'Results/Results_product_only/Results Adam/Results_hyper_ADAM_',
                lr,
                '_',
                max_iter,
                '_',
                batchsize,
                '_',
                tar,
                keri,
                'GA',
                leakage,
                setupp,
                '_',
                ifelse(noisy_sd,'noisy_sd_',''),id,
                '.rda'
              )
            )
          }
        }
      }else{
        
        if (NLOPTR) {
          if (validation_gradient) {
            save(
              list = "res",
              file = paste0(
                'Results/Results LBFGS/Results_hyper_VAL_',
                NLOPT_alg,
                '_',
                batchsize,
                '_',
                tar,
                keri,
                'GA',
                leakage,
                setupp,
                '_',
                ifelse(noisy_sd,'noisy_sd_',''),id,
                '.rda'
              )
            )
          } else{
            
            save(
              list = "res",
              file = paste0(
                'Results/Results LBFGS/Results_hyper_',
                NLOPT_alg,
                '_',
                batchsize,
                '_',
                tar,
                keri,
                'GA',
                leakage,
                setupp,
                '_',
                ifelse(noisy_sd,'noisy_sd_',''),id,
                '.rda'
              )
            )
          }
        } else{
          if (validation_gradient) {
            save(
              list = "res",
              file = paste0(
                'Results/Results Adam/Results_hyper_VAL_ADAM_',
                lr,
                '_',
                max_iter,
                '_',
                tar,
                keri,
                'GA',
                leakage,
                setupp,
                '_',
                ifelse(noisy_sd,'noisy_sd_',''),id,
                '.rda'
              )
            )
          } else{
            save(
              list = "res", 
              file = paste0(
                'Results/Results Adam/Results_hyper_ADAM_',
                lr,
                '_',
                max_iter,
                '_',
                batchsize,
                '_',
                tar,
                keri,
                'GA',
                leakage,
                setupp,
                '_',
                ifelse(noisy_sd,'noisy_sd_',''),id,
                '.rda'
              )
            )
          }
        }
      }
    }
    
  }
}


stopCluster(cl)

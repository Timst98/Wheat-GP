
id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print("version 26.6. 13.15")
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
  library(dtwclust)
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
  library(lhs)
  library(R.utils)
})
############################################################################################

keri <- as.integer(Sys.getenv("KERNEL"))     # which kernel combination
setupp <- ifelse(as.integer(Sys.getenv("SETUP"))==1,'1','2') # new variety or new environment
leakage <- ifelse(as.integer(Sys.getenv("LEAKAGE"))==0,'no','yes')
leak_spec <- 'random'
additive_only=1
product_only=0
GA_kernel=as.integer(Sys.getenv("GA_KERNEL"))

noisy_obs = TRUE
train_prop = 0.8
nr_outer_split = 1 # Do all of them seperately

num_cores=120
ngrid=1000
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

if(additive_only){
  grid= maximinLHS(ngrid,2)
  na_rows=sapply(1:nrow(grid),function(i){ifelse(sum(is.na(as.numeric(grid[i,]))>0),1,0)})
  grid=grid[which(na_rows==0),]
    
  write.csv(grid, file = paste0("GP_modeling/grids/grid_keri_ENV",keri,".csv"), row.names = FALSE)
  
  
}



print(nrow(grid))

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
  Sig <- -1 / (2 * sigma^2)
  
  # Dynamic programming to compute the log of the Global Alignment Kernel
  cur <- 2
  old <- 1
  
  for (i in 1:nX) {
    logM[, cur] <- LOG0  # Reset the current column
    for (j in 1:nY) {
      if (logTriangularCoefficients[abs(i - j) + 1] > LOG0) {
        # Compute the Gaussian kernel value
        diff_sq <- sum((seq1[i, ] - seq2[j, ])^2)
        gram <- logTriangularCoefficients[abs(i - j) + 1] + diff_sq * Sig
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
  Sig <- -1 / (2 * sigma^2)
  del_Sig=sigma^(-3)
  # Dynamic programming to compute the log of the Global Alignment Kernel
  cur <- 2
  old <- 1
  
  for (i in 1:nX) {
    logM[, cur] <- LOG0  # Reset the current column
    for (j in 1:nY) {
      if (logTriangularCoefficients[abs(i - j) + 1] > LOG0) {
        # Compute the Gaussian kernel value
        diff_sq <- sum((seq1[i, ] - seq2[j, ])^2)
        gram <- logTriangularCoefficients[abs(i - j) + 1] + diff_sq * Sig
        gram <- gram - log(2 - exp(gram))
        del_gram=diff_sq * del_Sig*(1+1/(2-exp(gram))*exp(gram))
        
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

varNames <- sub("_.*", "", colnames(meteo)[3:ncol(meteo)]) %>% unique()



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
###################################################################################
# create the kernel matrices (continuous in parameters, plus derivatives)


if(additive_only){
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
          del_theta1 = as.matrix(Kmat_meteo_deriv)
           
        )
      )
    }
  
}




meteo_list <- lapply(1:nrow(meteo), function(j1) {
  meteoM <- matrix(NA, nrow = 6, ncol = 7)
  for (jj in 1:7) {
    meteoM[1:length(meteo[j1, grep(varNames[jj], colnames(meteo))]), jj] <- 
      as.numeric(meteo[j1, grep(varNames[jj], colnames(meteo))])
  }
  return(meteoM)
})



###################################################################################
# FUNCTIONS
GP_test <- function(data, train, test, betahat, Kobs_inv, Kmat, tt, vv){
  
  ind_tt <- which(data$Env == test$Env[tt] & data$variety == test$variety[tt])[1]   
  # where is the current test data in data
  Kcross <- Kmat[ind_tt, vv] %>% as.numeric()
  
  meanGP <- betahat + t(Kcross)%*%(Kobs_inv%*%(
    train[[tar]] - betahat*matrix(1, nrow = nrow(train), ncol = 1)))
  varGP <- Kmat[ind_tt,ind_tt] - 
    t(Kcross)%*%Kobs_inv%*%Kcross + 
    (1-t(matrix(1, nrow = nrow(train), ncol = 1))%*%
       Kobs_inv%*%Kcross)^2/(t(matrix(1, nrow = nrow(train), ncol = 1))%*%
                               Kobs_inv%*%matrix(1, nrow = nrow(train), ncol = 1))
  
  pred <- meanGP
  sd <- sqrt(varGP)
  
  return(c(pred, sd))    
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

select_train_ind <- function(setupp, leakage, leak_spec){
  if(setupp == '1'){
    unique_envs <- unique(data$Env)
    split_envs <- sample(unique_envs)
    train_envs <- split_envs[1:floor(length(split_envs) * train_prop)]
    test_envs <- split_envs[(floor(length(split_envs) * train_prop) + 1):length(split_envs)]}
  
  if(setupp == '2'){
    unique_vars <- unique(data$variety)
    split_vars <- sample(unique_vars)
    train_vars <- split_vars[1:floor(length(split_vars) * train_prop)]
    test_vars <- split_vars[(floor(length(split_vars) * train_prop) + 1):length(split_vars)]}
  
  if(setupp == '1' & leakage == 'no'){
    ind_train <- which(data$Env %in% train_envs)  }
  
  if(setupp == '1' & leakage == 'yes' & leak_spec == 'standard'){
    ind_train <- which(data$Env %in% train_envs | data$variety_name %in% leak_name)   }
  
  if(setupp == '1' & leakage == 'yes' & leak_spec == 'random'){
    ind_train <- which(data$Env %in% train_envs)
    # one random leak per env
    ind_leak <- NaN
    for(ttt in 1:length(test_envs)){
      vvs <- data %>% filter(Env == test_envs[ttt]) %>% pull(variety_name)
      v_rand <- sample(vvs)[1]
      ind_leak[ttt] <- which(data$Env == test_envs[ttt] & data$variety_name == v_rand)}
    ind_train <- c(ind_train, ind_leak) %>% sort()  }
  
  if(setupp == '2' & leakage == 'no'){
    ind_train <- which(data$variety_name %in% train_vars)   }
  
  if(setupp == '2' & leakage == 'yes' & leak_spec == 'standard'){
    ind_train <- which(data$variety_name %in% train_vars | data$Env %in% leak_name)   }
  
  if(setupp == '2' & leakage == 'yes' & leak_spec == 'random'){
    ind_train <- which(data$variety_name %in% train_vars)
    # one random leak per env
    ind_leak <- NaN
    for(ttt in 1:length(test_vars)){
      vvs <- data %>% filter(variety_name == test_vars[ttt]) %>% pull(Env)
      v_rand <- sample(vvs)[1]
      ind_leak[ttt] <- which(data$variety_name == test_vars[ttt] & data$Env == v_rand)}
    ind_train <- c(ind_train, ind_leak) %>% sort()  }
  
  return(ind_train)
  
}


cl <- makeCluster(num_cores)  # Adjust the number of cores as needed
registerDoParallel(cl)

for(ss in 1:2){
  print(c("ss=",ss))
  tar <- c('yield', 'prot')[ss]
  
  hyperpp <- matrix(NA, 3, 10)
  
  tau <- mean(data[[tar]]*0.01)
  # 1 per cent obs noise, can also be fitted as hyperparameter
  
  res <- matrix(NA,nr_outer_split,8)
  #[1] 0.112 0.001 1.000 0.001
  
  
  for(i in 1:nr_outer_split){
    
    # split train and test
    set.seed(id)
    ind_train <- select_train_ind(setupp, leakage, leak_spec)
    train <- data[ind_train,]
    test <- data[-ind_train,]
    vv <- ind_train
    
    
    logLH <- foreach(j = 1:nrow(grid), .combine = 'c',
                 .export = c('grid', 'create_kernelmat_derivative', 'train',
                             'kernel_meteo', 'kernel_geno', 'product_only',
                             'additive_only', 'noisy_obs', 'tar'),
                 .packages = c("mvtnorm", "MASS", "R.utils", "doParallel",
                               "foreach")) %dopar% {
                   
                   tryCatch({
                     withTimeout({
                       
             
                      alpha <- grid[j, 2]
                      
                       
             
                      Kmat1=create_kernelmat_derivative(kernel_meteo,
                                                           kernel_geno, 
                                                           grid[j,1],
                                                           train)$Kmat
                      
                       
                       if(noisy_obs==TRUE){Kobs <-(alpha*Kmat1 + (1-alpha)*diag(nrow(train))) }
                       
                       
                       
                       Kobs_inv <- tryCatch(
                         chol2inv(chol(Kobs)),            # Attempt Cholesky-based inversion
                         error = function(e) {            # Handle errors
                           message("Cholesky decomposition failed. 
                                   Falling back to generalized inverse (ginv).")
                           ginv(as.matrix(Kobs))          # Use ginv() from MASS as a fallback
                         }
                       )
                       
                       Kobs_inv <-  matrix(as.numeric(Kobs_inv), nrow = nrow(Kobs_inv), ncol = ncol(Kobs_inv)) 
                       
                       betahat <- sum(Kobs_inv%*%train[[tar]])/sum(Kobs_inv)
                       v=as.numeric(1/nrow(train)*t(train[[tar]]-betahat)%*%Kobs_inv%*%(train[[tar]]-betahat))
                       
                       Kobs=v*Kobs
                       
                       mvtnorm::dmvnorm(train[[tar]], mean = rep(betahat, nrow(train)),  
                                        sigma = Kobs, log = TRUE)
                     }, timeout = 300, onTimeout = "error") 
                   }, error = function(e) {
                     message("logLH computation failed or timed out: ", e$message)
                     NaN
                   })
                   
                }
      
      print("finished grid search")
      print(c('Number of NaN vals: ',sum(is.na(logLH))))
      
      res1=res2=res3=res4=NaN
      while(sum(c(res1,res2,res3,res4)%>%is.nan)>0||
            sum(c(res1,res2,res3,res4)%>%is.na)>0){
        mm <- which(logLH == max(logLH[!is.na(logLH)]))[1]
        theta1 <- grid[mm,1]
        alpha <- grid[mm,2]
          initial_par=c(theta1,alpha)%>%unname
          Kmat1=create_kernelmat_derivative(kernel_meteo, kernel_geno, 
                                            theta1,
                                          
                                            data)$Kmat
          
        
        
     
      
      if(noisy_obs==FALSE){Kobs <- Kmat[vv,vv] }
      if(noisy_obs==TRUE){Kobs <- alpha*Kmat1[vv,vv] + (1-alpha)*diag(nrow(train)) }
      
      
      Kobs_inv <- tryCatch(
        chol2inv(chol(Kobs)),            # Attempt Cholesky-based inversion
        error = function(e) {            # Handle errors
          message("Cholesky decomposition failed. Falling back to generalized inverse (ginv).")
          ginv(as.matrix(Kobs))          # Use ginv() from MASS as a fallback
        }
      )
      
      Kobs_inv <-  matrix(as.numeric(Kobs_inv), nrow = nrow(Kobs_inv), ncol = ncol(Kobs_inv)) 
      betahat <- sum(Kobs_inv%*%train[[tar]])/sum(Kobs_inv)
      v=as.numeric(1/nrow(train)*t(train[[tar]]-betahat)%*%Kobs_inv%*%(train[[tar]]-betahat))
      
      Kobs_inv=1/v*Kobs_inv
      Kmat=v*alpha*Kmat1
      results <- sapply(1:nrow(test), function(tt)
        GP_test(data, train, test, betahat, Kobs_inv, Kmat, tt, vv))
      
      test$pred <- results[1,]
      test$sd <- (results[2,])
      # test <- test[complete.cases(test),]
      
      res1 = mse(test$pred, test[[tar]])
      res2 = mean(crps_norm(test[[tar]], test$pred, test$sd)) # include tau2!
      res3= mean(logs_norm(test[[tar]], test$pred, test$sd))
      res4=LogS(data, train, test, betahat, Kobs_inv,
                Kmat, 1:nrow(test), vv)
      
      logLH[mm]=NA
      }
    }
    
    Initial_res=list(ind_train,initial_par,-2*max(logLH),res1,res2,res3,res4,logLH)
    

                        
      
  if(GA_kernel){
       save(list = "Initial_res", 
           file = paste0('Results/Results_env_only/Results grid search/Initial_par_',
                         tar, keri,'GA', leakage,setupp,'_',id,'.rda'))
      
  }else{
    
      save(list = "Initial_res", 
           file = paste0('Results/Results_env_only/Results grid search/Initial_par_',
                         tar, keri, leakage,setupp,'_',id,'.rda'))
      
    
    
    
  }
  
  
  print(c("INITIAL_PAR: ", initial_par))  
  print(c("RESULTS: ",res1,res2,res3,res4))
}
stopCluster(cl)

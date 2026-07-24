
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

#num_cores=120
#ngrid=100
num_cores=1
ngrid=100
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


grid= maximinLHS(ngrid,5)
na_rows=sapply(1:nrow(grid),function(i){ifelse(sum(is.na(as.numeric(grid[i,]))>0),1,0)})
grid=grid[which(na_rows==0),]
  
if(kernel_geno == 'spe'){
  grid[,2]=as.integer(1+10* grid[,2])
}

if(!additive_only&&!product_only){
  write.csv(grid, file = paste0("GP_modeling/grids/grid_keri",keri,".csv"), row.names = FALSE)
  
}

if(additive_only){
  grid=grid[,-5]
  write.csv(grid, file = paste0("GP_modeling/grids/grid_keri_additive",keri,".csv"), row.names = FALSE)

}

if(product_only){
  grid=grid[,-c(4,5)]
  write.csv(grid, file = paste0("GP_modeling/grids/grid_keri_product",keri,".csv"), row.names = FALSE)
}

UB=grid%>%apply(2,max)%>%round
LB=grid%>%apply(2,min)%>%round


print(nrow(grid))

###################################################################################

create_kernelmat_derivative <-
  function(kernel_meteo,
           kernel_geno,
           params=c(theta1,theta2,k,w1,w2,alpha),
           mat) {
   
      theta1=params[1]
      if(kernel_geno=='spe'){
        k=params[2]
      }else{
        theta2=params[2]
      }
      if(additive_only){
        w1=params[3]
      }
      if(!additive_only&&!product_only){
        w1=params[3]
        w2=params[4]
      }
    
    
    
    
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
          .export = c('train',
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
    
    
    if (kernel_geno == 'ham') {
      int <-theta2 * max(geno_dist_hamming)   # seq(0.1, 1.5, length.out = 10)[i2]*max(geno_dist_hamming)
      Dmat <- geno_dist_hamming[mat$variety_name, mat$variety_name]
      Kmat_geno <- exp(-Dmat / (int))
      Kmat_geno_deriv <-
        Kmat_geno * Dmat / (int ^ 2) * max(geno_dist_hamming)
    }
    
    if(kernel_geno=='dos'){
      int <-theta2 * max(D_gblup)   # seq(0.1, 1.5, length.out = 10)[i2]*max(geno_dist_hamming)
      Dmat <- D_gblup[mat$variety_name, mat$variety_name]
      Kmat_geno <- exp(-Dmat^2 / (2*int^2))       # exp(-Dmat / (int)) sollte Gaussian kernel sein!
      Kmat_geno_deriv <-
        Kmat_geno * Dmat^2 / (int ^ 3) * max(D_gblup)     #  Kmat_geno * Dmat / (int ^ 2) * max(D_gblup)
    }
    
    # spectrum kernel
    if (kernel_geno == 'spe') {
      s <- k
      geno_spec1 = GENO_SPEC[[s]]
     
      geno_spec <-
        geno_spec1[mat$variety_name, mat$variety_name] / max(geno_spec1)
      rownames(geno_spec) <- colnames(geno_spec)
      Kmat_geno_deriv = Kmat_geno <- geno_spec
    }
    
    if(kernel_meteo=='ali'){
      if(kernel_geno=='spe'){
        if(product_only){
          return(
            list(
              Kmat = as.matrix(Kmat_meteo*Kmat_geno),
                Kernel_GA=Kernel_GA,
                Kernel_GA_deriv=Kernel_GA_deriv
            )
          )
        }else{
          if(additive_only){
            return(
              list(
                Kmat = as.matrix(w1*Kmat_meteo+(1-w1)*Kmat_geno),
                del_w1=as.matrix(Kmat_meteo-Kmat_geno),
                Kernel_GA=Kernel_GA,
                Kernel_GA_deriv=Kernel_GA_deriv
              ))
          }else{
            return(
              list(
                Kmat = as.matrix(w1*(1-w2)*Kmat_meteo+w1*w2*Kmat_geno+(1-w1)*Kmat_meteo*Kmat_geno),
                del_w1=as.matrix((1-w2)*Kmat_meteo+w2*Kmat_geno-Kmat_meteo*Kmat_geno),
                del_w2=as.matrix(-w1*Kmat_meteo+w1*Kmat_geno),
                Kernel_GA=Kernel_GA,
                Kernel_GA_deriv=Kernel_GA_deriv
              )
            )
          }
        }
      }else{
        if(product_only){
          return(
            list(
              Kmat = as.matrix(Kmat_meteo*Kmat_geno),
              del_theta2 = as.matrix(Kmat_geno_deriv* Kmat_meteo),
              Kernel_GA=Kernel_GA,
              Kernel_GA_deriv=Kernel_GA_deriv
            )
          )
        }else{
          if(additive_only){
            return(
              list(
                Kmat = as.matrix(w1*Kmat_meteo+(1-w1)*Kmat_geno),
                del_theta2 = as.matrix((1-w1)*Kmat_geno_deriv),
                del_w1=as.matrix(Kmat_meteo-Kmat_geno),
                Kernel_GA=Kernel_GA,
                Kernel_GA_deriv=Kernel_GA_deriv
              )
            )
          }else{
            return(
              list(
                Kmat = as.matrix(w1*(1-w2)*Kmat_meteo+w1*w2*Kmat_geno+(1-w1)*Kmat_meteo*Kmat_geno),
                del_theta2 = as.matrix(w1*(w2)*Kmat_geno_deriv+(1-w1)*Kmat_meteo*Kmat_geno_deriv),
                del_w1=as.matrix((1-w2)*Kmat_meteo+w2*Kmat_geno-Kmat_meteo*Kmat_geno),
                del_w2=as.matrix(-w1*Kmat_meteo+w1*Kmat_geno),
                Kernel_GA=Kernel_GA,
                Kernel_GA_deriv=Kernel_GA_deriv
              )
            )
          }
        }
      }
    }else{
      if(kernel_geno=='spe'){
        if(product_only){
          return(
            list(
              Kmat = as.matrix(Kmat_meteo*Kmat_geno),
              del_theta1 = as.matrix(Kmat_meteo_deriv*Kmat_geno)
            )
          )
        }else{
          if(additive_only){
            return(
              list(
                Kmat = as.matrix(w1*Kmat_meteo+(1-w1)*Kmat_geno),
                del_theta1 = as.matrix(w1*Kmat_meteo_deriv),
                del_w1=as.matrix(Kmat_meteo-Kmat_geno)
              )
            )
          }else{
            return(
              list(
                Kmat = as.matrix(w1*(1-w2)*Kmat_meteo+w1*w2*Kmat_geno+(1-w1)*Kmat_meteo*Kmat_geno),
                del_theta1 = as.matrix(w1*(1-w2)*Kmat_meteo_deriv+(1-w1)*Kmat_meteo_deriv*Kmat_geno),
                del_w1=as.matrix((1-w2)*Kmat_meteo+w2*Kmat_geno-Kmat_meteo*Kmat_geno),
                del_w2=as.matrix(-w1*Kmat_meteo+w1*Kmat_geno)
              )
            )
          }
        }
      }else{
        if(product_only){
          return(
            list(
              Kmat = as.matrix(Kmat_meteo*Kmat_geno),
              del_theta1 = as.matrix(Kmat_meteo_deriv*Kmat_geno),
              del_theta2 = as.matrix(Kmat_geno_deriv* Kmat_meteo)
            )
          )
        }else{
          if(additive_only){
            return(
              list(
                Kmat = as.matrix(w1*Kmat_meteo+(1-w1)*Kmat_geno),
                del_theta1 = as.matrix(w1*Kmat_meteo_deriv),
                del_theta2 = as.matrix((1-w1)*Kmat_geno_deriv),
                del_w1=as.matrix(Kmat_meteo-Kmat_geno)
              )
            )
          }else{
            return(
              list(
                Kmat = as.matrix(w1*(1-w2)*Kmat_meteo+w1*w2*Kmat_geno+(1-w1)*Kmat_meteo*Kmat_geno),
                del_theta1 = as.matrix(w1*(1-w2)*Kmat_meteo_deriv+(1-w1)*Kmat_meteo_deriv*Kmat_geno),
                del_theta2 = as.matrix(w1*(w2)*Kmat_geno_deriv+(1-w1)*Kmat_meteo*Kmat_geno_deriv),
                del_w1=as.matrix((1-w2)*Kmat_meteo+w2*Kmat_geno-Kmat_meteo*Kmat_geno),
                del_w2=as.matrix(-w1*Kmat_meteo+w1*Kmat_geno)
              )
            )
          }
        }
      }
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


logLH <- foreach(j = 1:nrow(grid), .combine="c", .export=c("grid","ind_train","data",
                           "create_kernelmat_derivative",
                           "kernel_meteo","kernel_geno","noisy_obs","tar",
                           "meteo","geno_dist_hamming","D_gblup","GENO_SPEC","meteo_list","varNames",
                           "additive_only","product_only"),
                 .packages=c("mvtnorm","MASS","R.utils")) %dopar% {
  train_local <-  data[ind_train, , drop = FALSE]

                 

                   
                   tryCatch({
                     withTimeout({
                       
                       kernel_params=grid[j,-ncol(grid)]
                       alpha=grid[j,ncol(grid)]
                       
                       
                       
                       Kmat1=create_kernelmat_derivative(kernel_meteo,
                                                         kernel_geno, 
                                                         kernel_params,
                                                         train_local)$Kmat
                       
                       if(noisy_obs==TRUE){Kobs <-(alpha*Kmat1 + (1-alpha)*diag(nrow(train_local))) }
                       
                       
                       
                       Kobs_inv <- tryCatch(
                         chol2inv(chol(Kobs)),            # Attempt Cholesky-based inversion
                         error = function(e) {            # Handle errors
                           message("Cholesky decomposition failed. Falling back to generalized inverse (ginv).")
                           ginv(as.matrix(Kobs))          # Use ginv() from MASS as a fallback
                         }
                       )
                       
                       Kobs_inv <-  matrix(as.numeric(Kobs_inv), nrow = nrow(Kobs_inv), ncol = ncol(Kobs_inv)) 
                       
                       betahat <- sum(Kobs_inv%*%train_local[[tar]])/sum(Kobs_inv)
                       v=as.numeric(1/nrow(train_local)*t(train_local[[tar]]-betahat)%*%Kobs_inv%*%(train_local[[tar]]-betahat))
                       
                       Kobs=v*Kobs
                       
                       mvtnorm::dmvnorm(train_local[[tar]], mean = rep(betahat, nrow(train_local)), 
                                        sigma = Kobs, log = TRUE)
                     }, timeout = 3600, onTimeout = "error") 
                   }, error = function(e) {
                     message("logLH computation failed or timed out: ", e$message)
                     NaN
                   })
                   
                 }

print("finished grid search")
print(c('Number of NaN vals: ',sum(is.na(logLH))))


mm <- which(logLH == max(logLH[!is.na(logLH)]))[1]
initial_par<<-grid[mm,]%>%unname
if(sum(is.na(initial_par))>0){initial_par=runif(length(initial_par))}

stopCluster(cl)

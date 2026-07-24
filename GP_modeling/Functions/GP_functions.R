source('GP_modeling/Functions/Start.R')

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
              Kmat = as.matrix(Kmat_meteo*Kmat_geno)
            )
          )
        }else{
          if(additive_only){
            return(
              list(
                Kmat = as.matrix(w1*Kmat_meteo+(1-w1)*Kmat_geno),
                del_w1=as.matrix(Kmat_meteo-Kmat_geno)
              ))
          }else{
            return(
              list(
                Kmat = as.matrix(w1*(1-w2)*Kmat_meteo+w1*w2*Kmat_geno+(1-w1)*Kmat_meteo*Kmat_geno),
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
              del_theta2 = as.matrix(Kmat_geno_deriv* Kmat_meteo)
            )
          )
        }else{
          if(additive_only){
            return(
              list(
                Kmat = as.matrix(w1*Kmat_meteo+(1-w1)*Kmat_geno),
                del_theta2 = as.matrix((1-w1)*Kmat_geno_deriv),
                del_w1=as.matrix(Kmat_meteo-Kmat_geno)
              )
            )
          }else{
            return(
              list(
                Kmat = as.matrix(w1*(1-w2)*Kmat_meteo+w1*w2*Kmat_geno+(1-w1)*Kmat_meteo*Kmat_geno),
                del_theta2 = as.matrix(w1*(w2)*Kmat_geno_deriv+(1-w1)*Kmat_meteo*Kmat_geno_deriv),
                del_w1=as.matrix((1-w2)*Kmat_meteo+w2*Kmat_geno-Kmat_meteo*Kmat_geno),
                del_w2=as.matrix(-w1*Kmat_meteo+w1*Kmat_geno)
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



################Continuous likelihood function##################################



log_likelihood = function(kernel_params,
                          alpha,
                          train) {
  
  Kmat1 <- create_kernelmat_derivative(kernel_meteo,
                                       kernel_geno,
                                       kernel_params,
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
     # message("Cholesky decomposition failed. Falling back to generalized inverse (ginv).")
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



likelihood_gradient = function(kernel_params,
                               alpha,
                               train) {
  Trace = function(x) {
    sum(diag(x))
  }
  
  kernel_matrices  <- create_kernelmat_derivative(kernel_meteo,
                                                  kernel_geno,
                                                  kernel_params,
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
     # message("Cholesky decomposition failed. Falling back to generalized inverse (ginv).")
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
  
  if(kernel_meteo=='ali'){
    if (kernel_geno == 'spe') {
      
      if(product_only){
        
        
        
        
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
           
          Trace(Kobs_inv %*% as.matrix(v * (
            Kmat1 - diag(nrow(train))
          ) + del_v_alpha * Kobs / v)) -
            2 * del_alpha_betahat %*% rep(1, nrow(train)) %*% Kobs_inv %*% (zn -
                                                                              betahat) +
            t(zn - betahat) %*% Kobs_inv_deriv_alpha %*% (zn - betahat)
        )
        
        
      }
      
      if(additive_only){
        
        
             
        
        #########################################################################################################
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(alpha * kernel_matrices$del_w1) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_w1_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        #v=as.numeric(1/nrow(train)*t(zn-betahat)%*%Kobs_inv%*%(zn-betahat))
        
        del_v_w1 = -2 * del_w1_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                             (zn - betahat)) +
          as.numeric(1 / nrow(train) * t(zn - betahat) %*% Kobs_inv_deriv %*%
                       (zn - betahat)) %>% as.numeric()
        
        
        Kobs_inv_deriv_w1 = as.matrix(1 / v * Kobs_inv_deriv - 1 / 
                                        (v ^2) * del_v_w1 * (Kobs_inv * v))
        
        
        ###########################################################################################################
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(Kmat1 - diag(nrow(train))) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_alpha_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        
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
            Kobs_inv %*% as.matrix(v * alpha * kernel_matrices$del_w1 + del_v_w1 *
                                     Kobs / v)
          ) -
            2 * del_w1_betahat %*% rep(1, nrow(train)) %*% (Kobs_inv %*% (zn -
                                                                            betahat)) -
            t(zn - betahat) %*% Kobs_inv_deriv_w1 %*% (zn - betahat),
          
          
          Trace(Kobs_inv %*% as.matrix(v * (
            Kmat1 - diag(nrow(train))
          ) + del_v_alpha * Kobs / v)) -
            2 * del_alpha_betahat %*% rep(1, nrow(train)) %*% Kobs_inv %*% (zn -
                                                                              betahat) +
            t(zn - betahat) %*% Kobs_inv_deriv_alpha %*% (zn - betahat)
        )
        
        
      }
      
      if(additive_only==0&&product_only==0){
        
        
        
        #########################################################################################################
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(alpha * kernel_matrices$del_w1) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_w1_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        #v=as.numeric(1/nrow(train)*t(zn-betahat)%*%Kobs_inv%*%(zn-betahat))
        
        del_v_w1 = -2 * del_w1_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                             (zn - betahat)) +
          as.numeric(1 / nrow(train) * t(zn - betahat) %*% Kobs_inv_deriv %*%
                       (zn - betahat)) %>% as.numeric()
        
        
        Kobs_inv_deriv_w1 = as.matrix(1 / v * Kobs_inv_deriv - 1 / 
                                        (v ^2) * del_v_w1 * (Kobs_inv * v))
        
        
        #########################################################################################################
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(alpha * kernel_matrices$del_w2) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_w2_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        #v=as.numeric(1/nrow(train)*t(zn-betahat)%*%Kobs_inv%*%(zn-betahat))
        
        del_v_w2 = -2 * del_w2_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                             (zn - betahat)) +
          as.numeric(1 / nrow(train) * t(zn - betahat) %*% Kobs_inv_deriv %*%
                       (zn - betahat)) %>% as.numeric()
        
        
        Kobs_inv_deriv_w2 = as.matrix(1 / v * Kobs_inv_deriv - 1 /
                                        (v ^2) * del_v_w2 * (Kobs_inv * v))
        
        
        
        
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
            Kobs_inv %*% as.matrix(v * alpha * kernel_matrices$del_w1 + del_v_w1 *
                                     Kobs / v)
          ) -
            2 * del_w1_betahat %*% rep(1, nrow(train)) %*% (Kobs_inv %*% (zn -
                                                                            betahat)) -
            t(zn - betahat) %*% Kobs_inv_deriv_w1 %*% (zn - betahat),
          
          Trace(
            Kobs_inv %*% as.matrix(v * alpha * kernel_matrices$del_w2 + del_v_w2 *
                                     Kobs / v)
          ) -
            2 * del_w2_betahat %*% rep(1, nrow(train)) %*% (Kobs_inv %*% (zn -
                                                                            betahat)) -
            t(zn - betahat) %*% Kobs_inv_deriv_w2 %*% (zn - betahat),
          
          
          Trace(Kobs_inv %*% as.matrix(v * (
            Kmat1 - diag(nrow(train))
          ) + del_v_alpha * Kobs / v)) -
            2 * del_alpha_betahat %*% rep(1, nrow(train)) %*% Kobs_inv %*%
            (zn -betahat) +
            t(zn - betahat) %*% Kobs_inv_deriv_alpha %*% (zn - betahat)
        )
        
        
      }
      
    } else{
      
      if(additive_only){
        
        
       
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(alpha * kernel_matrices$del_theta2) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_theta2_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        #v=as.numeric(1/nrow(train)*t(zn-betahat)%*%Kobs_inv%*%(zn-betahat))
        
        del_v_theta2 = -2 * del_theta2_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                                     (zn - betahat)) +
          as.numeric(1 / nrow(train) * t(zn - betahat) %*% Kobs_inv_deriv %*%
                       (zn - betahat)) %>% as.numeric()
        
        
        Kobs_inv_deriv_theta2 = as.matrix(1 / v * Kobs_inv_deriv - 1 / (v ^
                                                                          2) * del_v_theta2 * (Kobs_inv * v))
        
        #########################################################################################################
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(alpha * kernel_matrices$del_w1) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_w1_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        #v=as.numeric(1/nrow(train)*t(zn-betahat)%*%Kobs_inv%*%(zn-betahat))
        
        del_v_w1 = -2 * del_w1_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                             (zn - betahat)) +
          as.numeric(1 / nrow(train) * t(zn - betahat) %*% Kobs_inv_deriv %*%
                       (zn - betahat)) %>% as.numeric()
        
        
        Kobs_inv_deriv_w1 = as.matrix(1 / v * Kobs_inv_deriv - 1 / (v ^2) * 
                                        del_v_w1 * (Kobs_inv * v))
        
        
        ###########################################################################################################
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(Kmat1 - diag(nrow(train))) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_alpha_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        
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
            Kobs_inv %*% as.matrix(v * alpha * kernel_matrices$del_theta2 + del_v_theta2 *
                                     Kobs / v)
          ) -
            2 * del_theta2_betahat %*% rep(1, nrow(train)) %*% (Kobs_inv %*% (zn -
                                                                                betahat)) -
            t(zn - betahat) %*% Kobs_inv_deriv_theta2 %*% (zn - betahat),
          
          Trace(
            Kobs_inv %*% as.matrix(v * alpha * kernel_matrices$del_w1 + del_v_w1 *
                                     Kobs / v)
          ) -
            2 * del_w1_betahat %*% rep(1, nrow(train)) %*% (Kobs_inv %*% (zn -
                                                                            betahat)) -
            t(zn - betahat) %*% Kobs_inv_deriv_w1 %*% (zn - betahat),
          
          
          Trace(Kobs_inv %*% as.matrix(v * (
            Kmat1 - diag(nrow(train))
          ) + del_v_alpha * Kobs / v)) -
            2 * del_alpha_betahat %*% rep(1, nrow(train)) %*% Kobs_inv %*% (zn -
                                                                              betahat) +
            t(zn - betahat) %*% Kobs_inv_deriv_alpha %*% (zn - betahat)
          
          
          
        )
        
        
      }
      
      if(product_only){
       
        #########################################################################################################
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(alpha * kernel_matrices$del_theta2) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_theta2_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        #v=as.numeric(1/nrow(train)*t(zn-betahat)%*%Kobs_inv%*%(zn-betahat))
        
        del_v_theta2 = -2 * del_theta2_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                                     (zn - betahat)) +
          as.numeric(1 / nrow(train) * t(zn - betahat) %*% Kobs_inv_deriv %*%
                       (zn - betahat)) %>% as.numeric()
        
        
        Kobs_inv_deriv_theta2 = as.matrix(1 / v * Kobs_inv_deriv - 1 / (v ^
                                                                          2) * del_v_theta2 * (Kobs_inv * v))
        
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
            Kobs_inv %*% as.matrix(v * alpha * kernel_matrices$del_theta2 + del_v_theta2 *
                                     Kobs / v)
          ) -
            2 * del_theta2_betahat %*% rep(1, nrow(train)) %*% (Kobs_inv %*% (zn -
                                                                                betahat)) -
            t(zn - betahat) %*% Kobs_inv_deriv_theta2 %*% (zn - betahat),
          
          
          Trace(Kobs_inv %*% as.matrix(v * (
            Kmat1 - diag(nrow(train))
          ) + del_v_alpha * Kobs / v)) -
            2 * del_alpha_betahat %*% rep(1, nrow(train)) %*% Kobs_inv %*% (zn -
                                                                              betahat) +
            t(zn - betahat) %*% Kobs_inv_deriv_alpha %*% (zn - betahat)
          
          
          
        )
        
        
      }
      if(additive_only==0&&product_only==0){
      
          #########################################################################################################
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(alpha * kernel_matrices$del_theta2) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_theta2_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        #v=as.numeric(1/nrow(train)*t(zn-betahat)%*%Kobs_inv%*%(zn-betahat))
        
        del_v_theta2 = -2 * del_theta2_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                                     (zn - betahat)) +
          as.numeric(1 / nrow(train) * t(zn - betahat) %*% Kobs_inv_deriv %*%
                       (zn - betahat)) %>% as.numeric()
        
        
        Kobs_inv_deriv_theta2 = as.matrix(1 / v * Kobs_inv_deriv - 1 / (v ^
                                                                          2) * del_v_theta2 * (Kobs_inv * v))
        
        #########################################################################################################
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(alpha * kernel_matrices$del_w1) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_w1_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        #v=as.numeric(1/nrow(train)*t(zn-betahat)%*%Kobs_inv%*%(zn-betahat))
        
        del_v_w1 = -2 * del_w1_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                             (zn - betahat)) +
          as.numeric(1 / nrow(train) * t(zn - betahat) %*% Kobs_inv_deriv %*%
                       (zn - betahat)) %>% as.numeric()
        
        
        Kobs_inv_deriv_w1 = as.matrix(1 / v * Kobs_inv_deriv - 1 / (v ^2) * 
                                        del_v_w1 * (Kobs_inv * v))
        
        
        #########################################################################################################
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(alpha * kernel_matrices$del_w2) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_w2_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        #v=as.numeric(1/nrow(train)*t(zn-betahat)%*%Kobs_inv%*%(zn-betahat))
        
        del_v_w2 = -2 * del_w2_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                             (zn - betahat)) +
          as.numeric(1 / nrow(train) * t(zn - betahat) %*% Kobs_inv_deriv %*%
                       (zn - betahat)) %>% as.numeric()
        
        
        Kobs_inv_deriv_w2 = as.matrix(1 / v * Kobs_inv_deriv - 1 / (v ^ 2)
                                      * del_v_w2 * (Kobs_inv * v))
        
        
        
        
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
            Kobs_inv %*% as.matrix(v * alpha * kernel_matrices$del_theta2 + del_v_theta2 *
                                     Kobs / v)
          ) -
            2 * del_theta2_betahat %*% rep(1, nrow(train)) %*% (Kobs_inv %*% (zn -
                                                                                betahat)) -
            t(zn - betahat) %*% Kobs_inv_deriv_theta2 %*% (zn - betahat),
          
          Trace(
            Kobs_inv %*% as.matrix(v * alpha * kernel_matrices$del_w1 + del_v_w1 *
                                     Kobs / v)
          ) -
            2 * del_w1_betahat %*% rep(1, nrow(train)) %*% (Kobs_inv %*% (zn -
                                                                            betahat)) -
            t(zn - betahat) %*% Kobs_inv_deriv_w1 %*% (zn - betahat),
          
          Trace(
            Kobs_inv %*% as.matrix(v * alpha * kernel_matrices$del_w2 + del_v_w2 *
                                     Kobs / v)
          ) -
            2 * del_w2_betahat %*% rep(1, nrow(train)) %*% (Kobs_inv %*% (zn -
                                                                            betahat)) -
            t(zn - betahat) %*% Kobs_inv_deriv_w2 %*% (zn - betahat),
          
          
          Trace(Kobs_inv %*% as.matrix(v * (
            Kmat1 - diag(nrow(train))
          ) + del_v_alpha * Kobs / v)) -
            2 * del_alpha_betahat %*% rep(1, nrow(train)) %*% Kobs_inv %*% (zn -
                                                                              betahat) +
            t(zn - betahat) %*% Kobs_inv_deriv_alpha %*% (zn - betahat)
          
          
          
        )
        
        
      }
    }
  }else{
    if (kernel_geno == 'spe') {
      
      if(product_only){
        
        
        #########################################################################################################
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(alpha * kernel_matrices$del_theta1) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_theta1_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        #v=as.numeric(1/nrow(train)*t(zn-betahat)%*%Kobs_inv%*%(zn-betahat))
        
        del_v_theta1 = -2 * del_theta1_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                                     (zn - betahat)) +
          as.numeric(1 / nrow(train) * t(zn - betahat) %*% Kobs_inv_deriv %*%
                       (zn - betahat)) %>% as.numeric()
        
        
        Kobs_inv_deriv_theta1 = as.matrix(1 / v * Kobs_inv_deriv - 1 /
                                            (v ^2) * del_v_theta1 * (Kobs_inv * v))
        
        
        
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
        
        
      }
      
      if(additive_only){
        
        
        #########################################################################################################
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(alpha * kernel_matrices$del_theta1) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_theta1_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        
        del_v_theta1 = -2 * del_theta1_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                                     (zn - betahat)) +
          as.numeric(1 / nrow(train) * t(zn - betahat) %*% Kobs_inv_deriv %*%
                       (zn - betahat)) %>% as.numeric()
        
        
        Kobs_inv_deriv_theta1 = as.matrix(1 / v * Kobs_inv_deriv - 1 /
                                            (v ^2) * del_v_theta1 * (Kobs_inv * v))
        
        
        
        #########################################################################################################
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(alpha * kernel_matrices$del_w1) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_w1_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        #v=as.numeric(1/nrow(train)*t(zn-betahat)%*%Kobs_inv%*%(zn-betahat))
        
        del_v_w1 = -2 * del_w1_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                                     (zn - betahat)) +
          as.numeric(1 / nrow(train) * t(zn - betahat) %*% Kobs_inv_deriv %*%
                       (zn - betahat)) %>% as.numeric()
        
        
        Kobs_inv_deriv_w1 = as.matrix(1 / v * Kobs_inv_deriv - 1 / 
                                            (v ^2) * del_v_w1 * (Kobs_inv * v))
        
        
       ###########################################################################################################
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(Kmat1 - diag(nrow(train))) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_alpha_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        
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
          
          
          Trace(
            Kobs_inv %*% as.matrix(v * alpha * kernel_matrices$del_w1 + del_v_w1 *
                                     Kobs / v)
          ) -
            2 * del_w1_betahat %*% rep(1, nrow(train)) %*% (Kobs_inv %*% (zn -
                                                                                betahat)) -
            t(zn - betahat) %*% Kobs_inv_deriv_w1 %*% (zn - betahat),
          
          
          Trace(Kobs_inv %*% as.matrix(v * (
            Kmat1 - diag(nrow(train))
          ) + del_v_alpha * Kobs / v)) -
            2 * del_alpha_betahat %*% rep(1, nrow(train)) %*% Kobs_inv %*% (zn -
                                                                              betahat) +
            t(zn - betahat) %*% Kobs_inv_deriv_alpha %*% (zn - betahat)
        )
        
        
      }
      
      if(additive_only==0&&product_only==0){
        
        
        
        #########################################################################################################
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(alpha * kernel_matrices$del_theta1) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_theta1_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        #v=as.numeric(1/nrow(train)*t(zn-betahat)%*%Kobs_inv%*%(zn-betahat))
        
        del_v_theta1 = -2 * del_theta1_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                                     (zn - betahat)) +
          as.numeric(1 / nrow(train) * t(zn - betahat) %*% Kobs_inv_deriv %*%
                       (zn - betahat)) %>% as.numeric()
        
        
        Kobs_inv_deriv_theta1 = as.matrix(1 / v * Kobs_inv_deriv - 1 /
                                            (v ^2) * del_v_theta1 * (Kobs_inv * v))
        
        
        
        #########################################################################################################
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(alpha * kernel_matrices$del_w1) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_w1_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        #v=as.numeric(1/nrow(train)*t(zn-betahat)%*%Kobs_inv%*%(zn-betahat))
        
        del_v_w1 = -2 * del_w1_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                                     (zn - betahat)) +
          as.numeric(1 / nrow(train) * t(zn - betahat) %*% Kobs_inv_deriv %*%
                       (zn - betahat)) %>% as.numeric()
        
        
        Kobs_inv_deriv_w1 = as.matrix(1 / v * Kobs_inv_deriv - 1 / 
                                            (v ^2) * del_v_w1 * (Kobs_inv * v))
        
        
        #########################################################################################################
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(alpha * kernel_matrices$del_w2) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_w2_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        #v=as.numeric(1/nrow(train)*t(zn-betahat)%*%Kobs_inv%*%(zn-betahat))
        
        del_v_w2 = -2 * del_w2_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                                     (zn - betahat)) +
          as.numeric(1 / nrow(train) * t(zn - betahat) %*% Kobs_inv_deriv %*%
                       (zn - betahat)) %>% as.numeric()
        
        
        Kobs_inv_deriv_w2 = as.matrix(1 / v * Kobs_inv_deriv - 1 /
                                            (v ^2) * del_v_w2 * (Kobs_inv * v))
        
        
        
        
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
          
          
          Trace(
            Kobs_inv %*% as.matrix(v * alpha * kernel_matrices$del_w1 + del_v_w1 *
                                     Kobs / v)
          ) -
            2 * del_w1_betahat %*% rep(1, nrow(train)) %*% (Kobs_inv %*% (zn -
                                                                                betahat)) -
            t(zn - betahat) %*% Kobs_inv_deriv_w1 %*% (zn - betahat),
          
          Trace(
            Kobs_inv %*% as.matrix(v * alpha * kernel_matrices$del_w2 + del_v_w2 *
                                     Kobs / v)
          ) -
            2 * del_w2_betahat %*% rep(1, nrow(train)) %*% (Kobs_inv %*% (zn -
                                                                                betahat)) -
            t(zn - betahat) %*% Kobs_inv_deriv_w2 %*% (zn - betahat),
          
          
          Trace(Kobs_inv %*% as.matrix(v * (
            Kmat1 - diag(nrow(train))
          ) + del_v_alpha * Kobs / v)) -
            2 * del_alpha_betahat %*% rep(1, nrow(train)) %*% Kobs_inv %*%
            (zn -betahat) +
            t(zn - betahat) %*% Kobs_inv_deriv_alpha %*% (zn - betahat)
        )
        
        
      }
      
    } else{
    
      if(additive_only){
        
      
        ############################################################################
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(alpha * kernel_matrices$del_theta1) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_theta1_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        
        del_v_theta1 = -2 * del_theta1_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                                     (zn - betahat)) +
          as.numeric(1 / nrow(train) * t(zn - betahat) %*% Kobs_inv_deriv %*%
                       (zn - betahat)) %>% as.numeric()
        
        
        Kobs_inv_deriv_theta1 = as.matrix(1 / v * Kobs_inv_deriv - 1 / (v ^
                                                                          2) * del_v_theta1 * (Kobs_inv * v))
        
        
        #########################################################################################################
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(alpha * kernel_matrices$del_theta2) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_theta2_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        #v=as.numeric(1/nrow(train)*t(zn-betahat)%*%Kobs_inv%*%(zn-betahat))
        
        del_v_theta2 = -2 * del_theta2_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                                     (zn - betahat)) +
          as.numeric(1 / nrow(train) * t(zn - betahat) %*% Kobs_inv_deriv %*%
                       (zn - betahat)) %>% as.numeric()
        
        
        Kobs_inv_deriv_theta2 = as.matrix(1 / v * Kobs_inv_deriv - 1 / (v ^
                                                                          2) * del_v_theta2 * (Kobs_inv * v))
        
        #########################################################################################################
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(alpha * kernel_matrices$del_w1) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_w1_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        #v=as.numeric(1/nrow(train)*t(zn-betahat)%*%Kobs_inv%*%(zn-betahat))
        
        del_v_w1 = -2 * del_w1_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                                     (zn - betahat)) +
          as.numeric(1 / nrow(train) * t(zn - betahat) %*% Kobs_inv_deriv %*%
                       (zn - betahat)) %>% as.numeric()
        
        
        Kobs_inv_deriv_w1 = as.matrix(1 / v * Kobs_inv_deriv - 1 / (v ^2) * 
                                        del_v_w1 * (Kobs_inv * v))
        
        
       ###########################################################################################################
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(Kmat1 - diag(nrow(train))) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_alpha_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        
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
          
          Trace(
            Kobs_inv %*% as.matrix(v * alpha * kernel_matrices$del_theta2 + del_v_theta2 *
                                     Kobs / v)
          ) -
            2 * del_theta2_betahat %*% rep(1, nrow(train)) %*% (Kobs_inv %*% (zn -
                                                                                betahat)) -
            t(zn - betahat) %*% Kobs_inv_deriv_theta2 %*% (zn - betahat),
          
          Trace(
            Kobs_inv %*% as.matrix(v * alpha * kernel_matrices$del_w1 + del_v_w1 *
                                     Kobs / v)
          ) -
            2 * del_w1_betahat %*% rep(1, nrow(train)) %*% (Kobs_inv %*% (zn -
                                                                                betahat)) -
            t(zn - betahat) %*% Kobs_inv_deriv_w1 %*% (zn - betahat),
          
          
          Trace(Kobs_inv %*% as.matrix(v * (
            Kmat1 - diag(nrow(train))
          ) + del_v_alpha * Kobs / v)) -
            2 * del_alpha_betahat %*% rep(1, nrow(train)) %*% Kobs_inv %*% (zn -
                                                                              betahat) +
            t(zn - betahat) %*% Kobs_inv_deriv_alpha %*% (zn - betahat)
          
          
          
        )
        
        
      }
      
      if(product_only){
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
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(alpha * kernel_matrices$del_theta2) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_theta2_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        #v=as.numeric(1/nrow(train)*t(zn-betahat)%*%Kobs_inv%*%(zn-betahat))
        
        del_v_theta2 = -2 * del_theta2_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                                     (zn - betahat)) +
          as.numeric(1 / nrow(train) * t(zn - betahat) %*% Kobs_inv_deriv %*%
                       (zn - betahat)) %>% as.numeric()
        
        
        Kobs_inv_deriv_theta2 = as.matrix(1 / v * Kobs_inv_deriv - 1 / (v ^
                                                                          2) * del_v_theta2 * (Kobs_inv * v))
        
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
          
          Trace(
            Kobs_inv %*% as.matrix(v * alpha * kernel_matrices$del_theta2 + del_v_theta2 *
                                     Kobs / v)
          ) -
            2 * del_theta2_betahat %*% rep(1, nrow(train)) %*% (Kobs_inv %*% (zn -
                                                                                betahat)) -
            t(zn - betahat) %*% Kobs_inv_deriv_theta2 %*% (zn - betahat),
          
          
          Trace(Kobs_inv %*% as.matrix(v * (
            Kmat1 - diag(nrow(train))
          ) + del_v_alpha * Kobs / v)) -
            2 * del_alpha_betahat %*% rep(1, nrow(train)) %*% Kobs_inv %*% (zn -
                                                                              betahat) +
            t(zn - betahat) %*% Kobs_inv_deriv_alpha %*% (zn - betahat)
          
          
          
        )
        
        
      }
      if(additive_only==0&&product_only==0){
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
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(alpha * kernel_matrices$del_theta2) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_theta2_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        #v=as.numeric(1/nrow(train)*t(zn-betahat)%*%Kobs_inv%*%(zn-betahat))
        
        del_v_theta2 = -2 * del_theta2_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                                     (zn - betahat)) +
          as.numeric(1 / nrow(train) * t(zn - betahat) %*% Kobs_inv_deriv %*%
                       (zn - betahat)) %>% as.numeric()
        
        
        Kobs_inv_deriv_theta2 = as.matrix(1 / v * Kobs_inv_deriv - 1 / (v ^
                                                                          2) * del_v_theta2 * (Kobs_inv * v))
        
        #########################################################################################################
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(alpha * kernel_matrices$del_w1) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_w1_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        #v=as.numeric(1/nrow(train)*t(zn-betahat)%*%Kobs_inv%*%(zn-betahat))
        
        del_v_w1 = -2 * del_w1_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                                     (zn - betahat)) +
          as.numeric(1 / nrow(train) * t(zn - betahat) %*% Kobs_inv_deriv %*%
                       (zn - betahat)) %>% as.numeric()
        
        
        Kobs_inv_deriv_w1 = as.matrix(1 / v * Kobs_inv_deriv - 1 / (v ^2) * 
                                        del_v_w1 * (Kobs_inv * v))
        
        
        #########################################################################################################
        
        Kobs_inv_deriv = -as.matrix(v * Kobs_inv) %*% as.matrix(alpha * kernel_matrices$del_w2) %*%
          as.matrix(v * Kobs_inv)
        
        
        del_w2_betahat = sum(Kobs_inv_deriv %*% zn) / sum(v * Kobs_inv) -
          sum(v * Kobs_inv %*% zn) / (sum(v * Kobs_inv) ^ 2) * sum(Kobs_inv_deriv)
        
        #v=as.numeric(1/nrow(train)*t(zn-betahat)%*%Kobs_inv%*%(zn-betahat))
        
        del_v_w2 = -2 * del_w2_betahat / nrow(train) * sum(v * Kobs_inv %*%
                                                                     (zn - betahat)) +
          as.numeric(1 / nrow(train) * t(zn - betahat) %*% Kobs_inv_deriv %*%
                       (zn - betahat)) %>% as.numeric()
        
        
        Kobs_inv_deriv_w2 = as.matrix(1 / v * Kobs_inv_deriv - 1 / (v ^ 2)
                                      * del_v_w2 * (Kobs_inv * v))
        
        
        
        
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
          
          Trace(
            Kobs_inv %*% as.matrix(v * alpha * kernel_matrices$del_theta2 + del_v_theta2 *
                                     Kobs / v)
          ) -
            2 * del_theta2_betahat %*% rep(1, nrow(train)) %*% (Kobs_inv %*% (zn -
                                                                                betahat)) -
            t(zn - betahat) %*% Kobs_inv_deriv_theta2 %*% (zn - betahat),
          
          Trace(
            Kobs_inv %*% as.matrix(v * alpha * kernel_matrices$del_w1 + del_v_w1 *
                                     Kobs / v)
          ) -
            2 * del_w1_betahat %*% rep(1, nrow(train)) %*% (Kobs_inv %*% (zn -
                                                                                betahat)) -
            t(zn - betahat) %*% Kobs_inv_deriv_w1 %*% (zn - betahat),
          
          Trace(
            Kobs_inv %*% as.matrix(v * alpha * kernel_matrices$del_w2 + del_v_w2 *
                                     Kobs / v)
          ) -
            2 * del_w2_betahat %*% rep(1, nrow(train)) %*% (Kobs_inv %*% (zn -
                                                                                betahat)) -
            t(zn - betahat) %*% Kobs_inv_deriv_w2 %*% (zn - betahat),
          
          
          Trace(Kobs_inv %*% as.matrix(v * (
            Kmat1 - diag(nrow(train))
          ) + del_v_alpha * Kobs / v)) -
            2 * del_alpha_betahat %*% rep(1, nrow(train)) %*% Kobs_inv %*% (zn -
                                                                              betahat) +
            t(zn - betahat) %*% Kobs_inv_deriv_alpha %*% (zn - betahat)
          
          
          
        )
        
        
      }
    }
  }
  
  
  
  return(grad)
  
  
}


             
                      

`%||%` <- function(a,b) if (!is.null(a) && !is.na(a) && a != "") a else b

id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID") %||% 1)

MODEL <- 3#"kerlmm"
grid_search_params <- 0
keri       <- 5
additive_only <- 0
product_only  <- 0
long_Adam  <- 1
GA_kernel  <- 0
NLOPTR     <- 0

full_data=0
                         
leak_spec <- "random"
noisy_obs <- TRUE
train_prop <- 0.8
nr_outer_split <- 1

# Optim tuning
batchsize <- 0.5

if (long_Adam) {
  lr <- Lr <- c(0.01)
  max_iter <- 1000
} else {
  lr <- Lr <- 0.1
  max_iter <- 100
}
if (NLOPTR) Lr <- 1



# Load common functions/data
source("GP_modeling/Functions/Start.R", local = TRUE)
source("GP_modeling/Functions/GP_functions.R", local = TRUE)

data <- merge(data, meteo, by = "Env")


for (setupp in 1) {
  for (leakage in c("yes", "no")[2]) {
    for(tar in c("yield", "prot")){
      
        source("GP_modeling/Functions/Load_files_LMM_BGLR.R", local = TRUE) #loads the train split and corresponding optimal parameters learned through the GP
          if(full_data){
            if(tar=='prot'){opt_par=c(0.2449263,1.065338,0.8696127,0.3976185,0.9417286)
            }else{
            opt_par=c(0.2644644 ,0.78254 ,0.8542934,0.2599605,0.9478597) 
            } #these are optimal params for the setting no, 1

            train=data
            ind_train=vv=1:nrow(data)
          }  
      
          if(MODEL=='kerlmm'){
          X <- model.matrix(~ 1, data = train)
          
          Z <- matrix(0, nrow(train), nrow(data))
          for (iii in 1:nrow(train)) Z[iii, ind_train[iii]] <- 1
          
          alpha <- opt_par[length(opt_par)]
          
          Kmat1 <- create_kernelmat_derivative(
            kernel_meteo, kernel_geno, opt_par[-length(opt_par)], data
          )$Kmat
          
          noisy_obs_local <- 0
          
          if (noisy_obs_local == FALSE) {
            Kobs <- Kmat1[vv, vv]
            alpha <- 1
            Kmat <- Kmat1 + diag(1e-4, nrow = nrow(Kmat1))
          }
          
          if (noisy_obs_local == TRUE) {
            alpha <- opt_par[length(opt_par)]
            Kobs <- alpha * Kmat1[vv, vv] + (1 - alpha) * diag(nrow(train))
            Kobs_inv <- tryCatch(
              chol2inv(chol(Kobs)),
              error = function(e) ginv(as.matrix(Kobs))
            )
            Kobs_inv <- matrix(as.numeric(Kobs_inv), nrow = nrow(Kobs_inv), ncol = ncol(Kobs_inv))
            
            betahat <- sum(Kobs_inv %*% train[[tar]]) / sum(Kobs_inv)
            v <- as.numeric(1 / nrow(train) * t(train[[tar]] - betahat) %*% Kobs_inv %*% (train[[tar]] - betahat))
            
            Kobs_inv <- 1 / v * Kobs_inv
            Kmat <- Kmat1
            Kmat[vv, vv] <- v * alpha * Kobs + diag(1e-4, nrow = nrow(Kobs))
          }
          
          kerlmm <- mixed.solve(y = train[[tar]], K = Kmat, X = X, Z = Z, method = "REML", SE = TRUE)
          
          cat('Target: ',tar,',Leakage: ',leakage, ',Setup: ', setupp, '\n')
          print( summary(kerlmm))
          }else{
            
          theta1 <- opt_par[1]
          theta2 <- opt_par[2]

          meteo_dist_eucl <- as.matrix(dist(meteo[, -1], method = "euclidean"))

          rownames(meteo_dist_eucl) <- meteo$Env
          colnames(meteo_dist_eucl) <- meteo$Env
          Dmat <- meteo_dist_eucl[data$Env, data$Env]/max(meteo_dist_eucl)
          K_m <- exp(-Dmat / theta1)
                
          int <-theta2 * max(geno_dist_hamming)  
          Dmat <- geno_dist_hamming[data$variety_name, data$variety_name]
          K_g <- as.matrix(exp(-Dmat / (int)))
          K_mg <- K_m * K_g

          if (MODEL == 1) ETA <- list(list(K = K_g, model = "RKHS"),
                                      list(K = K_m, model = "RKHS"))
          if (MODEL == 2) ETA <- list(list(K = K_mg, model = "RKHS"))
          if (MODEL == 3) ETA <- list(list(K = K_g, model = "RKHS"),
                                      list(K = K_m, model = "RKHS"),
                                      list(K = K_mg, model = "RKHS"))

          response <- data[[tar]]
          #response[-ind_train] <- NA

          bglr_model <- BGLR(
            y = response, ETA = ETA,
            nIter = 10000, burnIn = 1000, thin = 2,
            verbose = FALSE, saveAt = "bglr/"
          )

            cat('Target: ',tar,', Leakage: ',leakage, ', Setup: ', setupp, '\n')
            varG   <- bglr_model$ETA[[1]]$varU
            varEnv <- bglr_model$ETA[[2]]$varU
            varGE  <- bglr_model$ETA[[3]]$varU
            varRes <- bglr_model$varE
            
            tot_all    <- varG + varEnv + varGE + varRes
            tot_signal <- varG + varEnv + varGE
            
            tab <- data.frame(
              Component = c("G", "E", "GEI", "Residual"),
              Variance  = c(varG, varEnv, varGE, varRes),
              Prop_total = c(varG, varEnv, varGE, varRes) / tot_all,
              Prop_signal = c(varG, varEnv, varGE, NA) / c(rep(tot_signal, 3), NA)
            )
            
            print(tab)
            save(
            list = "tab",
            file = paste0(
              'Results/Results BGLR/Summary_lmm_',ifelse(full_data,1,0.8),'_',
              lr,
              '_',
              max_iter,
              '_',
              MODEL,
              '_',
              tar,
              keri,
              leakage,
              setupp,
              '_',
              id,
              '.rda'
            )
          )
          }
      
      }
    
  }
}

# Logic:
# - If MODEL in 1-3: run BGLR model 1-3 (and kerlmm if MODEL='kerlmm')
# - If grid_search_params: also run LMM model 1-3 (or 4 if MODEL=4)
# - If MODEL=4: only train LMM model 4

`%||%` <- function(a,b) if (!is.null(a) && !is.na(a) && a != "") a else b

MODEL <- Sys.getenv("MODEL") %||% "kerlmm"
if (MODEL != "kerlmm") MODEL <- as.integer(MODEL)

grid_search_params <- as.integer(Sys.getenv("GS_PARAMS") %||% 0)
id         <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID") %||% 1)
keri       <- as.integer(Sys.getenv("KERNEL") %||% 5)
additive_only <- as.integer(Sys.getenv("ADDONLY") %||% 0)
product_only  <- as.integer(Sys.getenv("PRODONLY") %||% 0)
long_Adam  <- as.integer(Sys.getenv("LONG_ADAM") %||% 0)
GA_kernel  <- as.integer(Sys.getenv("GA_KERNEL") %||% 0)
NLOPTR     <- as.integer(Sys.getenv("NLOPTR") %||% 0)

leak_spec <- "random"
noisy_obs <- TRUE
train_prop <- 0.8
nr_outer_split <- 1

# Optim tuning
batchsize <- 0.5
NLOPT_alg <- "NLOPT_LD_LBFGS"

if (long_Adam) {
  lr <- Lr <- c(0.01)
  max_iter <- 1000
} else {
  lr <- Lr <- 0.1
  max_iter <- 100
}
if (NLOPTR) Lr <- 1

num_cores <- 8
learning_rate_decay <- 1

suppressPackageStartupMessages({
  library(doParallel)
  library(foreach)
})

# Load common functions/data
source("GP_modeling/Functions/Start.R", local = TRUE)
source("GP_modeling/Functions/GP_functions.R", local = TRUE)
lr_tar <- NA

cl <- makeCluster(num_cores)
registerDoParallel(cl)

data <- merge(data, meteo, by = "Env")

for (setupp in 1:2) {
  for (leakage in c("yes", "no")) {

    mm2 <- inds <- ntr <- train <- c()

    ########################
    # BGLR / kerlmm block
    ########################
    if (MODEL != 4) {

      Res <- foreach(
        tar = c("yield", "prot"),
        leakage_iter = rep(leakage,2),    # <-- pass to workers
        setupp_iter  = rep(setupp,2),     # <-- pass to workers
        .combine = "list",
        .packages = Packages,
        .export = unique(c(
          Export,
          "data", "meteo", "noisy_obs",
          "MODEL", "grid_search_params", "id", "keri",
          "additive_only", "product_only",
          "long_Adam", "GA_kernel", "NLOPTR",
          "batchsize", "NLOPT_alg", "lr", "max_iter"
        ))
      ) %dopar% {
        tic()
        # Make variables exactly as loader expects
        leakage <- leakage_iter[1]
        setupp  <- setupp_iter[1]

        ress <- numeric(26)

        # This defines: train,test,ind_train,vv,opt_par,save_results,... (branch-specific)
        source("GP_modeling/Functions/Load_files_LMM_BGLR.R", local = TRUE)

        if (MODEL == "kerlmm") {

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

          ind_test <- (1:nrow(data))[-ind_train]
          Ztest <- matrix(0, nrow(test), nrow(data))
          for (iii in 1:nrow(test)) Ztest[iii, ind_test[iii]] <- 1

          test$pred <- as.numeric(rep(kerlmm$beta, nrow(test)) + Ztest %*% kerlmm$u)

          ress[1] <- mse(test$pred, test[[tar]])
          ress[2] <- mean(abs(test$pred - test[[tar]]))

          mean_train <- mean(train[[tar]])
          ress[3] <- mse(test[[tar]], mean_train)
          ress[4] <- mean(abs(test[[tar]] - mean_train))

        } else {

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
          response[-ind_train] <- NA

          bglr_model <- BGLR(
            y = response, ETA = ETA,
            nIter = 10000, burnIn = 1000, thin = 2,
            verbose = FALSE, saveAt = "bglr/"
          )

          test$pred <- bglr_model$yHat[-ind_train]
          test$sd   <- bglr_model$SD.yHat[-ind_train]

          ress[1] <- mse(test$pred, test[[tar]])
          ress[2] <- mean(crps_norm(test[[tar]], test$pred, test$sd))

          mean_train <- mean(train[[tar]])
          ress[3] <- mse(test[[tar]], mean_train)
          ress[4] <- mean(abs(test[[tar]] - mean_train))

          ress[10] <- mean(logs_norm(test[[tar]], test$pred, test$sd))
        }

        # variety avg
        if (leakage == "yes" || as.character(setupp) == "1") {
          pred_var <- rep(NaN, nrow(test))
          for (j in 1:nrow(test)) {
            pred_var[j] <- train %>%
              dplyr::filter(variety_name == test$variety_name[j]) %>%
              dplyr::pull(tar) %>% mean() %>% as.numeric()
          }
          loca <- 1 - is.nan(pred_var)
          ress[5] <- mse(test[[tar]][loca == 1], pred_var[loca == 1])
          ress[6] <- mean(abs(test[[tar]][loca == 1] - pred_var[loca == 1]))
        }

        # env avg
        if (leakage == "yes" || as.character(setupp) == "2") {
          pred_loc <- rep(NaN, nrow(test))
          for (j in 1:nrow(test)) {
            pred_loc[j] <- train %>% dplyr::ungroup() %>%
              dplyr::filter(Env == test$Env[j]) %>%
              dplyr::pull(tar) %>% mean() %>% as.numeric()
          }
          loca <- 1 - is.nan(pred_loc)
          ress[7] <- mse(test[[tar]][loca == 1], pred_loc[loca == 1])
          ress[8] <- mean(abs(test[[tar]][loca == 1] - pred_loc[loca == 1]))
        }
        total_time=toc()
        ress[26]=total_time

        #list(ress)  # <- IMPORTANT for flat Res under .combine="c"
        ress
      }
      print(Res)
      source("GP_modeling/Functions/Save_BGLR.R")
      save_results(Res)
    }

    ########################
    # LMM grid-search block
    ########################
    if (grid_search_params) {

      Res2 <- foreach(
        tar = c("yield", "prot"),
        leakage_iter = leakage,
        setupp_iter  = setupp,
        .combine = "c",
        .multicombine = TRUE,
        .packages = Packages,
        .export = unique(c(Export, "data", "meteo", "noisy_obs", "MODEL"))
      ) %dopar% {
        tic()
        leakage <- leakage_iter
        setupp  <- setupp_iter

        ress <- numeric(26)

        source("GP_modeling/Functions/Load_files_LMM_BGLR.R", local = TRUE)

        meteo_var <- colnames(data)[5:ncol(data)]

        if (MODEL == 1) form <- paste(tar, "~ (1 | Env) + (1 | variety_name)")
        if (MODEL == 2) form <- paste(tar, "~", paste(meteo_var, collapse = " + "), "+ (1 | variety_name)")
        if (MODEL == 3) form <- paste(tar, "~", paste(meteo_var, collapse = " + "), "+ (1 | Env) + (1 | variety_name)")
        if (MODEL == 4) form <- paste(tar, "~", paste(meteo_var, collapse = " + "),
                                      "+ (1 | Env) + (1 | variety_name) + (1 | Env:variety_name)")

        model <- lmer(as.formula(form), data = train)
        test$pred <- predict(model, newdata = test, allow.new.levels = TRUE)

        ress[1] <- mse(test$pred, test[[tar]])
        ress[2] <- mean(abs(test$pred - test[[tar]]))

        mean_train <- mean(train[[tar]])
        ress[3] <- mse(test[[tar]], mean_train)
        ress[4] <- mean(abs(test[[tar]] - mean_train))

        if (leakage == "yes" || as.character(setupp) == "1") {
          pred_var <- rep(NaN, nrow(test))
          for (j in 1:nrow(test)) {
            pred_var[j] <- train %>%
              dplyr::filter(variety_name == test$variety_name[j]) %>%
              dplyr::pull(tar) %>% mean() %>% as.numeric()
          }
          loca <- 1 - is.nan(pred_var)
          ress[5] <- mse(test[[tar]][loca == 1], pred_var[loca == 1])
          ress[6] <- mean(abs(test[[tar]][loca == 1] - pred_var[loca == 1]))
        }

        if (leakage == "yes" || as.character(setupp) == "2") {
          pred_loc <- rep(NaN, nrow(test))
          for (j in 1:nrow(test)) {
            pred_loc[j] <- train %>% dplyr::ungroup() %>%
              dplyr::filter(Env == test$Env[j]) %>%
              dplyr::pull(tar) %>% mean() %>% as.numeric()
          }
          loca <- 1 - is.nan(pred_loc)
          ress[7] <- mse(test[[tar]][loca == 1], pred_loc[loca == 1])
          ress[8] <- mean(abs(test[[tar]][loca == 1] - pred_loc[loca == 1]))
        }
        total_time=toc()
        ress[26]=total_time

        list(ress)
      }

      # if your loader's save_results() uses Res2 internally, keep it visible here
      Res2 <<- Res2
      source("GP_modeling/Functions/Save_BGLR.R")
      save_results(Res2)
    }

  }
}

stopCluster(cl)

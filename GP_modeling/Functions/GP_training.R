mm2 = inds = ntr = train = c()
nr_outer_split=1
SMALL_TEST_SUBSET=0
Export <- unique(c(Export,'Export', "tar", "leakage", "setupp"))

read_int0 <- function(key, default=0L){
  v <- Sys.getenv(key, unset = "")
  x <- suppressWarnings(as.integer(v))
  if (!nzchar(v) || is.na(x)) default else x
}

id            <- read_int0("SLURM_ARRAY_TASK_ID", 1L)
keri          <- read_int0("KERNEL", 5L)
additive_only <- read_int0("ADDONLY", 0L)
product_only  <- read_int0("PRODONLY", 0L)
long_Adam     <- read_int0("LONG_ADAM", 0L)
GA_kernel     <- read_int0("GA_KERNEL", 0L)
NLOPTR        <- read_int0("NLOPTR", 0L)
UB[2]=100
res <- matrix(NA, nr_outer_split, 26)
  
  for (i in 1:nr_outer_split) {
    tic()
    if(GA_kernel){
      if(additive_only){
         inital_params_file=paste0('Results/Results_additive_only/Results grid search/Initial_par_',
                    tar, keri,'GA', leakage, setupp, '_', id,
                    '.rda')
          if(file.exists(inital_params_file)){
            load(inital_params_file)
          }else{
            print('Perform grid search for optimal initial parameters \n')
            source('GP_modeling/Functions/Initial_params.R')
            load(inital_params_file)
          }      
      }else{
        if(product_only){
          inital_params_file=paste0('Results/Results_product_only/Results grid search/Initial_par_',
                    tar, keri,'GA', leakage, setupp, '_', id,
                    '.rda')
          if(file.exists(inital_params_file)){
            load(inital_params_file)
          }else{
            print('Perform grid search for optimal initial parameters \n')
            source('GP_modeling/Functions/Initial_params.R')
            load(inital_params_file)
          }
        }else{
          inital_params_file=paste0('Results/Results grid search/Initial_par_',
                      tar, keri,'GA', leakage, setupp, '_', id,
                      '.rda')
          if(file.exists(inital_params_file)){
            load(inital_params_file)
          }else{
            print('Perform grid search for optimal initial parameters \n')
            source('GP_modeling/Functions/Initial_params.R')
            load(inital_params_file)
          }
        }
        
      }
      
      Kernel_GA = Initial_res[[9]]
      Kernel_GA_deriv=Initial_res[[10]]
    }else{
      if(additive_only){
         inital_params_file=paste0('Results/Results_additive_only/Results grid search/Initial_par_',
                    tar, keri,leakage, setupp, '_', id,
                    '.rda')
          if(file.exists(inital_params_file)){
            load(inital_params_file)
          }else{
            print('Perform grid search for optimal initial parameters \n')
            source('GP_modeling/Functions/Initial_params.R')
            load(inital_params_file)
          }      
      }else{
        if(product_only){
          inital_params_file=paste0('Results/Results_product_only/Results grid search/Initial_par_',
                    tar, keri, leakage, setupp, '_', id,
                    '.rda')
          if(file.exists(inital_params_file)){
            load(inital_params_file)
          }else{
            print('Perform grid search for optimal initial parameters \n')
            source('GP_modeling/Functions/Initial_params.R')
            load(inital_params_file)
          }
        }else{
          inital_params_file=paste0('Results/Results grid search/Initial_par_',
                      tar, keri, leakage, setupp, '_', id,
                      '.rda')
          if(file.exists(inital_params_file)){
            load(inital_params_file)
          }else{
            print('Perform grid search for optimal initial parameters \n')
            source('GP_modeling/Functions/Initial_params.R')
            load(inital_params_file)
          }
        }
        
      }
      
    }

    if(SMALL_TEST_SUBSET){
      Initial_params = Initial_res
      ind_train = Initial_params[[1]][1:10]
      train <- data[ind_train, ]
      test <- data[-ind_train, ][1:100,]
      vv <- ind_train
      ntr = nrow(train)
      }else{
        Initial_params = Initial_res
        ind_train = Initial_params[[1]]
        train <- data[ind_train, ]
        test <- data[-ind_train, ]
        vv <- ind_train
        ntr = nrow(train)
      }
    
    
    initial_par = Initial_params[[2]]
    
    
    if(  kernel_geno == 'spe'){
      if(GA_kernel){
        if (NLOPTR) {
          opt_par = nloptr(initial_par[-c(1, 2)],
                           function(x) {
                             set.seed(sum(x))
                             inds = matrix(sample(1:nrow(train)), ncol = 1 / batchsize)
                             apply(
                               foreach(
                                 ii = 1:(1 / batchsize),
                                 .combine = 'cbind',
                                 .export = c(
                                   Export
                                 ),
                                 .packages = Packages
                               ) %do% {
                                 log_likelihood(c(initial_par[1:2], x[-length(x)]), x[length(x)],
                                                train[inds[, ii],])
                               },
                               1,
                               mean
                             )
                           },
                           function(x) {
                             set.seed(sum(x))
                             inds = matrix(sample(1:nrow(train)), ncol = 1 / batchsize)
                             apply(
                               foreach(
                                 ii = 1:(1 / batchsize),
                                 .combine = 'cbind',
                                 .export = c(
                                   Export
                                 ),
                                 .packages = Packages
                               ) %do% {
                                 likelihood_gradient(c(initial_par[1:2], x[-length(x)]), x[length(x)],
                                                     train[inds[, ii],])
                               },
                               1,
                               mean
                             )
                             
                           }, opts = list(algorithm = NLOPT_alg, print_level = 1),
                           lb = LB[-c(1:2)],
                           ub = UB[-c(1:2)])$solution
          
          
        }else{
          opt_par = Adam(function(x) {
            log_likelihood(c(initial_par[1:2], x[-length(x)]), x[length(x)],
                           train)
          },function(x) {
            inds = matrix(sample(1:nrow(train)), ncol = 1 / batchsize)
            apply(
              foreach(
                ii = 1:(1 / batchsize),
                .combine = 'cbind',
                .export = c(
                  Export
                ),
                .packages = Packages
              ) %do% {
                likelihood_gradient(c(initial_par[1:2], x[-length(x)]), x[length(x)],
                                    train[inds[, ii],])
              },
              1,
              mean
            )
            
          },initial_par[-c(1, 2)], lr = lr, max_iter = max_iter, 
          beta1 = .6, beta2 = .95,
          lb = LB[-c(1:2)],
          ub = UB[-c(1:2)])
                    
        }
        opt_par=c(initial_par[1:2],opt_par)
      }else{
        if (NLOPTR) {
          opt_par = nloptr(initial_par[-2],
                           function(x) {
                             set.seed(sum(x))
                             inds = matrix(sample(1:nrow(train)), ncol = 1 / batchsize)
                             apply(
                               foreach(
                                 ii = 1:(1 / batchsize),
                                 .combine = 'cbind',
                                 .export = c(
                                   Export
                                 ),
                                 .packages = Packages
                               ) %do% {
                                 log_likelihood(c(x[1],initial_par[2],
                                                  x[-c(1,length(x))]),
                                                  x[length(x)],
                                                train[inds[, ii],])
                               },
                               1,
                               mean
                             )
                           },
                           function(x) {
                             set.seed(sum(x))
                             inds = matrix(sample(1:nrow(train)), ncol = 1 / batchsize)
                             apply(
                               foreach(
                                 ii = 1:(1 / batchsize),
                                 .combine = 'cbind',
                                 .export = c(
                                   Export
                                 ),
                                 .packages = Packages
                               ) %do% {
                                 likelihood_gradient(c(x[1],initial_par[2],
                                                  x[-c(1,length(x))]),
                                                  x[length(x)],
                                                train[inds[, ii],])
                               },
                               1,
                               mean
                             )
                             
                           }, opts = list(algorithm = NLOPT_alg, print_level = 1),
                           lb = LB[-c(2)],
                           ub = UB[-c(2)])$solution
          
          
        }else{
           opt_par = Adam(function(x) {
             log_likelihood(c(x[1],initial_par[2],
                              x[-c(1,length(x))]),
                            x[length(x)],
                            train)
          },function(x) {
            inds = matrix(sample(1:nrow(train)), ncol = 1 / batchsize)
            apply(
            foreach(
              ii = 1:(1 / batchsize),
              .combine = 'cbind',
              .export = c(
                Export
              ),
              .packages = Packages
            ) %do% {
              likelihood_gradient(c(x[1],initial_par[2],
                                    x[-c(1,length(x))]),
                                    x[length(x)],
                                  train[inds[, ii],])
            },
            1,
            mean
          )
          
        },initial_par[-c(2)], lr = lr, max_iter = max_iter, 
          beta1 = .6, beta2 = .95,
          lb = LB[-c(2)],
          ub = UB[-c(2)])
           
          
          
          
        }
        opt_par=c(opt_par[1],initial_par[2],opt_par[2:length(opt_par)])
        }
    }else{
      if(GA_kernel){
        if (NLOPTR) {
        opt_par = nloptr(initial_par[-c(1)],
                         function(x) {
                           set.seed(sum(x))
                           inds = matrix(sample(1:nrow(train)), ncol = 1 / batchsize)
                           apply(
                             foreach(
                               ii = 1:(1 / batchsize),
                               .combine = 'cbind',
                               .export = c(
                                 Export
                               ),
                               .packages = Packages
                             ) %do% {
                               log_likelihood(c(initial_par[1],x[-length(x)]), x[length(x)],
                                              train[inds[, ii],])
                             },
                             1,
                             mean
                           )
                         },
                         function(x) {
                           set.seed(sum(x))
                           inds = matrix(sample(1:nrow(train)), ncol = 1 / batchsize)
                           apply(
                             foreach(
                               ii = 1:(1 / batchsize),
                               .combine = 'cbind',
                               .export = c(
                                 Export
                               ),
                               .packages = Packages
                             ) %do% {
                               likelihood_gradient(c(initial_par[1],x[-length(x)]), x[length(x)],
                                                   train[inds[, ii],])
                             },
                             1,
                             mean
                           )
                           
                         }, opts = list(algorithm = NLOPT_alg, print_level = 1),
                         lb = LB[-c(1)],
                         ub = UB[-c(1)])$solution
        
        
      }else{
        opt_par = Adam(function(x) {
          log_likelihood(c(initial_par[1],x[-length(x)]), x[length(x)],
                         train)
        },function(x) {
          inds = matrix(sample(1:nrow(train)), ncol = 1 / batchsize)
          apply(
            foreach(
              ii = 1:(1 / batchsize),
              .combine = 'cbind',
              .export = c(
                Export
              ),
              .packages = Packages
            ) %do% {
              likelihood_gradient(c(initial_par[1],x[-length(x)]), x[length(x)],
                                  train[inds[, ii],])
            },
            1,
            mean
          )
          
        },initial_par[-c(1)], lr = lr, max_iter = max_iter, 
        beta1 = .6, beta2 = .95,
        lb = LB[-c(1)],
        ub = UB[-c(1)])
        
      }
      opt_par=c(initial_par[1],opt_par)
    }else{
      if (NLOPTR) {
        opt_par = nloptr(initial_par,
                         function(x) {
                           set.seed(sum(x))
                           inds = matrix(sample(1:nrow(train)), ncol = 1 / batchsize)
                           apply(
                             foreach(
                               ii = 1:(1 / batchsize),
                               .combine = 'cbind',
                               .export = c(
                                 Export
                               ),
                               .packages = Packages
                             ) %do% {
                               log_likelihood(x[-length(x)], x[length(x)],
                                              train[inds[, ii],])
                             },
                             1,
                             mean
                           )
                         },
                         function(x) {
                           set.seed(sum(x))
                           inds = matrix(sample(1:nrow(train)), ncol = 1 / batchsize)
                           apply(
                             foreach(
                               ii = 1:(1 / batchsize),
                               .combine = 'cbind',
                               .export = c(
                                 Export
                               ),
                               .packages = Packages
                             ) %do% {
                               likelihood_gradient( x[-length(x)], x[length(x)],
                                                    train[inds[, ii],])
                             },
                             1,
                             mean
                           )
                           
                         }, opts = list(algorithm = NLOPT_alg, print_level = 1),
                         lb = LB,
                         ub = UB)$solution
        
        
      }else{
        opt_par = Adam(function(x) {
          log_likelihood(x[-length(x)], x[length(x)],
                         train)
        },function(x) {
          inds = matrix(sample(1:nrow(train)), ncol = 1 / batchsize)
          apply(
            foreach(
              ii = 1:(1 / batchsize),
              .combine = 'cbind',
              .export = c(
                Export
              ),
              .packages = Packages
            ) %do% {
              likelihood_gradient(x[-length(x)], x[length(x)],
                                  train[inds[, ii],])
            },
            1,
            mean
          )
          
        },initial_par, lr = lr, max_iter = max_iter, 
        beta1 = .6, beta2 = .95,
        lb = LB,
        ub = UB)
        
      }
      
    }
      }
      optim_time=toc()
      tic()
      alpha=opt_par[length(opt_par)]
      
      Kmat1 <-create_kernelmat_derivative(kernel_meteo, 
                                          kernel_geno, 
                                          opt_par[-length(opt_par)],
                                          data)$Kmat
      
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
      
        res[i, 9] = log_likelihood(opt_par[-length(opt_par)],
                                   opt_par[length(opt_par)],
                                   train)
        res[i,10]=mean(logs_norm(test[[tar]], test$pred, test$sd)) #mean(test$LogS)
        res[i,11]=LogS(data, train, test, betahat, Kobs_inv,
                       Kmat, 1:nrow(test), vv)
        res[i, 13:(13+length(opt_par)-1)] = opt_par
        opt_par = initial_par
        res[i, 12] = log_likelihood(opt_par[-length(opt_par)],
                                    opt_par[length(opt_par)],
                                    train)
        res[i, 19:(19+length(opt_par)-1)] = opt_par
        inference_time=toc()
        res[i, 24]=optim_time
        res[i, 25]=inference_time
        res[i, 26]=optim_time+inference_time 
               
      
      
  }

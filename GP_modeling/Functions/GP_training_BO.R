
nr_outer_split=1
SMALL_TEST_SUBSET=0
Export <- unique(c(Export,'Export', "tar", "leakage", "setupp"))

UB[2]=100

if(params_opt){
  tic()
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
        if(batchsize!=1){
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
            #apply(sapply(1:(1 / batchsize),
            #             function(ii){ 
            #               likelihood_gradient(
            #                 x[-length(x)], x[length(x)],
            #                 train[inds[, ii],])}),
            #     1,mean)
            
          },initial_par, lr = lr, max_iter = max_iter, 
          beta1 = .6, beta2 = .95,
          lb = LB,
          ub = UB)
          
        }else{
          opt_par = Adam(function(x) {
            log_likelihood(x[-length(x)], x[length(x)],
                           train)
          },function(x) {
            likelihood_gradient(x[-length(x)], x[length(x)],
                                train)
          }
          
          ,initial_par, lr = lr, max_iter = max_iter, 
          beta1 = .6, beta2 = .95,
          lb = LB,
          ub = UB)
        }
      }
    }
  }
  optim_time=toc()
  print(optim_time)
}
opt_par<<-opt_par


alpha<-opt_par[length(opt_par)]

theta1=opt_par[1]
theta2=opt_par[2]
w1 <- opt_par[3]
w2 <- opt_par[4]


# Meteo kernel
meteo_dist_eucl <- as.matrix(stats::dist(meteo[, -1, drop = FALSE], method = "euclidean"))
rownames(meteo_dist_eucl) <- meteo$Env
colnames(meteo_dist_eucl) <- meteo$Env
Dmat_m <- meteo_dist_eucl[data$Env, data$Env, drop = FALSE] / max(meteo_dist_eucl)
K_m <- exp(-(Dmat_m^2) / (2 * theta1^2))

# Genomic kernel
Dmat_g <- D_gblup[data$variety_name, data$variety_name] / max(D_gblup)
K_g <- exp(-(Dmat_g^2) / (2 * theta2^2))
K_mg <- K_m * K_g

# Combined kernel (constant across reps)
Kmat1 <- as.matrix(w1 * (1 - w2) * K_m + w1 * w2 * K_g + (1 - w1) * K_mg)

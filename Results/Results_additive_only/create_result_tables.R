product_only=0
for(keri in 1:9){
  for(optimizer in 1:3){
    for(leakage in c('yes', 'no')){
      for(setup in 1:2){
        for(tar in c('yield','prot')){
          optimizers=c('Results grid search/Initial_par',
                       'Results LBFGS/Results_hyper_NLOPT_LD_LBFGS_0.5',
                      # 'Results Adam/Results_hyper_ADAM_0.1_100_0.5',
                       'Results Adam/Results_hyper_ADAM_0.001_1000_0.5')
          
          tryCatch({
          results <- vector("list", length(l))
          for (id in l) {
            if(keri<7){
            file_name <- paste0(
              optimizers[optimizer],
              '_',
              tar,
              keri,
              leakage,
              setup,
              "_",
              id,
              ".rda")
            }else{
               file_name <- paste0(
              optimizers[optimizer],
              '_',
              tar,
              keri,'GA',
              leakage,
              setup,
              "_",
              id,
              ".rda")
            }
            
            if (file.exists(file_name)) {
              load(file_name)
              if(optimizer==1){
                results[[id]] <- c(Initial_res[[4]],Initial_res[[5]],Initial_res[[6]])
              }else{
                results[[id]] <- res[c(1,2,10)]
              }
            } else {
              cat("File not found:", file_name, "\n")
            }
          }
          results <- matrix(unlist(results),ncol=3,byrow=1)
          colnames(results)=c('MSE','CRPS','LogS')
          
          save(list = "results", 
               file = paste0('Result_tables/Result_table_',
                             c('Grid_search','LBFGS',
                               'Adam_0.1_100',
                               'Adam_0.001_1000')[optimizer],
                             '_',
                             tar,'_', keri, leakage,setup,'.rda'))
            },error = function(e) {
                     message("Results not yet existing", e$message)
                     NaN
                   })
        }
      }
    }
  }
}


for(MODEL in 1:3){
  for(leakage in c('yes', 'no')){
    for(setup in 1:2){
      for(tar in c('yield','prot')){
        if(MODEL==4){
          if(product_only){
            optimizers=c('grid_search_',
                         'NLOPT_LD_LBFGS',
                        # 'ADAM_0.1_100_',
                         'ADAM_0.001_1000_')
            for(optimizer in 1:3){
              results <- vector("list", length(l))
              
              for (id in l) {
                file_name <- paste0(
                  'Results BGLR/Results_hyper_',
                  optimizers[optimizer],
                  MODEL,
                  '_',
                  tar,
                  1,
                  leakage,
                  setup,
                  '_',
                  id,
                  '.rda'
                )
                
                if (file.exists(file_name)) {
                  load(file_name)
                  results[[id]] <- res[c(1,2,10)]
                } else {
                  cat("File not found:", file_name, "\n")
                }
              }
              results <- do.call(rbind, results)
              colnames(results)=c('MSE','CRPS','LogS')
              
              save(list = "results", 
                   file = paste0('Result_tables/Result_table_BLGR_',
                                 MODEL,
                                 '_',
                                 optimizers[optimizer],
                                 '_',
                                 tar, '_',leakage,setup,'.rda'))
              
          }
        
          }else{
              break
          }
        }else{
          
        
          results <- vector("list", length(l))
          for (id in l) {
            file_name <- paste0(
              'Results LMM/Results_LMM',
              MODEL,
              '_',
              tar,
              1,
              leakage,
              setup,
              '_',
              id,
              '.rda'
            )
            if (file.exists(file_name)) {
              load(file_name)
              results[[id]] <- res[1:2]
            } else {
              cat("File not found:", file_name, "\n")
            }
          }
          results <- do.call(rbind, results)
          colnames(results)=c('MSE','CRPS')
          save(list = "results", 
               file = paste0('Result_tables/Result_table_LMM_',MODEL,
                             '_',
                             tar,'_', leakage,setup,'.rda'))
          
          
          optimizers=c('grid_search_',
                       'NLOPT_LD_LBFGS',
                     #  'ADAM_0.1_100_',
                       'ADAM_0.001_1000_')
          for(optimizer in 1:3){
            results <- vector("list", length(l))
            
            for (id in l) {
              file_name <- paste0(
                'Results BGLR/Results_hyper_',
                optimizers[optimizer],
                MODEL,
                '_',
                tar,
                1,
                leakage,
                setup,
                '_',
                id,
                '.rda'
              )
              
              if (file.exists(file_name)) {
                load(file_name)
                results[[id]] <- res[c(1,2,10)]
              } else {
                cat("File not found:", file_name, "\n")
              }
            }
            results <- do.call(rbind, results)
            colnames(results)=c('MSE','CRPS','LogS')
            
            save(list = "results", 
                 file = paste0('Result_tables/Result_table_BGLR_',
                               MODEL,
                               '_',
                               optimizers[optimizer],
                               tar, '_',leakage,setup,'.rda'))
            
            
        }
        
        }
      }
    }
  }
}

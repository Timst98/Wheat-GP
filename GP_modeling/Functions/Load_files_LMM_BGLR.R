if(additive_only){

load(paste0('Results/Results_additive_only/Results grid search/Initial_par_',
                  tar, keri, leakage, setupp, '_', id,
                  '.rda'))
      
      Initial_params = Initial_res
      ind_train = Initial_params[[1]]
      train <- data[ind_train, ]
      test <- data[-ind_train, ]
      vv <- ind_train
      ntr = nrow(train)

if(grid_search_params){
        opt_par= Initial_params[[2]]
      }else{
        if(NLOPTR){
          load(paste0(
            'Results/Results_additive_only/Results LBFGS/Results_hyper_',
            NLOPT_alg,
            '_',
            batchsize,
            '_',
            tar,
            keri,
            leakage,
            setupp,
            '_',
            id,
            '.rda'
          ))
          opt_par=res[13:14]
        }else{
          load(paste0(
            'Results/Results_additive_only/Results Adam/Results_hyper_ADAM_',
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
            id,
            '.rda'
          ))
          opt_par=res[13:16]
        }  
      }

  save_results=function(Res){
    for (i in 1:2) {
        tar = c('yield', 'prot')[i]
        res = as.numeric(Res[[i]])
        
        if(grid_search_params){
          if(MODEL!=4){
            save(
              list = "res",
              file = paste0(
                'Results/Results_additive_only/Results BGLR/Results_hyper_grid_search_',MODEL,
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
          res=as.numeric(Res2[[i]])
          save(
            list = "res",
            file = paste0(
              'Results/Results_additive_only/Results LMM/Results_LMM',
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
        }else{
          if (NLOPTR) {
            
            save(
              list = "res",
              file = paste0(
                'Results/Results_additive_only/Results BGLR/Results_hyper_',
                NLOPT_alg,
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
            
          }else{
            
            save(
              list = "res",
              file = paste0(
                'Results/Results_additive_only/Results BGLR/Results_hyper_ADAM_',
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
    cat('Saved results \n')
  }

  }else{
  if(product_only){
    load(paste0('Results/Results_product_only/Results grid search/Initial_par_',
                          tar, keri, leakage, setupp, '_', id,
                          '.rda'))
              
              Initial_params = Initial_res
              ind_train = Initial_params[[1]]
              train <- data[ind_train, ]
              test <- data[-ind_train, ]
              vv <- ind_train
              ntr = nrow(train)
              
              if(grid_search_params){
                opt_par= Initial_params[[2]]
              }else{
                if(NLOPTR){
                  load(paste0(
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
                    id,
                    '.rda'
                  ))
                  opt_par=res[13:14]
                }else{
                  load(paste0(
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
                    id,
                    '.rda'
                  ))
                  opt_par=res[13:15]
                }  
              }


    save_results=function(Res){
      for (i in 1:2) {
      tar = c('yield', 'prot')[i]
      res = as.numeric(Res[[i]])
      
      if(grid_search_params){
        if(MODEL!=4){
        save(
          list = "res",
          file = paste0(
            'Results/Results_product_only/Results BGLR/Results_hyper_grid_search_',MODEL,
            '_',
            tar,
            keri,
            leakage,
            setupp,
            '_',
            id,
            '.rda'
          )
        )}
       # if(MODEL!=4){
          res=as.numeric(Res2[[i]])
          save(
            list = "res",
            file = paste0(
              'Results/Results_product_only/Results LMM/Results_LMM',
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
        #}
      }else{
        if (NLOPTR) {
          
          save(
            list = "res",
            file = paste0(
              'Results/Results_product_only/Results BGLR/Results_hyper_',
              NLOPT_alg,
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
          
        }else{
          
          save(
            list = "res",
            file = paste0(
              'Results/Results_product_only/Results BGLR/Results_hyper_ADAM_',
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
      cat('Saved results \n')
      }
    
  }else{
    load(paste0('Results/Results grid search/Initial_par_',
                    tar, keri, leakage, setupp, '_', id,
                    '.rda'))
        
        Initial_params = Initial_res
        ind_train = Initial_params[[1]]
        train <- data[ind_train, ]
        test <- data[-ind_train, ]
        vv <- ind_train
        ntr = nrow(train)
        
        if(grid_search_params){
          opt_par= Initial_params[[2]]
        }else{
          if(NLOPTR){
            load(paste0(
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
              id,
              '.rda'
            ))
            opt_par=res[13:14]
          }else{
            load(paste0(
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
              id,
              '.rda'
            ))
            opt_par=res[13:17]
          }  
        }
     save_results=function(Res){
       for (i in 1:2) {
      tar = c('yield', 'prot')[i]
      res = as.numeric(Res[[i]])
      
      if(grid_search_params){
        if(MODEL!=4){
        save(
          list = "res",
          file = paste0(
            'Results/Results BGLR/Results_hyper_grid_search_',MODEL,
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
        res=as.numeric(Res2[[i]])
        save(
          list = "res",
          file = paste0(
            'Results/Results LMM/Results_LMM',
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
      }else{
        if (NLOPTR) {
          
          save(
            list = "res",
            file = paste0(
              'Results/Results BGLR/Results_hyper_',
              NLOPT_alg,
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
          
        }else{
          
          save(
            list = "res",
            file = paste0(
              'Results/Results BGLR/Results_hyper_ADAM_',
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
       cat('Saved results \n')
     }
  }
  
  }

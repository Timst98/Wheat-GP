if(additive_only){
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
          res=as.numeric(Res[[i]])
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
          res=as.numeric(Res[[i]])
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
        res=as.numeric(Res[[i]])
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

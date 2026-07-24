library(dplyr)
l=(1:30)

for(setup in 1:2){
  for(leakage in c('no','yes')){
    
    
    RESULTS=matrix(nrow=15,ncol=6)
    colnames(RESULTS)=c('MSE Y','CRPS Y','LogS Y',
                        'MSE P','CRPS P','LogS P')
    rownames(RESULTS)=c(
      'ker 1 GS',
      'ker 2 GS',
      'ker 3 GS',
      'ker 4 GS',
      'ker 5 GS',
      'ker 6 GS',
      'ker 7 GS',
      'ker 8 GS',
      'ker 9 GS',
      'LMM1',
      'LMM2',
      'LMM3',
      'BGLR 1 GS',
      'BGLR 2 GS',
      'BGLR 3 GS'
    )
    
    for(keri in 1:9){
      tar='yield'
      results <- vector("list", length(l))
      for (id in l) {
        if(keri<7){
          file_name <- paste0(
            "Results/Results_additive_only/Results grid search/Initial_par_",tar,keri,leakage,setup,"_",
            id, ".rda")
        }else{
          file_name <- paste0(
            "Results/Results_additive_only/Results grid search/Initial_par_",tar,keri,'GA',leakage,setup,"_",
            id, ".rda")
        }
        if (file.exists(file_name)) {
          load(file_name)
          results[[id]] <- Initial_res
          a=sapply(4:6,function(i)Initial_res[[i]])
          if(sum(is.nan(a))+sum(is.na(a))>0){
           cat("NaN values in Results: ", file_name, "\n")
         }
        } else {
          cat("File not found:", file_name, "\n")
        }
      }
      results <- do.call(rbind, results)
      RESULTS[keri,1]=(mean(results[,4]%>%as.numeric)%>%round(2))
      RESULTS[keri,2]=(mean(results[,5]%>%as.numeric)%>%round(2))
      RESULTS[keri,3]=(mean(results[,6]%>%as.numeric)%>%round(2))
      
      
      tar='prot'
      results2 <- vector("list", length(l))
      for (id in l) {
        if(keri<7){
          file_name <- paste0(
            "Results/Results_additive_only/Results grid search/Initial_par_",tar,keri,leakage,setup,"_",
            id, ".rda")
        }else{
          file_name <- paste0(
            "Results/Results_additive_only/Results grid search/Initial_par_",tar,keri,'GA',leakage,setup,"_",
            id, ".rda")
        }
        if (file.exists(file_name)) {
          load(file_name)
          results2[[id]] <-  Initial_res
          a=sapply(4:6,function(i)Initial_res[[i]])
          if(sum(is.nan(a))+sum(is.na(a))>0){
           cat("NaN values in Results: ", file_name, "\n")
         }
        } else {
          cat("File not found:", file_name, "\n")
        }
      }
      results2 <- do.call(rbind, results2)
      
      RESULTS[keri,4]=(mean(results2[,4]%>%as.numeric)%>%round(2))
      RESULTS[keri,5]=(mean(results2[,5]%>%as.numeric)%>%round(2))
      RESULTS[keri,6]=(mean(results2[,6]%>%as.numeric)%>%round(2))
      
    }
     
       
      
      keri=1
      # LMM SCORES
      for(MODEL in 1:3){
        j=0
        for(tar in c('yield','prot')){
          results <- vector("list", length(l))
          for (id in l) {
            file_name <- paste0(
              'Results/Results_additive_only/Results LMM/Results_LMM',
              MODEL,
              '_',
              tar,
              keri,
              leakage,
              setup,
              '_',
              id,
              '.rda'
            )
            if (file.exists(file_name)) {
              load(file_name)
              results[[id]] <- res
            } else {
              cat("File not found:", file_name, "\n")
            }
          }
          results <- do.call(rbind, results)
          
          RESULTS[(10)+MODEL-1,((1:2)+j)]=round(apply(results[,1:2],2,mean),2)
          
          j=j+3
        }
        
      }
      
      Alg=c('grid_search_')#,'NLOPT_LD_LBFGS',Adam_alg2)
      ind=c(13:15)
            for(MODEL in 1:3){
              j=0
              for (alg in Alg){
                jj=0
                for(tar in c('yield','prot')){
                  results <- vector("list", length(l))
                  for (id in l) {
                    file_name <- paste0(
                      'Results/Results_additive_only/Results BGLR/Results_hyper_',
                      alg,
                      MODEL,
                      '_',
                      tar,
                      keri,
                      leakage,
                      setup,
                      '_',
                      id,
                      '.rda'
                    )
                    
                    if (file.exists(file_name)) {
                      load(file_name)
                      results[[id]] <- res
                    } else {
                      cat("File not found:", file_name, "\n")
                    }
                  }
                  results <- do.call(rbind, results)
                  
                  RESULTS[ind[MODEL]+j,(1:2)+jj]=round(apply(results[,1:2],2,mean),2)
                  
                  RESULTS[ind[MODEL]+j,3+jj]=mean(results[,10]%>%as.numeric)%>%round(2)
                  
                  jj=3
                }
                j=j+1
                
              }
              
            }
            print(c('SETUP=',setup,'Leakage=',leakage))
            
            print(RESULTS)
            
            
    }
    
  }

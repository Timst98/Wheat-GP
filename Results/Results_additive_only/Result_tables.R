library(dplyr)
l=(1:30)

Adam_alg="Results/Results_additive_only/Results Adam/Results_hyper_ADAM_0.1_100_0.5"
Adam_alg="Results/Results_additive_only/Results Adam/Results_hyper_ADAM_0.001_1000_0.5"

Adam_alg2='ADAM_0.1_100_'
Adam_alg2='ADAM_0.001_1000_'

for(setup in 1:2){
  for(leakage in c('no','yes')){
    
    RESULTS=matrix(nrow=35+ifelse(leakage=='yes',1,0),ncol=6)
    colnames(RESULTS)=c('MSE Y','CRPS Y','LogS Y',
                        'MSE P','CRPS P','LogS P')
    names=c(
      'ker 1 GS',
      'ker 1 LBGFS',
      'ker 1 Adam',
      'ker 2 GS',
      'ker 2 LBGFS',
      'ker 2 Adam',
      'ker 3 GS',
      'ker 3 LBGFS',
      'ker 3 Adam',
      'ker 4 GS',
      'ker 4 LBGFS',
      'ker 4 Adam',
      'ker 5 GS',
      'ker 5 LBGFS',
      'ker 5 Adam',
      'ker 6 GS',
      'ker 6 LBGFS',
      'ker 6 Adam',
      'ker 7 GS',
      'ker 8 GS',
      'ker 9 GS',
      'Global',
      'Variety',
      'Env',
      'LMM1',
      'LMM2',
      'LMM3',
      'BGLR 1 GS',
      'BGLR 1 LBGFS',
      'BGLR 1 Adam',
      'BGLR 2 GS',
      'BGLR 2 LBGFS',
      'BGLR 2 Adam',
      'BGLR 3 GS',
      'BGLR 3 LBGFS',
      'BGLR 3 Adam'
    )
    if(leakage=='no'){
      rownames(RESULTS)=names[-24]
    }else{
      rownames(RESULTS)=names
    }
    for(keri in 1:6){
      results <- vector("list", length(l))
      for (id in l) {
        file_name <- paste0(
          
          #"Results_hyper_yield1no1_"
          'Results/Results_additive_only/Results grid search/Initial_par_yield',keri,leakage,setup,"_",
          id, ".rda")
        if (file.exists(file_name)) {
          load(file_name)
          results[[id]] <- Initial_res
        } else {
      #    cat("File not found:", file_name, "\n")
        }
      }
      results <- do.call(rbind, results)
      RESULTS[3*(keri-1)+1,1]=(mean(results[,4]%>%as.numeric)%>%round(2))
      RESULTS[3*(keri-1)+1,2]=(mean(results[,5]%>%as.numeric)%>%round(2))
      RESULTS[3*(keri-1)+1,3]=(mean(results[,6]%>%as.numeric)%>%round(2))
      
      j=2
      for(alg in c('Results/Results_additive_only/Results LBFGS/Results_hyper_NLOPT_LD_LBFGS_0.5',Adam_alg)){
        
        results <- vector("list", length(l))
        for (id in l) {
          file_name <- paste0(
            alg,"_yield",keri,leakage,setup,"_",
            id, ".rda")
          if (file.exists(file_name)) {
            load(file_name)
            results[[id]] <- res
          } else {
            cat("File not found:", file_name, "\n")
          }
        }
        results <- do.call(rbind, results)
        
        
        RESULTS[3*(keri-1)+j,1:2]=(round(apply(results[,1:2],2,mean),2))
        RESULTS[3*(keri-1)+j,3]=mean(results[,10]%>%as.numeric)%>%round(2)
        
        j=j+1
      }
      
      
      results2 <- vector("list", length(l))
      for (id in l) {
        file_name <- paste0(
          'Results/Results_additive_only/Results grid search/Initial_par_prot',keri,leakage,setup,"_",
          id, ".rda")
        if (file.exists(file_name)) {
          load(file_name)
          results2[[id]] <-  Initial_res
        } else {
      #    cat("File not found:", file_name, "\n")
        }
      }
      results2 <- do.call(rbind, results2)
      
      RESULTS[3*(keri-1)+1,4]=(mean(results2[,4]%>%as.numeric)%>%round(2))
      RESULTS[3*(keri-1)+1,5]=(mean(results2[,5]%>%as.numeric)%>%round(2))
      RESULTS[3*(keri-1)+1,6]=(mean(results2[,6]%>%as.numeric)%>%round(2))
      
      j=2
      for(alg in  c('Results/Results_additive_only/Results LBFGS/Results_hyper_NLOPT_LD_LBFGS_0.5',Adam_alg)){
        
        results <- vector("list", length(l))
        for (id in l) {
          file_name <- paste0(
            alg,"_prot",keri,leakage,setup,"_",
            id, ".rda")
          if (file.exists(file_name)) {
            load(file_name)
            results[[id]] <- res
          } else {
         #   cat("File not found:", file_name, "\n")
          }
        }
        results <- do.call(rbind, results)
        
        RESULTS[3*(keri-1)+j,4:5]=(round(apply(results[,1:2],2,mean),2))
        RESULTS[3*(keri-1)+j,6]=mean(results[,10]%>%as.numeric)%>%round(2)
        
        j=j+1
      }
      
    } 

    for(k in 1:3){
      keri=(7:9)[k]
      results <- vector("list", length(l))
      for (id in l) {
        file_name <- paste0(
          
          #"Results_hyper_yield1no1_"
          'Results/Results_additive_only/Results grid search/Initial_par_yield',keri,'GA',leakage,setup,"_",
          id, ".rda")
        if (file.exists(file_name)) {
          load(file_name)
          results[[id]] <- Initial_res
        } else {
        #  cat("File not found:", file_name, "\n")
        }
      }
      results <- do.call(rbind, results)
      RESULTS[18+k,1]=(mean(results[,4]%>%as.numeric)%>%round(2))
      RESULTS[18+k,2]=(mean(results[,5]%>%as.numeric)%>%round(2))
      RESULTS[18+k,3]=(mean(results[,6]%>%as.numeric)%>%round(2))
      
      results2 <- vector("list", length(l))
      for (id in l) {
        file_name <- paste0(
          'Results/Results_additive_only/Results grid search/Initial_par_prot',keri,'GA',leakage,setup,"_",
          id, ".rda")
        if (file.exists(file_name)) {
          load(file_name)
          results2[[id]] <-  Initial_res
        } else {
        #  cat("File not found:", file_name, "\n")
        }
      }
      results2 <- do.call(rbind, results2)
      
      RESULTS[18+k,4]=(mean(results2[,4]%>%as.numeric)%>%round(2))
      RESULTS[18+k,5]=(mean(results2[,5]%>%as.numeric)%>%round(2))
      RESULTS[18+k,6]=(mean(results2[,6]%>%as.numeric)%>%round(2))
      
    }
    ## EXTRA SCORES
    jj=0
    for (tar in c('yield','prot')){
      R=c()
      
      for(keri in 1:6){
        results <- vector("list", length(l))
        for (id in l) {
           if(keri<7){
          file_name <- paste0(
            alg,'_',tar,keri,leakage,setup,"_",
            id, ".rda")
        }else{
          file_name <- paste0(
            alg,tar,keri,'GA',leakage,setup,"_",
            id, ".rda")
        }
          if (file.exists(file_name)) {
            load(file_name)
            results[[id]] <- res
          } else {
            cat("File not found:", file_name, "\n")
          }
        }
        results <- as.matrix(do.call(rbind, results)[,1:8])
        R=rbind(R,results)
        
      }
      
      #print(round(apply(as.matrix(R),2,mean),2))
      res=round(apply(as.matrix(R),2,mean),2)
      RESULTS[22,1:2+jj]=res[3:4]
#      if (leakage == 'no' && setup == '1') {
        RESULTS[23,1:2+jj]=res[5:6]
        
  #    }
      # environmental average
      if (leakage == 'no' && setup == 2) {
        RESULTS[23,1:2+jj]=res[7:8]
        rownames(RESULTS)[23]='Env'
      }
      
      if(leakage=='yes'){
        RESULTS[24,1:2+jj]=res[7:8]
      }
      jj=3
    }
    
    
    
   
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
            1,
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
        
       RESULTS[(24+ifelse(leakage=='yes',1,0))+MODEL-1,((1:2)+j)]=round(apply(results[,1:2],2,mean),2)
       
         j=j+3
      }
     
    }
    
    Alg=c('grid_search_','NLOPT_LD_LBFGS',Adam_alg2)
    ind=c(27,30,33)+ifelse(leakage=='yes',1,0)
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
              1,
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
       #       cat("File not found:", file_name, "\n")
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
  
  
  
   



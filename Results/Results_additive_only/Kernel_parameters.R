library(dplyr)
l=(1:30)

for(setup in 1:2){
  for(leakage in c('no','yes')){
    print(c('SETUP=',setup,'Leakage=',leakage))
    
    RESULTS=matrix(nrow=12,ncol=4)
    colnames(RESULTS)=c('alpha1 Y ','alpha2 Y',
                        'alpha1 P','alpha2 P')
    rownames(RESULTS)=c(
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
      'ker 4 Adam')
    for(keri in 1:4){
      if(keri%%2==0){
        alpha_par_ind=2:3
      }else{
        alpha_par_ind=3:4
      }
      results <- vector("list", length(l))
      for (id in l) {
        file_name <- paste0(
          
          #"Results_hyper_yield1no1_"
          "Results/Results_additive_only/Results grid search/Initial_par_yield",keri,leakage,setup,"_",
          id, ".rda")
        if (file.exists(file_name)) {
          load(file_name)
          results[[id]] <- Initial_res
        } else {
          cat("File not found:", file_name, "\n")
        }
      }
      results <- do.call(rbind, results)
      RESULTS[3*(keri-1)+1,1:2]=t(apply(matrix(results[,2]%>%unlist,ncol=5,
                                               byrow=1)[,alpha_par_ind],
                                        1,function(x){exp(x)/sum(exp(x))})
      )%>%apply(2,mean)%>%round(2)
      
      j=2
      for(alg in c("Results/Results_additive_only/Results LBFGS/Results_hyper_NLOPT_LD_LBFGS_0.5",
                   "Results/Results_additive_only/Results Adam/Results_hyper_ADAM_0.1_100_0.5")){
        
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
        
        
        RESULTS[3*(keri-1)+j,1:2]= t(apply(results[,13:18][,alpha_par_ind],1, function(x){
          exp(as.numeric(unlist(x)))/sum(exp(as.numeric(unlist(x))))})
        )%>%apply(2,mean)%>%round(2)
        
        j=j+1
      }
      
      
      results <- vector("list", length(l))
      for (id in l) {
        file_name <- paste0(
          "Results/Results_additive_only/Results grid search/Initial_par_prot",keri,leakage,setup,"_",
          id, ".rda")
        if (file.exists(file_name)) {
          load(file_name)
          results[[id]] <-  Initial_res
        } else {
          cat("File not found:", file_name, "\n")
        }
      }
      results <- do.call(rbind, results)
      
      RESULTS[3*(keri-1)+1,3:4]=t(apply(matrix(results[,2]%>%unlist,ncol=5,
                                               byrow=1)[,alpha_par_ind],
                                        1,function(x){exp(x)/sum(exp(x))})
      )%>%apply(2,mean)%>%round(2)
      j=2
      for(alg in  c("Results/Results_additive_only/Results LBFGS/Results_hyper_NLOPT_LD_LBFGS_0.5",
                    "Results/Results_additive_only/Results Adam/Results_hyper_ADAM_0.1_100_0.5")){
        
        results <- vector("list", length(l))
        for (id in l) {
          file_name <- paste0(
            alg,"_prot",keri,leakage,setup,"_",
            id, ".rda")
          if (file.exists(file_name)) {
            load(file_name)
            results[[id]] <- res
          } else {
            cat("File not found:", file_name, "\n")
          }
        }
        results <- do.call(rbind, results)
        
        RESULTS[3*(keri-1)+j,3:4]=
          t(apply(results[,13:18][,alpha_par_ind],1, function(x){
            exp(as.numeric(unlist(x)))/sum(exp(as.numeric(unlist(x))))})
          )%>%apply(2,mean)%>%round(2)
        j=j+1
      }
      
    } 
    
    print(RESULTS)
  }
}





### HISTOGRAM#### 

for(keri in 1:4){
  par(mfrow=c(4,5))
  for(setup in 1:2){
    Alg=c("Results/Results_additive_only/Results LBFGS/Results_hyper_NLOPT_LD_LBFGS_0.5",
          "Results/Results_additive_only/Results Adam/Results_hyper_ADAM_0.1_100_0.5")
    
    for(leakage in c('no','yes')){
      if(keri%%2==0){
        alpha_par_ind=2:3
        
        title=c(expression(theta[Exp],alpha[1],alpha[2],alpha,"spectral k"))
        results <- vector("list", length(l))
        for (id in l) {
          file_name <- paste0(
            
            #"Results_hyper_yield1no1_"
            "Results/Results_additive_only/Results grid search/Initial_par_yield",keri,leakage,setup,"_",
            id, ".rda")
          if (file.exists(file_name)) {
            load(file_name)
            results[[id]] <- Initial_res
          } else {
            cat("File not found:", file_name, "\n")
          }
        }
        results <- do.call(rbind, results)
        results_GS=matrix(results[,2]%>%unlist,ncol=5,
                          byrow=1)
        
        results_GS[,alpha_par_ind]= t(apply(results_GS[,alpha_par_ind],1, function(x){
          exp(as.numeric(unlist(x)))/sum(exp(as.numeric(unlist(x))))})
        )
        
        alg=Alg[1]
        
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
        results_NLOPT=results[,13:18]
        results_NLOPT[,alpha_par_ind]= t(apply(results[,13:18][,alpha_par_ind],1, function(x){
          exp(as.numeric(unlist(x)))/sum(exp(as.numeric(unlist(x))))})
        )
        
        alg=Alg[2]
        
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
        results_Adam=results[,13:18]
        
        results_Adam[,alpha_par_ind]= t(apply(results[,13:18][,alpha_par_ind],1, function(x){
          exp(as.numeric(unlist(x)))/sum(exp(as.numeric(unlist(x))))})
        )
        
        col_green <- rgb(0, 1, 0, 0.5)  
        col_blue  <- rgb(0, 0, 1, 0.5)  
        
        
        for(i in 1:(length(title)-1)){
          if (length(unique(results_GS[,i])) <= 2) {
            if(i==1){
              if(setup==1&&leakage=="no"){
                par(mar=c(4,4,4,1))
              }else{
                par(mar=c(4,4,2,1))
              }
              hist(results_NLOPT[,i],col=col_blue,main=title[i],xlab='',
                   ylab=paste0("setup: ",setup,", leakage: ",leakage),cex.lab=1.5,
                   xlim=range(c(results_GS[,i],results_NLOPT[,i],results_Adam[,i])),breaks=10)
            }else{
              hist(results_NLOPT[,i],col=col_blue,main=title[i],xlab='',ylab='',
                   #ylab=paste0("setup: ",setup,", leakage: ",leakage),
                   xlim=range(c(results_GS[,i],results_NLOPT[,i],results_Adam[,i])),breaks=10)
              
            }
            abline(v=results_GS[1,i], col='grey', lwd=10)
            hist(results_Adam[,i],col=col_green,add=TRUE,breaks=10)
          }else{
            if(i==1){
              if(setup==1&&leakage=="no"){
                par(mar=c(4,4,4,1))
              }else{
                par(mar=c(4,4,2,1))
              }
              hist(results_GS[,i],main=title[i],xlab='',
                   ylab=paste0("setup: ",setup,", leakage: ",leakage),cex.lab=1.5,
                   xlim=range(c(results_GS[,i],results_NLOPT[,i],results_Adam[,i])),breaks=10)
            }else{
              hist(results_GS[,i],main=title[i],xlab='',ylab='',
                   #ylab=paste0("setup: ",setup,", leakage: ",leakage),
                   xlim=range(c(results_GS[,i],results_NLOPT[,i],results_Adam[,i])),breaks=10)
              
            }
            hist(results_NLOPT[,i],col=col_blue,add=TRUE,breaks=10)
            hist(results_Adam[,i],col=col_green,add=TRUE,breaks=10)
          }
          
        }
        
        i=5
        if (length(unique(results_GS[,i])) == 1) {
          hist(0,main=title[i],xlab='',ylab='',xlim=c(1,10),ylim=c(0,30),
               breaks=10)
          abline(v=results_GS[1,i], col='grey', lwd=10)
        }else{
          hist(results_GS[,i],main=title[i],xlab='',ylab=''
               ,breaks=10)
        }
        mtext(paste0("Kernel ", keri), outer = TRUE, cex = 1.5, line = -1)
        
        
        
      }else{
        alpha_par_ind=3:4
        title=c(expression(theta[Exp],theta[Ham],alpha[1],alpha[2],alpha))
        results <- vector("list", length(l))
        for (id in l) {
          file_name <- paste0(
            
            #"Results_hyper_yield1no1_"
            "Results/Results_additive_only/Results grid search/Initial_par_yield",keri,leakage,setup,"_",
            id, ".rda")
          if (file.exists(file_name)) {
            load(file_name)
            results[[id]] <- Initial_res
          } else {
            cat("File not found:", file_name, "\n")
          }
        }
        results <- do.call(rbind, results)
        results_GS=matrix(results[,2]%>%unlist,ncol=5,
                          byrow=1)
        
        results_GS[,alpha_par_ind]= t(apply(results_GS[,alpha_par_ind],1, function(x){
          exp(as.numeric(unlist(x)))/sum(exp(as.numeric(unlist(x))))})
        )
        
        alg=Alg[1]
        
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
        results_NLOPT=results[,13:18]
        results_NLOPT[,alpha_par_ind]= t(apply(results_NLOPT[,alpha_par_ind],1, function(x){
          exp(as.numeric(unlist(x)))/sum(exp(as.numeric(unlist(x))))})
        )
        
        alg=Alg[2]
        
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
        results_Adam=results[,13:18]
        
        results_Adam[,alpha_par_ind]= t(apply(results_Adam[,alpha_par_ind],1, function(x){
          exp(as.numeric(unlist(x)))/sum(exp(as.numeric(unlist(x))))})
        )
        
        col_green <- rgb(0, 1, 0, 0.5)  # green with 50% opacity
        col_blue  <- rgb(0, 0, 1, 0.5)  # blue with 50% opacity
        for(i in 1:length(title)){
          if (length(unique(results_GS[,i])) == 1) {
            if(i==1){
              if(setup==1&&leakage=="no"){
                par(mar=c(4,4,4,1))
              }else{
                par(mar=c(4,4,2,1))
              }
              hist(results_NLOPT[,i],col=col_blue,main=title[i],xlab='',
                   ylab=paste0("setup: ",setup,", leakage: ",leakage),cex.lab=1.5,
                   xlim=range(c(results_GS[,i],results_NLOPT[,i],results_Adam[,i])),breaks=10)
            }else{
              hist(results_NLOPT[,i],col=col_blue,main=title[i],xlab='',ylab='',
                   #ylab=paste0("setup: ",setup,", leakage: ",leakage),
                   xlim=range(c(results_GS[,i],results_NLOPT[,i],results_Adam[,i])),breaks=10)
              
            }
            abline(v=results_GS[1,i], col='grey', lwd=10)
            hist(results_Adam[,i],col=col_green,add=TRUE,breaks=10)
          }else{
            if(i==1){
              if(setup==1&&leakage=="no"){
                par(mar=c(4,4,4,1))
              }else{
                par(mar=c(4,4,2,1))
              }
              hist(results_GS[,i],main=title[i],xlab='',
                   ylab=paste0("setup: ",setup,", leakage: ",leakage),cex.lab=1.5,
                   xlim=range(c(results_GS[,i],results_NLOPT[,i],results_Adam[,i])),breaks=10)
            }else{
              hist(results_GS[,i],main=title[i],xlab='',ylab='',
                   #ylab=paste0("setup: ",setup,", leakage: ",leakage),
                   xlim=range(c(results_GS[,i],results_NLOPT[,i],results_Adam[,i])),breaks=10)
              
            }
            hist(results_NLOPT[,i],col=col_blue,add=TRUE,breaks=10)
            hist(results_Adam[,i],col=col_green,add=TRUE,breaks=10)
          }
          
        }
        mtext(paste0("Kernel ", keri), outer = TRUE, cex = 1.5, line = -1)
        
        
      }
      
      
    }
  }
}




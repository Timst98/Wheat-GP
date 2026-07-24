library(dplyr)
l=(1:30)


Algs=c("Results/Results_additive_only/Results grid search/Initial_par") 


cols=c('gray','blue','green','darkgreen')
Height=1024
tar='yield'

#What to plot
params=1
scores=ifelse(params,0,1)
### Boxplots #### 
for(tar in c('yield','prot')){
  for(params in c(0,1)){
    
    if(params){
      for( i in 1:3){
        for(j in 1:3){
          keri=3*(i-1)+j
          
          if (keri == 1) {
            kernel_geno <- 'dos'
            kernel_meteo <- 'rbf'
            param_names=c(expression(theta[GAU],theta[GBLUP],alpha[1],alpha[2],alpha)) 
            
          }
          if (keri == 2) {
            kernel_geno <- 'ham'
            kernel_meteo <- 'rbf'
            param_names=c(expression(theta[GAU],theta[HAM],alpha[1],alpha[2],alpha)) 
          }
          
          if (keri == 3) {
            kernel_geno <- 'spe'
            kernel_meteo <- 'rbf'
            param_names=c(expression(theta[Gau],alpha[1],alpha[2],alpha,k))
          }
          
          if (keri == 4) {
            kernel_geno <- 'dos'
            kernel_meteo <- 'exp'
            param_names=c(expression(theta[EXP],theta[GBLUP],alpha[1],alpha[2],alpha)) 
            
          }
          if (keri == 5) {
            kernel_geno <- 'ham'
            kernel_meteo <- 'exp'
            param_names=c(expression(theta[EXP],theta[HAM],alpha[1],alpha[2],alpha)) 
            
          }
          
          if (keri == 6) {
            kernel_geno <- 'spe'
            kernel_meteo <- 'exp'
            param_names=c(expression(theta[Exp],alpha[1],alpha[2],k))
          }
          if (keri == 7) {
            kernel_geno <- 'dos'
            GA_kernel=1
            kernel_meteo <- 'ali'
            param_names=c(expression(theta[GAK],theta[GBLUP],alpha[1],alpha[2],alpha)) 
            
          }
          if (keri == 8) {
            kernel_geno <- 'ham'
            GA_kernel=1
            kernel_meteo <- 'ali'
            param_names=c(expression(theta[GAK],theta[HAM],alpha[1],alpha[2],alpha)) 
            
          }
          
          if (keri == 9) {
            kernel_geno <- 'spe'
            GA_kernel=1
            kernel_meteo <- 'ali'
            param_names=c(expression(theta[GAK],alpha[1],alpha[2],alpha,k))
          }
          
          png(filename =
                paste0("Boxplots_grid_search/Parameters Sum model Kernels ",
                       c('1-3','4-6','7-9')[i],'_',tar,'.png'),
              height=Height,width=Height*3/4)
          par(mfrow=c(4,3))
          for(setup in 1:2){
            
            for(leakage in c('no','yes')){
              
              results <- vector("list", length(l))
              for (id in l) {
                file_name <- paste0(
                  Algs[1],"_",tar,keri,leakage,setup,"_",
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
              
              if(setup==1&&leakage=="no"){
                par(mar=c(2,4,4,.1))
              }else{
                par(mar=c(2,4,2,.1))
              }
              if(j==1){
                boxplot(
                  results_GS,
                  main = paste0('Ker Meteo: ', kernel_meteo, '  Ker Geno: ', kernel_geno),
                  ylab = paste0("setup: ", setup, ", leakage: ", leakage),
                  xaxt = "n"  # suppress x-axis so we can add math labels manually
                )
                
                # Add x-axis with math labels
                axis(1, at = 1:length(param_names), labels = param_names)
              }else{
                boxplot(
                  results_GS,
                  main = paste0('Kernel Meteo: ', kernel_meteo, '  Kernel Geno: ', kernel_geno),
                  ylab = '',
                  xaxt = "n"  # suppress x-axis so we can add math labels manually
                )
                
                # Add x-axis with math labels
                axis(1, at = 1:length(param_names), labels = param_names)
                
              }
             }
            
            
          }
          mtext(paste0("Kernels ", c('1-3','4-6','7-9')[i], ", ", tar,", Sum"), outer = TRUE, 
                cex = 1.5, line = -1.5)
        }
        dev.off()
      }
    }else{
      
      for( i in 1:3){
        for(j in 1:3){
          keri=3*(i-1)+j
          
          if (keri == 1) {
            kernel_geno <- 'dos'
            kernel_meteo <- 'rbf'
           
            
          }
          if (keri == 2) {
            kernel_geno <- 'ham'
            kernel_meteo <- 'rbf'
           
          }
          
          if (keri == 3) {
            kernel_geno <- 'spe'
            kernel_meteo <- 'rbf'
           
          }
          
          if (keri == 4) {
            kernel_geno <- 'dos'
            kernel_meteo <- 'exp'
                     
          }
          if (keri == 5) {
            kernel_geno <- 'ham'
            kernel_meteo <- 'exp'
                   }
          
          if (keri == 6) {
            kernel_geno <- 'spe'
            kernel_meteo <- 'exp'
         
          }
          if (keri == 7) {
            kernel_geno <- 'dos'
            GA_kernel=1
            kernel_meteo <- 'ali'
            
            
          }
          if (keri == 8) {
            kernel_geno <- 'ham'
            GA_kernel=1
            kernel_meteo <- 'ali'
            
          }
          
          if (keri == 9) {
            kernel_geno <- 'spe'
            GA_kernel=1
            kernel_meteo <- 'ali'
          }
          param_names=c('MSE','CRPS','LogS')
          png(filename =
                paste0("Boxplots_grid_search/Scores Sum model Kernels ",
                       c('1-3','4-6','7-9')[i],'_',tar,'.png'),
              height=Height,width=Height*3/4)
          par(mfrow=c(4,3))
          for(setup in 1:2){
            
            for(leakage in c('no','yes')){
              
              results <- vector("list", length(l))
              for (id in l) {
                file_name <- paste0(
                  Algs[1],"_",tar,keri,leakage,setup,"_",
                  id, ".rda")
                if (file.exists(file_name)) {
                  load(file_name)
                  results[[id]] <- Initial_res
                } else {
                  cat("File not found:", file_name, "\n")
                }
              }
              results <- do.call(rbind, results)[,4:6]%>%unlist
              results_GS=matrix(results,ncol=3)
              
              if(setup==1&&leakage=="no"){
                par(mar=c(2,4,4,.1))
              }else{
                par(mar=c(2,4,2,.1))
              }
              if(j==1){
                boxplot(
                  results_GS,
                  main = paste0('Ker Meteo: ', kernel_meteo, '  Ker Geno: ', kernel_geno),
                  ylab = paste0("setup: ", setup, ", leakage: ", leakage),
                  xaxt = "n"  # suppress x-axis so we can add math labels manually
                )
                
                # Add x-axis with math labels
                axis(1, at = 1:length(param_names), labels = param_names)
              }else{
                boxplot(
                  results_GS,
                  main = paste0('Kernel Meteo: ', kernel_meteo, '  Kernel Geno: ', kernel_geno),
                  ylab = '',
                  xaxt = "n"  # suppress x-axis so we can add math labels manually
                )
                
                # Add x-axis with math labels
                axis(1, at = 1:length(param_names), labels = param_names)
                
              }
            }
            
            
          }
          mtext(paste0("Kernels ", c('1-3','4-6','7-9')[i], ", ", tar,", Sum"), outer = TRUE, 
                cex = 1.5, line = -1.5)
        }
        dev.off()
      }
      
      
    }
  
  
  
}}



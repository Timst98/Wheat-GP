library(dplyr)

Adam_algs=c("Results/Results_additive_only/Results Adam/Results_hyper_ADAM_0.01_1000_0.5",
	    "Results/Results_additive_only/Results Adam/Results_hyper_ADAM_0.001_1000_0.5",
	    "Results/Results_additive_only/Results Adam/Results_hyper_ADAM_1e-04_1000_0.5",
            "Results/Results_additive_only/Results Adam/Results_hyper_ADAM_0.1_100_0.5")

Algs=c("Results/Results_additive_only/Results grid search/Initial_par",
       "Results/Results_additive_only/Results LBFGS/Results_hyper_NLOPT_LD_LBFGS_0.5",
       Adam_algs)[-c(4,5)]

cols=c('gray','blue',"#3CB371", 
             'green',"#ADFF2F",'darkgreen')[-c(4,5)]
Height=1024
insety=-3.67
#What to plot
params=1
scores=ifelse(params,0,1)
### Boxplots #### 
for(tar in c('yield','prot')){
  for(params in c(0,1)){

if(params){
  
  for(keri in (1:9)){
	  if (keri == 1) {
            kernel_geno <- 'dos'
            kernel_meteo <- 'rbf'
            param_names=c(expression(theta[GAU],theta[GBLUP],omega[1],alpha)) 
            
          }
          if (keri == 2) {
            kernel_geno <- 'ham'
            kernel_meteo <- 'rbf'
            param_names=c(expression(theta[GAU],theta[HAM],omega[1],alpha)) 
	    
          }
          
          if (keri == 3) {
            kernel_geno <- 'spe'
            kernel_meteo <- 'rbf'
            param_names=c(expression(theta[Gau],omega[1],alpha,k))
            
          }
          
          if (keri == 4) {
            kernel_geno <- 'dos'
            kernel_meteo <- 'exp'
            param_names=c(expression(theta[EXP],theta[GBLUP],omega[1],alpha)) 
            
          }
          if (keri == 5) {
            kernel_geno <- 'ham'
            kernel_meteo <- 'exp'
            param_names=c(expression(theta[EXP],theta[HAM],omega[1],alpha)) 
	    
            
          }
          
          if (keri == 6) {
            kernel_geno <- 'spe'
            kernel_meteo <- 'exp'
            param_names=c(expression(theta[Exp],omega[1],alpha,k))
	    
          }
          if (keri == 7) {
            kernel_geno <- 'dos'
            GA_kernel=1
            kernel_meteo <- 'ali'
            param_names=c(expression(theta[GAK],theta[GBLUP],omega[1],alpha)) 
            
          }
          if (keri == 8) {
            kernel_geno <- 'ham'
            GA_kernel=1
            kernel_meteo <- 'ali'
            param_names=c(expression(theta[GAK],theta[HAM],omega[1],alpha)) 
            
          }
          
          if (keri == 9) {
            kernel_geno <- 'spe'
            GA_kernel=1
            kernel_meteo <- 'ali'
            param_names=c(expression(theta[GAK],omega[1],alpha,k))
          	  
          }
    title=param_names
    png(filename = paste0("Boxplots/Kernel ", keri,"/", "Parameters Sum model ",tar,'.png'),
        height=Height,width=Height*length(title)/4)
    par(mfrow = c(4, length(title)),              
        oma = c(1, 1, 3, 1)         
    )          

    for(setup in 1:2){
      
      for(leakage in c('no','yes')){
        
          results <- vector("list", length(l))
          for (id in l) {
            file_name <- paste0(
              Algs[1],"_",tar,keri,ifelse(keri<7,'','GA'),leakage,setup,"_",
              id, ".rda")
            if (file.exists(file_name)) {
              load(file_name)
              results[[id]] <- Initial_res
            } else {
              cat("File not found:", file_name, "\n")
            }
          }
          results <- do.call(rbind, results)
          results_GS=matrix(results[,2]%>%unlist,ncol=4,
                            byrow=1)
          
          
          
          Res_matrices=list(results_GS)
          for(i in 2:length(Algs)){
            alg=Algs[i]
            
            results <- vector("list", length(l))
            for (id in l) {
              file_name <- paste0(
                alg,"_",tar,keri,ifelse(keri<7,'','GA'),leakage,setup,"_",
                id, ".rda")
              if (file.exists(file_name)) {
                load(file_name)
                results[[id]] <- res
              } else {
                cat("File not found:", file_name, "\n")
              }
            }
            results <- do.call(rbind, results)[,13:17]
            
            Res_matrices[[i]]=results
            
          }
          
          
          
          for(i in 1:length(title)){
            Param_matrix=c()
            for(j in 1:length(Algs)){
              Param_matrix=cbind(Param_matrix,Res_matrices[[j]][,i])
            }
            if(setup==1&&leakage=="no"){
              par(mar=c(2,4,4,.1))
            }else{
              par(mar=c(2,4,2,.1))
            }
            if(i==1){
              boxplot(Param_matrix,main=title[i],
                      ylab=paste0("setup: ",setup,", leakage: ",leakage),col=cols)
               
            }else{
              boxplot(Param_matrix,main=title[i],xlab='',col=cols)
            }
          }
        
        
      }
      
    }
   mtext(paste0("Kernel ", keri, ", ", tar,", Sum"), outer = TRUE, 
              cex = 1.5, line = 1)
        
   par(xpd = NA)  # allow drawing outside plot region
        
   legend("topright", inset = c(0, insety),
               legend = c("Grid Search", "LBFGS",
                          "Adam (lr,maxiter):",
                          "  0.01, 1000",
                          "  0.001, 1000",
                          "  0.0001, 1000",
                          "  0.1, 100")[-(5:6)],  # dropping 0.001 and 0.0001
               ncol = 7,
               fill = c(cols[1:2], NA, cols[3:6]),  # matching dropped items
               border = NA,
               cex = 1.2,
               bty = "n")
        
    dev.off()
  }
  
}else{
  title=c('Train N2ll','MSE','CRPS','LogS')
  for(keri in (1:9)){
   
    png(filename = paste0("Boxplots/Kernel ", keri,"/","Scores Sum model ", tar,'.png'),
        height=Height,width=Height*length(title)/4)
   par(mfrow = c(4, length(title)),              
        oma = c(1, 1, 3, 1)         
    )       
    for(setup in 1:2){
      for(leakage in c('no','yes')){
          results <- vector("list", length(l))
          for (id in l) {
            file_name <- paste0(
              Algs[1],"_",tar,keri,ifelse(keri<7,'','GA'),leakage,setup,"_",
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
          
          Res_matrices=list(results_GS)
          for(i in 2:length(Algs)){
            alg=Algs[i]
            
            results <- vector("list", length(l))
            for (id in l) {
              file_name <- paste0(
                alg,"_",tar,keri,ifelse(keri<7,'','GA'),leakage,setup,"_",
                id, ".rda")
              if (file.exists(file_name)) {
                load(file_name)
                results[[id]] <- res
              } else {
                cat("File not found:", file_name, "\n")
              }
            }
            results <- do.call(rbind, results)%>%unlist
            
            Res_matrices[[i]]=results[,c(9,1,2,10)]
	    if(i==2){
	     Res_matrices[[1]]=cbind(results[,12],Res_matrices[[1]])
	    }
            
          }
          
          for(i in 1:length(title)){
            Param_matrix=c()
            for(j in 1:length(Algs)){
              Param_matrix=cbind(Param_matrix,Res_matrices[[j]][,i])
            }
            if(setup==1&&leakage=="no"){
              par(mar=c(2,4,4,.1))
            }else{
              par(mar=c(2,4,2,.1))
            }
            if(i==1){
              boxplot(Param_matrix,main=title[i],
                      ylab=paste0("setup: ",setup,", leakage: ",leakage),col=cols)
	     
            }else{
              boxplot(Param_matrix,main=title[i],xlab='',col=cols)
            }
          }
          
         
          
          
      }
    }
 mtext(paste0("Kernel ", keri, ", ", tar,", Sum"), outer = TRUE, 
              cex = 1.5, line = 1)
        
   par(xpd = NA)  # allow drawing outside plot region
        
   legend("topright", inset = c(0, insety),
               legend = c("Grid Search", "LBFGS",
                          "Adam (lr,maxiter):",
                          "  0.01, 1000",
                          "  0.001, 1000",
                          "  0.0001, 1000",
                          "  0.1, 100")[-(5:6)],  # dropping 0.001 and 0.0001
               ncol = 7,
               fill = c(cols[1:2], NA, cols[3:6]),  # matching dropped items
               border = NA,
               cex = 1.2,
               bty = "n")
          
    dev.off()
  }
  
}


}}

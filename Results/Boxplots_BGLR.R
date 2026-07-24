library(dplyr)
l=(1:30)

cols=c('gray','blue','green','darkgreen')
Height=1024

### Boxplots #### 
keri=1
Algs=c('grid_search_','NLOPT_LD_LBFGS','ADAM_0.001_1000_','ADAM_0.1_100_')
for(tar in c('yield','prot')){
  title=c('MSE','CRPS','LogS')
  for(MODEL in 1:3){
   
    png(filename = paste0("Boxplots/BGLR/Scores full model ",MODEL,' ', tar,'.png'),
        height=Height,width=Height*length(title)/4)
    par(mfrow=c(4,3))
    for(setup in 1:2){
      for(leakage in c('no','yes')){
         Res_matrices=list()
         
         i=0
         for(alg in Algs){
         i=i+1
         results <- vector("list", length(l))
          for (id in l) {
            file_name <- paste0(
              'Results/Results BGLR/Results_hyper_',
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
              cat("File not found:", file_name, "\n")
            }
          }
          Res_matrices[[i]] <- (do.call(rbind, results)%>%as.matrix)[,c(1,2,10)]
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
          
          mtext(paste0("BGLR ", MODEL, ", ", tar,", Full"), outer = TRUE, 
                cex = 1.5, line = -1.5)
          
          
      }
    }
    dev.off()


    
  
}


}


#LMM 


for(tar in c('yield','prot')){
  title=c('MSE','CRPS')
  png(filename = paste0("Boxplots/LMM/Scores model ", tar,'.png'),
        height=Height,width=Height*length(title)/4)
  par(mfrow=c(4,2))
    for(setup in 1:2){
      for(leakage in c('no','yes')){

         Res_matrices=list()
         
         i=0
         for(MODEL in 1:3){
         i=i+1
         results <- vector("list", length(l))
        for (id in l) {
          file_name <- paste0(
            'Results/Results LMM/Results_LMM',
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
       
        
          Res_matrices[[i]] <- (do.call(rbind, results)%>%as.matrix)[,c(1,2)]
 }
        
        for(i in 1:length(title)){
          
            Param_matrix=c()
            for(j in 1:3){
              Param_matrix=cbind(Param_matrix,Res_matrices[[j]][,i])
            }
            if(setup==1&&leakage=="no"){
              par(mar=c(2,4,4,.1))
            }else{
              par(mar=c(2,4,2,.1))
            }
            if(i==1){
              boxplot(Param_matrix,main=title[i],
                      ylab=paste0("setup: ",setup,", leakage: ",leakage),col=2:4)
              legend('topright',c('Model 1','Model 2', 'Model 3'),col=2:4,fill=2:4)
            }else{
              boxplot(Param_matrix,main=title[i],xlab='',col=2:4)
            }
          }
          
          mtext(paste0("LMM ", MODEL, ", ", tar), outer = TRUE, 
                cex = 1.5, line = -1.5)

        
      }
    }
  
  dev.off()

}

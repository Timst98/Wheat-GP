#in case you want to run it locally
`%||%` <- function(a,b) if (!is.null(a) && !is.na(a) && a != "") a else b 

  
id=as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))%||%1
keri=as.integer(Sys.getenv("KERNEL"))%||%5     # which kernel combination
additive_only=as.integer(Sys.getenv("ADDONLY"))%||%0
product_only=as.integer(Sys.getenv("PRODONLY"))%||%0
long_Adam=as.integer(Sys.getenv("LONG_ADAM"))%||%0
GA_kernel=as.integer(Sys.getenv("GA_KERNEL"))%||%0
NLOPTR=as.integer(Sys.getenv("NLOPTR"))%||%0 #set to 0 for ADAM
grid_search_params=as.integer(Sys.getenv("GS_PARAMS"))%||%0
#setupp=ifelse(as.integer(Sys.getenv("SETUP"))==1,'1','2') # new variety or new environment
#we do it all parallel
#leakage=ifelse(as.integer(Sys.getenv("LEAKAGE"))==0,'no','yes') 

leak_spec='random'
noisy_obs = TRUE
noisy_sd= TRUE
train_prop = 0.8
nr_outer_split = 1

############################################################################################
# Optimization tuning parameters
batchsize=.5 # Batch the training data for 8 time speedup
NLOPT_alg='NLOPT_LD_LBFGS' 
if(long_Adam){
  lr=Lr=c(0.01)
  max_iter=1000
}else{
  lr=Lr=.1
  max_iter=100}
if(NLOPTR){Lr=1}

num_cores = length(Lr)*8
learning_rate_decay = 1

####################### Load kernels and GP training functions #############################################################
source('GP_modeling/Functions/Start.R')
source('GP_modeling/Functions/GP_functions.R')

################################################################################
cl=makeCluster(num_cores)  
registerDoParallel(cl)
# do parallel training for all 8 settings of (leakage,setup,target)
lr_tar=expand.grid(Lr,1:2,1:2,1:2)
Res <- foreach(
  ind = seq_len(nrow(lr_tar)),
  .combine  = rbind,
  .packages = Packages,
  .export = c(Export,'noisy_sd'),
  .errorhandling = "pass"
) %dopar% {
  
  tar     <- c("yield","prot")[lr_tar[ind,2]]
  leakage <- c("yes","no")[lr_tar[ind,3]]
  setupp  <- lr_tar[ind,4]
  lr      <- lr_tar[ind,1]
  
  Export <- character(0)   
  source("GP_modeling/Functions/GP_training.R", local=TRUE)
  res
}

stopCluster(cl)

##### Save results

result_directory=function(additive_only, product_only){
  if (additive_only) "Results_additive_only"
  else if (product_only) "Results_product_only"
  else "Results"
}
algorithm_directory=function(NLOPTR) if (NLOPTR) "Results LBFGS" else "Results Adam"

filename=function(res, tar, leakage, setupp, lr, GA_kernel, additive_only, product_only,
                  NLOPTR, NLOPT_alg, batchsize, max_iter, keri, id) {
  base=result_directory(additive_only, product_only)
  adir=algorithm_directory(NLOPTR)
  tag =if (GA_kernel != 0) "GA" else ""
  head=if (NLOPTR) paste0("Results_hyper_", NLOPT_alg, "_", batchsize, "_")
          else       paste0("Results_hyper_ADAM_", lr, "_", max_iter, "_", batchsize, "_")
  file.path(base, adir, paste0(head, tar, keri, tag, leakage, setupp, "_",ifelse(noisy_sd,'noisy_sd_',''),id, ".rda"))
}

for (ind in seq_len(nrow(lr_tar))) {
  tar    =c("yield","prot")[lr_tar[ind,2]]
  leakage=c("yes","no")[lr_tar[ind,3]]
  setupp =lr_tar[ind,4]
  lr     =lr_tar[ind,1]
  res    =as.numeric(Res[ind,])
  file_name=filename(res, tar, leakage, setupp, lr, GA_kernel,
                              additive_only, product_only, NLOPTR,
                              NLOPT_alg, batchsize, max_iter, keri, id)
  save(list="res", file=file_name)
  cat('\n Saved results under: ',file_name)
}

source('GP_modeling/Functions/GP_summary_output.R')

`%||%` <- function(a,b) if (!is.null(a) && !is.na(a) && a != "") a else b

id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID") %||% 1)

MODEL <- 3#"kerlmm"
grid_search_params <- 0
keri       <- 5
additive_only <- 0
product_only  <- 0
long_Adam  <- 1
GA_kernel  <- 0
NLOPTR     <- 0

full_data=0
                         
leak_spec <- "random"
noisy_obs <- TRUE
train_prop <- 0.8
nr_outer_split <- 1

# Optim tuning
batchsize <- 0.5

if (long_Adam) {
  lr <- Lr <- c(0.01)
  max_iter <- 1000
} else {
  lr <- Lr <- 0.1
  max_iter <- 100
}
if (NLOPTR) Lr <- 1



# Load common functions/data
source("GP_modeling/Functions/Start.R", local = TRUE)
source("GP_modeling/Functions/GP_functions.R", local = TRUE)

data <- merge(data, meteo, by = "Env")


for (setupp in 1) {
  for (leakage in c("yes", "no")[2]) {
    for(tar in c("yield", "prot")){
         Tab=data.frame(
              Component = c("G", "E", "GEI", "Residual"),
              Variance  = rep(0,4),
              Prop_total = rep(0,4),
              Prop_signal = rep(0,4)
            )
         count=0
          for(id in 1:30){
            file = paste0(
              'Results/Results BGLR/Summary_lmm_',ifelse(full_data,1,0.8),'_',
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
           if(file.exists(file)){
              load(file)
              Tab[,-1]=Tab[,-1]+tab[,-1]
            }else{
             cat('file does not exist: ',file,'\n')
             count=count+1}
            }
          cat('Target: ',tar,',Leakage: ',leakage, ',Setup: ', setupp, '\n')
          Tab[,-1]=Tab[,-1]/(30-count)
          print(Tab)
          }
      
      }
    
  }


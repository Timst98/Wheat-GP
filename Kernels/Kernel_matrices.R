
library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(dtw)
library(dbscan)
library(FNN)
library(caTools)
library(Metrics)
library(rsample)
library(Rfast)
library(lubridate)
library(patchwork)
library(tibble)
library(kernlab)
library(dtwclust)
library(Rcpp)

####################### KERNELS 

data <- read.table("Data/DATA.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
meteo  <- read.table("Data/METEO.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
geno   <- read.table("Data/GENO.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

data$yield <- (data$yield - mean(data$yield))/sd(data$yield)
data$prot <- (data$prot - mean(data$prot))/sd(data$prot)

# normalize the different weather variables
normalize_zscore <- function(x) {
  y <- (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE) 
  if(sd(x, na.rm=TRUE)==0){y <- 0}
  return(y)}
normalize_min_max <- function(x) {
  y <- (x - min(x, na.rm=TRUE)) / (max(x, na.rm=TRUE) - min(x, na.rm=TRUE))
  if(max(x, na.rm=TRUE)==0 & min(x, na.rm=TRUE)){y <- 0}
  return(y)}

c2 <- ncol(meteo)
meteo[,-c2] <- lapply(meteo[,-c2], normalize_zscore)
varNames <- sub("_.*", "", colnames(meteo)[3:ncol(meteo)]) %>% unique()
meteo1 <- meteo %>%
  dplyr::select(matches(paste(varNames[1:7], collapse = "|")))
meteo <- cbind(meteo$Env, meteo1)
colnames(meteo)[1] <- 'Env'

#################################################################################
# HAMMING for GENO

# count the number of equal or NA entries
hamming_distance_matrix <- function(maat) {
  n <- nrow(maat)
  hamming_matrix <- outer(1:n, 1:n, Vectorize(function(i, j) {
    v <- (maat[i, ] == maat[j, ]) | is.na(maat[i, ]) | is.na(maat[j, ]) 
    return((ncol(maat) - sum(v)) / ncol(maat)) 
  }))
  return(hamming_matrix)
}
 geno_dist_hamming <- hamming_distance_matrix(geno[,-1])
 rownames(geno_dist_hamming) <- geno$variety_name
 colnames(geno_dist_hamming) <- geno$variety_name
 write.table(geno_dist_hamming, file = "Kernels/GENO_HAMMING.txt", sep = "\t", row.names = FALSE, quote = FALSE)

geno_dist_hamming   <- read.table("Kernels/GENO_HAMMING.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rownames(geno_dist_hamming) <- colnames(geno_dist_hamming)

Dmat <- geno_dist_hamming[data$variety_name, data$variety_name]/max(geno_dist_hamming)
Dmat <- Dmat %>% round(digits=2)
valu <- Dmat %>% as.matrix() %>% c() %>% unique() %>% sort()
vario <- matrix(NA, length(valu), 3)

for(ii in 1:length(valu)){
  ind <- which(Dmat == valu[ii], arr.ind=TRUE)
  varioo <- matrix(NA, nrow(ind), 2)
  dista <- valu[ii]
  for(jj in 1:nrow(ind)){
    varioo[jj,1] <- 0.5*(data$yield[ind[jj,1]] - data$yield[ind[jj,2]])^2
    varioo[jj,2] <- 0.5*(data$prot[ind[jj,1]] - data$prot[ind[jj,2]])^2   }
  vario[ii,1] <- dista
  vario[ii,2] <- mean(varioo[,1], na.rm=TRUE)
  vario[ii,3] <- mean(varioo[,2], na.rm=TRUE)   }

rangeY <- 0.4
sillY <- 1
rangeP <- 0.4
sillP <- 1

expY <- (sillY - 0)* (exp(0) - exp(-valu/(rangeY/3))) + 0
expP <- (sillP - 0)* (exp(0) - exp(-valu/(rangeP/3))) + 0



#################################################################################
## Exponential for Meteo
# EUCLIDEAN for METEO

meteo_dist_eucl <- as.matrix(dist(meteo[,-1], method = "euclidean"))
rownames(meteo_dist_eucl) <- meteo$Env
colnames(meteo_dist_eucl) <- meteo$Env

Dmat <- meteo_dist_eucl[data$Env, data$Env]/max(meteo_dist_eucl)
Dmat <- Dmat %>% round(digits=2)
valu <- Dmat %>% as.matrix() %>% c() %>% unique() %>% sort()
vario <- matrix(NA, length(valu), 3)

for(ii in 1:length(valu)){
  ind <- which(Dmat == valu[ii], arr.ind=TRUE)
  varioo <- matrix(NA, nrow(ind), 2)
  dista <- valu[ii]
  for(jj in 1:nrow(ind)){
    varioo[jj,1] <- 0.5*(data$yield[ind[jj,1]] - data$yield[ind[jj,2]])^2
    varioo[jj,2] <- 0.5*(data$prot[ind[jj,1]] - data$prot[ind[jj,2]])^2   }
  vario[ii,1] <- dista
  vario[ii,2] <- mean(varioo[,1], na.rm=TRUE)
  vario[ii,3] <- mean(varioo[,2], na.rm=TRUE)   }

rangeY <- 0.4
sillY <- 1
rangeP <- 0.4
sillP <- 1

expY <- (sillY - 0)* (exp(0) - exp(-valu/(rangeY/3))) + 0
expP <- (sillP - 0)* (exp(0) - exp(-valu/(rangeP/3))) + 0


###################################################################################
## sum of 2
Dmat <- meteo_dist_eucl[data$Env, data$Env]/max(meteo_dist_eucl) + 
  geno_dist_hamming[data$variety_name, data$variety_name]/max(geno_dist_hamming)
Dmat <- Dmat %>% round(digits=2)
valu <- Dmat %>% as.matrix() %>% c() %>% unique() %>% sort()
vario <- matrix(NA, length(valu), 3)

for(ii in 1:length(valu)){
  ind <- which(Dmat == valu[ii], arr.ind=TRUE)
  varioo <- matrix(NA, nrow(ind), 2)
  dista <- valu[ii]
  for(jj in 1:nrow(ind)){
    varioo[jj,1] <- 0.5*(data$yield[ind[jj,1]] - data$yield[ind[jj,2]])^2
    varioo[jj,2] <- 0.5*(data$prot[ind[jj,1]] - data$prot[ind[jj,2]])^2   }
  vario[ii,1] <- dista
  vario[ii,2] <- mean(varioo[,1], na.rm=TRUE)
  vario[ii,3] <- mean(varioo[,2], na.rm=TRUE)   }

rangeY <- 0.8
sillY <- 1
rangeP <- 0.9
sillP <- 1

expY <- (sillY - 0)* (exp(0) - exp(-valu/(rangeY/3))) + 0
expP <- (sillP - 0)* (exp(0) - exp(-valu/(rangeP/3))) + 0



###################################################################################
# STRING KERNELS for GENO
geno_list <- vector("list", nrow(geno))
for(g in 1:nrow(geno)){
  geno_list[[g]] <- as.character(geno[g,-1]) 
}

string_kernel <- stringdot(type = "spectrum", length = 2)
geno_spec <- kernelMatrix(string_kernel, geno_list) 
rownames(geno_spec) <- geno$variety_name
colnames(geno_spec) <- geno$variety_name
write.table(geno_spec, file = "Kernels/GENO_SPEC2.txt", sep = "\t", row.names = FALSE, quote = FALSE)

string_kernel <- stringdot(type = "spectrum", length = 3)
geno_spec <- kernelMatrix(string_kernel, geno_list) 
rownames(geno_spec) <- geno$variety_name
colnames(geno_spec) <- geno$variety_name
write.table(geno_spec, file = "Kernels/GENO_SPEC3.txt", sep = "\t", row.names = FALSE, quote = FALSE)

string_kernel <- stringdot(type = "spectrum", length = 4)
geno_spec <- kernelMatrix(string_kernel, geno_list) 
rownames(geno_spec) <- geno$variety_name
colnames(geno_spec) <- geno$variety_name
write.table(geno_spec, file = "Kernels/GENO_SPEC4.txt", sep = "\t", row.names = FALSE, quote = FALSE)

string_kernel <- stringdot(type = "spectrum", length = 5)
geno_spec <- kernelMatrix(string_kernel, geno_list) 
rownames(geno_spec) <- geno$variety_name
colnames(geno_spec) <- geno$variety_name
write.table(geno_spec, file = "Kernels/GENO_SPEC5.txt", sep = "\t", row.names = FALSE, quote = FALSE)

string_kernel <- stringdot(type = "spectrum", length = 6)
geno_spec <- kernelMatrix(string_kernel, geno_list) 
rownames(geno_spec) <- geno$variety_name
colnames(geno_spec) <- geno$variety_name
write.table(geno_spec, file = "Kernels/GENO_SPEC6.txt", sep = "\t", row.names = FALSE, quote = FALSE)

string_kernel <- stringdot(type = "spectrum", length = 7)
geno_spec <- kernelMatrix(string_kernel, geno_list) 
rownames(geno_spec) <- geno$variety_name
colnames(geno_spec) <- geno$variety_name
write.table(geno_spec, file = "Kernels/GENO_SPEC7.txt", sep = "\t", row.names = FALSE, quote = FALSE)

string_kernel <- stringdot(type = "spectrum", length = 8)
geno_spec <- kernelMatrix(string_kernel, geno_list) 
rownames(geno_spec) <- geno$variety_name
colnames(geno_spec) <- geno$variety_name
write.table(geno_spec, file = "Kernels/GENO_SPEC8.txt", sep = "\t", row.names = FALSE, quote = FALSE)

string_kernel <- stringdot(type = "spectrum", length = 10)
geno_spec <- kernelMatrix(string_kernel, geno_list) 
rownames(geno_spec) <- geno$variety_name
colnames(geno_spec) <- geno$variety_name
write.table(geno_spec, file = "Kernels/GENO_SPEC10.txt", sep = "\t", row.names = FALSE, quote = FALSE)

string_kernel <- stringdot(type = "spectrum", length = 9)
geno_spec <- kernelMatrix(string_kernel, geno_list) 
rownames(geno_spec) <- geno$variety_name
colnames(geno_spec) <- geno$variety_name
write.table(geno_spec, file = "Kernels/GENO_SPEC9.txt", sep = "\t", row.names = FALSE, quote = FALSE)

string_kernel <- stringdot(type = "spectrum", length = 1)
geno_spec <- kernelMatrix(string_kernel, geno_list) 
rownames(geno_spec) <- geno$variety_name
colnames(geno_spec) <- geno$variety_name
write.table(geno_spec, file = "Kernels/GENO_SPEC1.txt", sep = "\t", row.names = FALSE, quote = FALSE)



#################################################################################
# Alignment for meteo  (from Marco Cuturi's webpage)
# Useful constants
LOG0 <- -10000  # log(0)

# LOGP function: stable computation of log(exp(x) + exp(y))
LOGP <- function(x, y) {
  if (x > y) {
    return(x + log1p(exp(y - x)))
  } else {
    return(y + log1p(exp(x - y)))
  }
}

# Global Alignment Kernel function
logGAK <- function(seq1, seq2, sigma, triangular = 0) {
  nX <- nrow(seq1)  # Number of rows in seq1
  nY <- nrow(seq2)  # Number of rows in seq2
  dimvect <- ncol(seq1)  # Dimension of the time series vectors
  
  # Length of a column for dynamic programming
  cl <- nY + 1
  
  # Initialize logM to store two successive columns (dynamic programming table)
  logM <- matrix(LOG0, nrow = cl, ncol = 2)
  logM[1, 1] <- 0  # Initialize the lower-left cell with log(1) = 0
  
  # Maximum of abs(i - j) when 1 <= i <= nX and 1 <= j <= nY
  trimax <- max(nX - 1, nY - 1)
  
  # Triangular coefficients initialization
  logTriangularCoefficients <- rep(0, trimax + 1)
  if (triangular > 0) {
    for (i in seq_len(min(trimax + 1, triangular)) - 1) {
      logTriangularCoefficients[i + 1] <- log(1 - i / triangular)
    }
  }
  
  # Sigma factor for Gaussian kernel
  Sig <- -1 / (2 * sigma^2)
  
  # Dynamic programming to compute the log of the Global Alignment Kernel
  cur <- 2
  old <- 1
  
  for (i in 1:nX) {
    logM[, cur] <- LOG0  # Reset the current column
    for (j in 1:nY) {
      if (logTriangularCoefficients[abs(i - j) + 1] > LOG0) {
        # Compute the Gaussian kernel value
        diff_sq <- sum((seq1[i, ] - seq2[j, ])^2)
        gram <- logTriangularCoefficients[abs(i - j) + 1] + diff_sq * Sig
        gram <- gram - log(2 - exp(gram))
        
        # Update logM for dynamic programming
        frompos1 <- logM[j + 1, old]
        frompos2 <- logM[j, cur]
        frompos3 <- logM[j, old]
        aux <- LOGP(frompos1, frompos2)
        logM[j + 1, cur] <- LOGP(aux, frompos3) + gram
      }
    }
    # Swap cur and old
    cur <- 3 - cur
    old <- 3 - old
  }
  
  # Return the final result
  return(logM[nY + 1, old])
}

varNames <- sub("_.*", "", colnames(meteo)[3:ncol(meteo)]) %>% unique()

for(s in seq(2,8,length.out=10)){

Kernel_GA <- matrix(NA, nrow=nrow(meteo), ncol=nrow(meteo))
sigma <- matrix(NA, nrow=nrow(meteo), ncol=nrow(meteo))
for(j1 in 1:nrow(meteo)){
  for(j2 in 1:j1){
    
    meteoM1 <- matrix(NA, nrow = 6, ncol = 7)
    meteoM2 <- matrix(NA, nrow = 6, ncol = 7)
    for(jj in 1:7){
      mm1 <- as.numeric(meteo[j1,grep(varNames[jj], colnames(meteo))])
      mm2 <- as.numeric(meteo[j2,grep(varNames[jj], colnames(meteo))])
      meteoM1[1:length(mm1),jj] <- mm1
      meteoM2[1:length(mm2),jj] <- mm2
    }
    
    t1 <- logGAK(meteoM1, meteoM2, sigma = s, triangular = 0)
    t2 <- logGAK(meteoM2, meteoM2, sigma = s, triangular = 0)
    t3 <- logGAK(meteoM1, meteoM1, sigma = s, triangular= 0)
    
    Kernel_GA[j1,j2] <- exp(t1-.5*(t2+t3))
    Kernel_GA[j2,j1] <- exp(t1-.5*(t2+t3))
    
  }}

rownames(Kernel_GA) <- meteo$Env
colnames(Kernel_GA) <- meteo$Env

write.table(Kernel_GA, file = paste0("Kernels/METEO_GA_sigma", s, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
}

# seq(2,8,length.out=10)






library(ggplot2)
library(dplyr)
library(tidyr)

setwd("C:/Users/leafr/OneDrive - Universitaet Bern/Desktop/Projects/Agroscope/1_DEF_Tim/Results")

prot <- matrix(NA, 4, 8)
yield <- matrix(NA, 4, 8)
for(kk in 1:4){
  
  pp <- read.table(paste0("Results_hyper_prot", kk, "yes.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  prot[kk,] <-  pp %>% colMeans()
  yy <- read.table(paste0("Results_hyper_yield", kk, "yes.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  yield[kk,] <-  yy %>% colMeans()
}

print(round(yield[,1:2],digits=2))
print(round(yield[,2],digits=2))
print(yield)

print(round(prot[,1:2],digits=2))
print(round(prot[,2],digits=2))
print(prot)






pp <- read.table(paste0("Results_hyper_yield1no2.txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
pp %>% colMeans()















prot$row_id <- paste0('ker', seq(1,4,1))
prot <- pivot_longer(prot, cols = -row_id, names_to = "method", values_to = "value")
prot <- prot %>% filter(row_id == 'ker1' | method == 'kernel') %>% arrange(value)

yield <- as.data.frame(yield[,-2])
rownames(yield) <- paste0('ker', seq(1,4,1))
colnames(yield) <- c('kernel', 'global', 'variety', 'local')
yield$row_id <- paste0('ker', seq(1,4,1))
yield <- pivot_longer(yield, cols = -row_id, names_to = "method", values_to = "value")
yield <- yield %>% filter(row_id == 'ker1' | method == 'kernel') %>% arrange(value)




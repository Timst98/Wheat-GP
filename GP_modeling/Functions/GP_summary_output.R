cat("\n \n \n ==================== GP RUN ====================\n")

kernel_mode <- if (additive_only) {
  "additive_only"
} else if (product_only) {
  "product_only"
} else {
  "full"
}

cat("Kernel setup:",
    kernel_mode,
    "| Kernel:", keri,
    "| optimizer:",
    if (NLOPTR) {
      paste0("LBFGS (", NLOPT_alg, ")")
    } else {
      paste0("Adam (lr=", lr, ", iters=", max_iter, ")")
    },
    "| batchsize:", batchsize,
    "\n------------------------------------------------\n", sep=" ")

fmt_vec <- function(x, digits=6) {
  x <- as.numeric(x)
  x[is.na(x)] <- NA_real_
  paste0("[", paste(ifelse(is.na(x), "NA",
                           formatC(x, digits=digits, format="fg")),
                    collapse=", "), "]")
}

# column blocks as in GP_training.R
opt_cols  <- 13:17
init_cols <- 19:23

for (ind in seq_len(nrow(lr_tar))) {
  tar     <- c("yield","prot")[lr_tar[ind,2]]
  leakage <- c("yes","no")[lr_tar[ind,3]]
  setupp  <- lr_tar[ind,4]
  lr      <- lr_tar[ind,1]
  
  r <- as.numeric(Res[ind, ])
  
  mse  <- r[1]
  crps <- r[2]
  logs <- r[10]
  t_total=r[26]
  
  init_par <- r[init_cols]
  opt_par  <- r[opt_cols]
  
  init_par <- init_par[!is.na(init_par)]
  opt_par  <- opt_par[!is.na(opt_par)]
  cat(sprintf("\n[target=%s | leakage=%s | setup=%s | lr=%s]\n",
              tar, leakage, setupp, lr))
  cat(sprintf("  MSE  = %.6f\n", mse))
  cat(sprintf("  CRPS = %.6f\n", crps))
  cat(sprintf("  LogS = %.6f\n", logs))
  cat(sprintf("  Total time (optimization+inference) = %.6f\n", t_total))
  cat("Initial parameters chosen by grid search:", fmt_vec(init_par), "\n")
  cat("Optimized parameters :", fmt_vec(opt_par),  "\n")
}


cat("\n================================================\n")

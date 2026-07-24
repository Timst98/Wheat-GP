suppressPackageStartupMessages({
  library(dplyr)
})

dir.create("Results", showWarnings = FALSE, recursive = TRUE)

IDS <- 1:30
AGG <- function(x) median(as.numeric(x), na.rm = TRUE)

Adam_full <- "Results/Results Adam/Results_hyper_ADAM_0.01_1000_0.5"
Adam_add  <- "Results/Results_additive_only/Results Adam/Results_hyper_ADAM_0.01_1000_0.5"
Adam_prod <- "Results/Results_product_only/Results Adam/Results_hyper_ADAM_0.01_1000_0.5"
Adam_env  <- "Results/Results_env_only/Results Adam/Results_hyper_ADAM_0.01_1000_0.5"
Adam_gen  <- "Results/Results_gen_only/Results Adam/Results_hyper_ADAM_0.01_1000_0.5"

alg <- Adam_full

DIR_BGLR <- "Results/Results BGLR"
DIR_LMM  <- "Results/Results LMM"

load_env <- function(path) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  e
}

as_numvec <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.matrix(x) || is.data.frame(x)) return(as.numeric(x[1, ]))
  as.numeric(x)
}

# robustly extract the results vector from an .rda env
extract_results_vec <- function(env, tar = c("yield","prot")) {
  tar <- match.arg(tar)

  if (exists("Res", envir = env, inherits = FALSE)) {
    R <- get("Res", envir = env)
    if (is.list(R) && length(R) >= 2) {
      block <- if (tar == "yield") 1 else 2
      return(as_numvec(R[[block]]))
    }
  }
  if (exists("res", envir = env, inherits = FALSE)) {
    return(as_numvec(get("res", envir = env)))
  }
  nms <- ls(envir = env)
  if (length(nms) == 0) return(NULL)
  as_numvec(get(nms[1], envir = env))
}

summarise_from_vecs <- function(vecs) {
  vecs <- vecs[!vapply(vecs, is.null, logical(1))]
  if (length(vecs) == 0) return(c(MSE=NA_real_, CRPS=NA_real_, LogS=NA_real_))
  mlen <- max(vapply(vecs, length, integer(1)))
  M <- do.call(rbind, lapply(vecs, function(v) { length(v) <- mlen; v }))
  out <- c(MSE=NA_real_, CRPS=NA_real_, LogS=NA_real_)
  if (ncol(M) >= 2) { out["MSE"] <- AGG(M[,1]); out["CRPS"] <- AGG(M[,2]) }
  if (ncol(M) >= 10) out["LogS"] <- AGG(M[,10])
  out
}

col_idx <- function(tar, metric, leakage) {
  base <- if (tar == "yield") 0 else 6
  off_metric <- switch(metric, MSE=0, CRPS=2, LogS=4)
  off_leak <- if (leakage == "no") 1 else 2
  base + off_metric + off_leak
}
fill_metrics <- function(TAB, row, tar, leakage, vals) {
  TAB[row, col_idx(tar,"MSE",  leakage)]  <- vals["MSE"]
  TAB[row, col_idx(tar,"CRPS", leakage)]  <- vals["CRPS"]
  TAB[row, col_idx(tar,"LogS", leakage)]  <- vals["LogS"]
  TAB
}

kernel_path <- function(prefix, tar, keri, leakage, setup, id) {
  paste0(prefix, "_", tar, keri, leakage, setup, "_", id, ".rda")
}
read_kernel_vec <- function(prefix, tar, keri, leakage, setup, id) {
  fp <- kernel_path(prefix, tar, keri, leakage, setup, id)
  if (!file.exists(fp)) return(NULL)
  e <- load_env(fp)
  if (exists("res", envir = e, inherits = FALSE)) return(as_numvec(get("res", envir = e)))
  extract_results_vec(e, tar)
}

compute_averages <- function(tar, leakage, setup) {
  R <- NULL
  for (keri in 1:6) {
    mats <- lapply(IDS, function(id) {
      fp <- paste0(alg, "_", tar, keri, leakage, setup, "_", id, ".rda")
      if (!file.exists(fp)) return(NULL)
      e <- load_env(fp)
      v <- if (exists("res", envir = e, inherits = FALSE)) as_numvec(get("res", envir = e)) else extract_results_vec(e, tar)
      if (is.null(v) || length(v) < 8) return(NULL)
      v[1:8]
    })
    mats <- mats[!vapply(mats, is.null, logical(1))]
    if (length(mats) == 0) next
    R <- rbind(R, do.call(rbind, mats))
  }
  if (is.null(R)) {
    return(list(Global=c(MSE=NA_real_,CRPS=NA_real_),
                Variety=c(MSE=NA_real_,CRPS=NA_real_),
                Env=c(MSE=NA_real_,CRPS=NA_real_)))
  }
  vv <- apply(R, 2, AGG)
  list(
    Global  = c(MSE=vv[3], CRPS=vv[4]),
    Variety = c(MSE=vv[5], CRPS=vv[6]),
    Env     = c(MSE=vv[7], CRPS=vv[8])
  )
}

compute_lmm <- function(MODEL, tar, leakage, setup) {
  vecs <- lapply(IDS, function(id) {
    fp <- paste0(DIR_LMM, "/Results_LMM", MODEL, "_", tar, 1, leakage, setup, "_", id, ".rda")
    if (!file.exists(fp)) return(NULL)
    e <- load_env(fp)
    extract_results_vec(e, tar)
  })
  summarise_from_vecs(vecs)
}

bglr1_path <- function(tar, keri, leakage, setup, id) {
  file.path(DIR_BGLR, paste0("Results_hyper_ADAM_0.01_1000_3_", tar, keri, leakage, setup, "_", id, ".rda"))
}
kerlmm_path <- function(tar, keri, leakage, setup, id) {
  file.path(DIR_BGLR, paste0("Results_hyper_ADAM_0.01_1000_kerlmm_", tar, keri, leakage, setup, "_", id, ".rda"))
}

compute_bglr1 <- function(tar, keri, leakage, setup) {
  vecs <- lapply(IDS, function(id) {
    fp <- bglr1_path(tar, keri, leakage, setup, id)
    if (!file.exists(fp)) return(NULL)
    e <- load_env(fp)
    extract_results_vec(e, tar)
  })
  summarise_from_vecs(vecs)
}

compute_kerlmm <- function(tar, keri, leakage, setup) {
  vecs <- lapply(IDS, function(id) {
    fp <- kerlmm_path(tar, keri, leakage, setup, id)
    if (!file.exists(fp)) return(NULL)
    e <- load_env(fp)
    extract_results_vec(e, tar)
  })
  summarise_from_vecs(vecs)
}

naify_kerlmm_logs <- function(TAB) {
  rn <- "Kernel LMM"
  if (!rn %in% rownames(TAB)) return(TAB)

  cols <- c(
    "Yield LogS | no leak","Yield LogS | leak",
    "Prot LogS | no leak","Prot LogS | leak"
  )
  cols <- cols[cols %in% colnames(TAB)]
  if (length(cols) == 0) return(TAB)

  z <- TAB[rn, cols]
  TAB[rn, cols] <- ifelse(is.na(z), NA_real_, ifelse(z == 0, NA_real_, z))
  TAB
}

build_table_for_setup <- function(setup) {

  rows <- c(
    "K1 full","K2 full","K3 full","K4 full","K5 full",
    "K5 additive","K5 product","K5 env-only","K5 gen-only",
    "LMM model 1","BGLR model 3","Kernel LMM",
    "Global average","Environmental average","Variety average"
  )

  cols <- c(
    "Yield MSE | no leak","Yield MSE | leak",
    "Yield CRPS | no leak","Yield CRPS | leak",
    "Yield LogS | no leak","Yield LogS | leak",
    "Prot MSE | no leak","Prot MSE | leak",
    "Prot CRPS | no leak","Prot CRPS | leak",
    "Prot LogS | no leak","Prot LogS | leak"
  )

  TAB <- matrix(NA_real_, nrow=length(rows), ncol=length(cols),
                dimnames=list(rows, cols))

  for (leakage in c("no","yes")) {

    for (k in 1:5) {
      y_vecs <- lapply(IDS, function(id) read_kernel_vec(Adam_full, "yield", k, leakage, setup, id))
      p_vecs <- lapply(IDS, function(id) read_kernel_vec(Adam_full, "prot",  k, leakage, setup, id))
      TAB <- fill_metrics(TAB, paste0("K",k," full"), "yield", leakage, summarise_from_vecs(y_vecs))
      TAB <- fill_metrics(TAB, paste0("K",k," full"), "prot",  leakage, summarise_from_vecs(p_vecs))
    }

    y_vecs <- lapply(IDS, function(id) read_kernel_vec(Adam_add,  "yield", 5, leakage, setup, id))
    p_vecs <- lapply(IDS, function(id) read_kernel_vec(Adam_add,  "prot",  5, leakage, setup, id))
    TAB <- fill_metrics(TAB, "K5 additive", "yield", leakage, summarise_from_vecs(y_vecs))
    TAB <- fill_metrics(TAB, "K5 additive", "prot",  leakage, summarise_from_vecs(p_vecs))

    y_vecs <- lapply(IDS, function(id) read_kernel_vec(Adam_prod, "yield", 5, leakage, setup, id))
    p_vecs <- lapply(IDS, function(id) read_kernel_vec(Adam_prod, "prot",  5, leakage, setup, id))
    TAB <- fill_metrics(TAB, "K5 product", "yield", leakage, summarise_from_vecs(y_vecs))
    TAB <- fill_metrics(TAB, "K5 product", "prot",  leakage, summarise_from_vecs(p_vecs))

    k_envgen <- 1
    y_vecs <- lapply(IDS, function(id) read_kernel_vec(Adam_env, "yield", k_envgen, leakage, setup, id))
    p_vecs <- lapply(IDS, function(id) read_kernel_vec(Adam_env, "prot",  k_envgen, leakage, setup, id))
    TAB <- fill_metrics(TAB, "K5 env-only", "yield", leakage, summarise_from_vecs(y_vecs))
    TAB <- fill_metrics(TAB, "K5 env-only", "prot",  leakage, summarise_from_vecs(p_vecs))

    y_vecs <- lapply(IDS, function(id) read_kernel_vec(Adam_gen, "yield", k_envgen, leakage, setup, id))
    p_vecs <- lapply(IDS, function(id) read_kernel_vec(Adam_gen, "prot",  k_envgen, leakage, setup, id))
    TAB <- fill_metrics(TAB, "K5 gen-only", "yield", leakage, summarise_from_vecs(y_vecs))
    TAB <- fill_metrics(TAB, "K5 gen-only", "prot",  leakage, summarise_from_vecs(p_vecs))

    TAB <- fill_metrics(TAB, "LMM model 1", "yield", leakage, compute_lmm(1, "yield", leakage, setup))
    TAB <- fill_metrics(TAB, "LMM model 1", "prot",  leakage, compute_lmm(1, "prot",  leakage, setup))

    TAB <- fill_metrics(TAB, "BGLR model 3", "yield", leakage, compute_bglr1("yield", 5, leakage, setup))
    TAB <- fill_metrics(TAB, "BGLR model 3", "prot",  leakage, compute_bglr1("prot",  5, leakage, setup))

    TAB <- fill_metrics(TAB, "Kernel LMM", "yield", leakage, compute_kerlmm("yield", 5, leakage, setup))
    TAB <- fill_metrics(TAB, "Kernel LMM", "prot",  leakage, compute_kerlmm("prot",  5, leakage, setup))

    avY <- compute_averages("yield", leakage, setup)
    avP <- compute_averages("prot",  leakage, setup)

    TAB["Global average", col_idx("yield","MSE",  leakage)] <- avY$Global["MSE"]
    TAB["Global average", col_idx("yield","CRPS", leakage)] <- avY$Global["CRPS"]
    TAB["Global average", col_idx("prot","MSE",   leakage)] <- avP$Global["MSE"]
    TAB["Global average", col_idx("prot","CRPS",  leakage)] <- avP$Global["CRPS"]

    if (leakage == "yes" || setup == 1) {
      TAB["Variety average", col_idx("yield","MSE",  leakage)] <- avY$Variety["MSE"]
      TAB["Variety average", col_idx("yield","CRPS", leakage)] <- avY$Variety["CRPS"]
      TAB["Variety average", col_idx("prot","MSE",   leakage)] <- avP$Variety["MSE"]
      TAB["Variety average", col_idx("prot","CRPS",  leakage)] <- avP$Variety["CRPS"]
    }

    if (leakage == "yes" || setup == 2) {
      TAB["Environmental average", col_idx("yield","MSE",  leakage)] <- avY$Env["MSE"]
      TAB["Environmental average", col_idx("yield","CRPS", leakage)] <- avY$Env["CRPS"]
      TAB["Environmental average", col_idx("prot","MSE",   leakage)] <- avP$Env["MSE"]
      TAB["Environmental average", col_idx("prot","CRPS",  leakage)] <- avP$Env["CRPS"]
    }
  }

  TAB <- naify_kerlmm_logs(TAB)
  TAB <- round(TAB, 3)

  cat("\n========================================\n")
  cat("Results summary (setup =", setup, "), medians over ids 1..30\n")
  cat("========================================\n\n")
  print(as.data.frame(TAB, check.names = FALSE))

  invisible(TAB)
}

# ---------- labels ----------
method_label_md <- function(rn) {
  mp <- c(
    "K1 full"      = "GP⁓₁",
    "K2 full"      = "GP⁓₂",
    "K3 full"      = "GP⁓₃",
    "K4 full"      = "GP⁓₄",
    "K5 full"      = "GP⁓₅",
    "K5 additive"  = "GP⁺₅",
    "K5 product"   = "GP×₅",
    "K5 gen-only"  = "GPᴳ₅",
    "K5 env-only"  = "GPᴱ₅",
    "LMM model 1"  = "LMM⁓₁",
    "BGLR model 3" = "BGLR⁓₅",
    "Kernel LMM"   = "kerLMM⁓₅",
    "Global average"        = "GLO_A",
    "Variety average"       = "VAR_A",
    "Environmental average" = "ENV_A"
  )
  if (!is.na(mp[rn])) mp[rn] else rn
}

# ---------- formatting ----------
fmt_num_md <- function(x, digits = 2, na_str = "—") {
  if (is.na(x)) return(na_str)
  formatC(x, format = "f", digits = digits)
}
fmt_pair_md <- function(a, b, digits = 2, na_str = "—") {
  paste0(fmt_num_md(a, digits, na_str), " | ", fmt_num_md(b, digits, na_str))
}

# "two-line in one cell" via <br> (GitHub markdown)
fmt_pair_md_two_lines <- function(non_a, non_b, noisy_a, noisy_b, digits = 2, na_str = "—") {
  paste0(
    fmt_pair_md(non_a, non_b, digits, na_str),
    "<br>(",
    fmt_pair_md(noisy_a, noisy_b, digits, na_str),
    ")"
  )
}

# ---------- Markdown writer (prevents accidental newlines in cells) ----------
write_md_table <- function(df, file) {
  # remove any literal newlines inside cells (breaks markdown tables)
  df[] <- lapply(df, function(x) gsub("[\r\n]+", " ", x))

  # escape literal pipes inside cells
  esc <- function(x) gsub("\\|", "\\\\|", x)
  df[] <- lapply(df, esc)

  header <- paste(names(df), collapse = " | ")
  sep    <- paste(rep("---", ncol(df)), collapse = " | ")
  rows   <- apply(df, 1, function(r) paste(r, collapse = " | "))

  writeLines(c(header, sep, rows), file)
}

## =========================================================
## Noise-aware GP reader (for "noisy_sd_" files)
## =========================================================
kernel_path2 <- function(prefix, tar, keri, leakage, setup, id, noisy_sd = FALSE) {
  paste0(prefix, "_", tar, keri, leakage, setup, "_",
         ifelse(noisy_sd, "noisy_sd_", ""),
         id, ".rda")
}
read_kernel_vec2 <- function(prefix, tar, keri, leakage, setup, id, noisy_sd = FALSE) {
  fp <- kernel_path2(prefix, tar, keri, leakage, setup, id, noisy_sd)
  if (!file.exists(fp)) return(NULL)
  e <- load_env(fp)
  if (exists("res", envir = e, inherits = FALSE)) return(as_numvec(get("res", envir = e)))
  extract_results_vec(e, tar)
}
compute_averages2 <- function(tar, leakage, setup, noisy_sd = FALSE) {
  R <- NULL
  for (keri in 1:6) {
    mats <- lapply(IDS, function(id) {
      fp <- kernel_path2(alg, tar, keri, leakage, setup, id, noisy_sd)
      if (!file.exists(fp)) return(NULL)
      e <- load_env(fp)
      v <- if (exists("res", envir = e, inherits = FALSE)) as_numvec(get("res", envir = e)) else extract_results_vec(e, tar)
      if (is.null(v) || length(v) < 8) return(NULL)
      v[1:8]
    })
    mats <- mats[!vapply(mats, is.null, logical(1))]
    if (length(mats) == 0) next
    R <- rbind(R, do.call(rbind, mats))
  }
  if (is.null(R)) {
    return(list(Global=c(MSE=NA_real_,CRPS=NA_real_),
                Variety=c(MSE=NA_real_,CRPS=NA_real_),
                Env=c(MSE=NA_real_,CRPS=NA_real_)))
  }
  vv <- apply(R, 2, AGG)
  list(
    Global  = c(MSE=vv[3], CRPS=vv[4]),
    Variety = c(MSE=vv[5], CRPS=vv[6]),
    Env     = c(MSE=vv[7], CRPS=vv[8])
  )
}

build_gp_table_for_setup <- function(setup, noisy_sd = FALSE) {

  rows <- c(
    "K1 full","K2 full","K3 full","K4 full","K5 full",
    "K5 additive","K5 product","K5 env-only","K5 gen-only",
    "Global average","Environmental average","Variety average"
  )

  cols <- c(
    "Yield MSE | no leak","Yield MSE | leak",
    "Yield CRPS | no leak","Yield CRPS | leak",
    "Yield LogS | no leak","Yield LogS | leak",
    "Prot MSE | no leak","Prot MSE | leak",
    "Prot CRPS | no leak","Prot CRPS | leak",
    "Prot LogS | no leak","Prot LogS | leak"
  )

  TAB <- matrix(NA_real_, nrow=length(rows), ncol=length(cols),
                dimnames=list(rows, cols))

  for (leakage in c("no","yes")) {

    for (k in 1:5) {
      y_vecs <- lapply(IDS, function(id) read_kernel_vec2(Adam_full, "yield", k, leakage, setup, id, noisy_sd))
      p_vecs <- lapply(IDS, function(id) read_kernel_vec2(Adam_full, "prot",  k, leakage, setup, id, noisy_sd))
      TAB <- fill_metrics(TAB, paste0("K",k," full"), "yield", leakage, summarise_from_vecs(y_vecs))
      TAB <- fill_metrics(TAB, paste0("K",k," full"), "prot",  leakage, summarise_from_vecs(p_vecs))
    }

    y_vecs <- lapply(IDS, function(id) read_kernel_vec2(Adam_add,  "yield", 5, leakage, setup, id, noisy_sd))
    p_vecs <- lapply(IDS, function(id) read_kernel_vec2(Adam_add,  "prot",  5, leakage, setup, id, noisy_sd))
    TAB <- fill_metrics(TAB, "K5 additive", "yield", leakage, summarise_from_vecs(y_vecs))
    TAB <- fill_metrics(TAB, "K5 additive", "prot",  leakage, summarise_from_vecs(p_vecs))

    y_vecs <- lapply(IDS, function(id) read_kernel_vec2(Adam_prod, "yield", 5, leakage, setup, id, noisy_sd))
    p_vecs <- lapply(IDS, function(id) read_kernel_vec2(Adam_prod, "prot",  5, leakage, setup, id, noisy_sd))
    TAB <- fill_metrics(TAB, "K5 product", "yield", leakage, summarise_from_vecs(y_vecs))
    TAB <- fill_metrics(TAB, "K5 product", "prot",  leakage, summarise_from_vecs(p_vecs))

    k_envgen <- 1
    y_vecs <- lapply(IDS, function(id) read_kernel_vec2(Adam_env, "yield", k_envgen, leakage, setup, id, noisy_sd))
    p_vecs <- lapply(IDS, function(id) read_kernel_vec2(Adam_env, "prot",  k_envgen, leakage, setup, id, noisy_sd))
    TAB <- fill_metrics(TAB, "K5 env-only", "yield", leakage, summarise_from_vecs(y_vecs))
    TAB <- fill_metrics(TAB, "K5 env-only", "prot",  leakage, summarise_from_vecs(p_vecs))

    y_vecs <- lapply(IDS, function(id) read_kernel_vec2(Adam_gen, "yield", k_envgen, leakage, setup, id, noisy_sd))
    p_vecs <- lapply(IDS, function(id) read_kernel_vec2(Adam_gen, "prot",  k_envgen, leakage, setup, id, noisy_sd))
    TAB <- fill_metrics(TAB, "K5 gen-only", "yield", leakage, summarise_from_vecs(y_vecs))
    TAB <- fill_metrics(TAB, "K5 gen-only", "prot",  leakage, summarise_from_vecs(p_vecs))

    avY <- compute_averages2("yield", leakage, setup, noisy_sd)
    avP <- compute_averages2("prot",  leakage, setup, noisy_sd)

    TAB["Global average", col_idx("yield","MSE",  leakage)] <- avY$Global["MSE"]
    TAB["Global average", col_idx("yield","CRPS", leakage)] <- avY$Global["CRPS"]
    TAB["Global average", col_idx("prot","MSE",   leakage)] <- avP$Global["MSE"]
    TAB["Global average", col_idx("prot","CRPS",  leakage)] <- avP$Global["CRPS"]

    if (leakage == "yes" || setup == 1) {
      TAB["Variety average", col_idx("yield","MSE",  leakage)] <- avY$Variety["MSE"]
      TAB["Variety average", col_idx("yield","CRPS", leakage)] <- avY$Variety["CRPS"]
      TAB["Variety average", col_idx("prot","MSE",   leakage)] <- avP$Variety["MSE"]
      TAB["Variety average", col_idx("prot","CRPS",  leakage)] <- avP$Variety["CRPS"]
    }
    if (leakage == "yes" || setup == 2) {
      TAB["Environmental average", col_idx("yield","MSE",  leakage)] <- avY$Env["MSE"]
      TAB["Environmental average", col_idx("yield","CRPS", leakage)] <- avY$Env["CRPS"]
      TAB["Environmental average", col_idx("prot","MSE",   leakage)] <- avP$Env["MSE"]
      TAB["Environmental average", col_idx("prot","CRPS",  leakage)] <- avP$Env["CRPS"]
    }
  }

  round(TAB, 3)
}


save_readme_table_md_inline_noisy <- function(TAB_full, setup, out_dir="Results", digits=2, na_str="—") {

  row_order <- c(
    "K1 full","K2 full","K4 full","K5 full",
    "K5 additive","K5 product","K5 gen-only","K5 env-only",
    "LMM model 1","BGLR model 3","Kernel LMM",
    "Global average","Variety average","Environmental average"
  )

  keep <- row_order[row_order %in% rownames(TAB_full)]
  TAB_full <- TAB_full[keep, , drop = FALSE]

  # noisy results exist only for GP rows (and maybe your averages, but you asked GP models)
  TAB_noisy_gp <- build_gp_table_for_setup(setup, noisy_sd = TRUE)

  gp_rows <- c("K1 full","K2 full","K3 full","K4 full","K5 full",
               "K5 additive","K5 product","K5 gen-only","K5 env-only")
  gp_rows <- intersect(gp_rows, rownames(TAB_full))
  gp_rows <- intersect(gp_rows, rownames(TAB_noisy_gp))

  mk_cell <- function(r, col_non, col_noisy) {
    if (r %in% gp_rows) {
      fmt_pair_md_two_lines(
        non_a   = TAB_full[r, col_non[1]],  non_b   = TAB_full[r, col_non[2]],
        noisy_a = TAB_noisy_gp[r, col_noisy[1]], noisy_b = TAB_noisy_gp[r, col_noisy[2]],
        digits = digits, na_str = na_str
      )
    } else {
      fmt_pair_md(TAB_full[r, col_non[1]], TAB_full[r, col_non[2]], digits = digits, na_str = na_str)
    }
  }

  rows <- rownames(TAB_full)

  df <- data.frame(
    Method  = vapply(rows, method_label_md, character(1)),
    `MSE Y`  = vapply(rows, \(r) mk_cell(r,
                                        c("Yield MSE | no leak","Yield MSE | leak"),
                                        c("Yield MSE | no leak","Yield MSE | leak")), character(1)),
    `CRPS Y` = vapply(rows, \(r) mk_cell(r,
                                        c("Yield CRPS | no leak","Yield CRPS | leak"),
                                        c("Yield CRPS | no leak","Yield CRPS | leak")), character(1)),
    `logS Y` = vapply(rows, \(r) mk_cell(r,
                                        c("Yield LogS | no leak","Yield LogS | leak"),
                                        c("Yield LogS | no leak","Yield LogS | leak")), character(1)),
    `MSE P`  = vapply(rows, \(r) mk_cell(r,
                                        c("Prot MSE | no leak","Prot MSE | leak"),
                                        c("Prot MSE | no leak","Prot MSE | leak")), character(1)),
    `CRPS P` = vapply(rows, \(r) mk_cell(r,
                                        c("Prot CRPS | no leak","Prot CRPS | leak"),
                                        c("Prot CRPS | no leak","Prot CRPS | leak")), character(1)),
    `logS P` = vapply(rows, \(r) mk_cell(r,
                                        c("Prot LogS | no leak","Prot LogS | leak"),
                                        c("Prot LogS | no leak","Prot LogS | leak")), character(1)),
    check.names = FALSE
  )

  out_file <- file.path(out_dir, paste0("README_table_setup", setup, ".md"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  write_md_table(df, out_file)
  message("Saved README table: ", out_file)
  invisible(out_file)
}

## ----------------------------
## Run
## ----------------------------
TAB1 <- build_table_for_setup(1)
TAB2 <- build_table_for_setup(2)

save_readme_table_md_inline_noisy(TAB1, setup = 1)
save_readme_table_md_inline_noisy(TAB2, setup = 2)

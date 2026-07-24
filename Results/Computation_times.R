suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

setup   <- 1
leakage <- "no"
tars    <- c("yield", "prot")
IDS     <- 1:30

method_levels_gp <- c(
  "GP 1 full","GP 2 full","GP 4 full","GP 5 full",
  "GP 5 additive","GP 5 product","GP 5 env-only","GP 5 gen-only"
)
method_levels_all <- c(method_levels_gp, "BGLR 3", "kerLMM")

load_res <- function(fp) {
  if (!file.exists(fp)) return(NULL)
  e <- new.env(parent = emptyenv())
  load(fp, envir = e)
  if (!exists("res", envir = e, inherits = FALSE)) return(NULL)
  as.numeric(get("res", envir = e))
}
mmean <- function(x) if (length(x) == 0) NA_real_ else mean(x, na.rm = TRUE)

dir.create("Results/TimePlots", showWarnings = FALSE, recursive = TRUE)

rows_gp <- list()

for (tar in tars) {

  # ---- GP full kernels 1,2,4,5 ----
  for (keri in c(1, 2, 4, 5)) {
    opt <- inf <- c()
    for (id in IDS) {
      fp <- paste0(
        "Results/Results Adam/Results_hyper_ADAM_0.01_1000_0.5_",
        tar, keri, leakage, setup, "_noisy_sd_", id, ".rda"
      )
      res <- load_res(fp)
      if (is.null(res) || length(res) < 26) next
      opt <- c(opt, res[24])
      inf <- c(inf, res[25])
    }
    rows_gp[[length(rows_gp) + 1]] <- data.frame(
      setup = setup, leakage = leakage, tar = tar,
      method = paste0("GP ", keri, " full"),
      Optimization = mmean(opt), Inference = mmean(inf),
      stringsAsFactors = FALSE
    )
  }

  # ---- GP 5 additive/product/env/gen ----
  add_paths <- list(
    "GP 5 additive" = "Results/Results_additive_only/Results Adam",
    "GP 5 product"  = "Results/Results_product_only/Results Adam",
    "GP 5 env-only" = "Results/Results_env_only/Results Adam",
    "GP 5 gen-only" = "Results/Results_gen_only/Results Adam"
  )

  for (m in names(add_paths)) {
    opt <- inf <- c()
    for (id in IDS) {
      fp <- paste0(
        add_paths[[m]], "/Results_hyper_ADAM_0.01_1000_0.5_",
        tar, 5, leakage, setup, "_noisy_sd_", id, ".rda"
      )
      res <- load_res(fp)
      if (is.null(res) || length(res) < 26) next
      opt <- c(opt, res[24])
      inf <- c(inf, res[25])
    }
    rows_gp[[length(rows_gp) + 1]] <- data.frame(
      setup = setup, leakage = leakage, tar = tar,
      method = m,
      Optimization = mmean(opt), Inference = mmean(inf),
      stringsAsFactors = FALSE
    )
  }
}

df_gp <- bind_rows(rows_gp) %>%
  mutate(
    method  = factor(method, levels = method_levels_gp),
    tar     = factor(tar, levels = c("yield", "prot")),
    leakage = factor(leakage, levels = c("no")),
    setup   = factor(setup, levels = c("1"))
  )

write.csv(df_gp, "Results/TimePlots/time_components_mean_GPonly_setup1_leakno.csv", row.names = FALSE)

df_gp_long <- df_gp %>%
  pivot_longer(c(Optimization, Inference), names_to = "phase", values_to = "time") %>%
  mutate(phase = factor(phase, levels = c("Optimization", "Inference")))

p_gp <- ggplot(df_gp_long, aes(x = method, y = time, fill = phase)) +
  geom_col(width = 0.78, na.rm = TRUE) +
  facet_wrap(~ tar, nrow = 1, scales = "free_y") +
  labs(
    title = "Computation time by GP model (setup 1, leakage no)",
    x = NULL, y = "Time (seconds)", fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    panel.spacing = unit(1.1, "lines")
  )

ggsave("Results/TimePlots/time_setup1_no_GP.png",       p_gp, width = 14, height = 5, dpi = 260)
ggsave("Results/TimePlots/time_setup1_no_GP_white.png", p_gp, width = 14, height = 5, dpi = 260, bg = "white")
ggsave("Results/TimePlots/time_setup1_no_GP.pdf",       p_gp, width = 14, height = 5)

print(p_gp)

# ============================================================
# (B) All models: single bars, colored like "Inference"
# Rule:
#   - GP models: use ONLY inference time (res[25])
#   - BGLR / kerLMM: keep total time (res[26])
# Fill: force everything to phase="Inference" (same color)
# ============================================================
rows_all <- list()

for (tar in tars) {

  # ---- GP full kernels 1,2,4,5 : inference only ----
  for (keri in c(1, 2, 4, 5)) {
    inf <- c()
    for (id in IDS) {
      fp <- paste0(
        "Results/Results Adam/Results_hyper_ADAM_0.01_1000_0.5_",
        tar, keri, leakage, setup, "_noisy_sd_", id, ".rda"
      )
      res <- load_res(fp)
      if (is.null(res) || length(res) < 26) next
      inf <- c(inf, res[25])
    }
    rows_all[[length(rows_all) + 1]] <- data.frame(
      setup = setup, leakage = leakage, tar = tar,
      method = paste0("GP ", keri, " full"),
      time = mmean(inf),
      phase = "Inference",
      stringsAsFactors = FALSE
    )
  }

  # ---- GP 5 additive/product/env/gen : inference only ----
  add_paths <- list(
    "GP 5 additive" = "Results/Results_additive_only/Results Adam",
    "GP 5 product"  = "Results/Results_product_only/Results Adam",
    "GP 5 env-only" = "Results/Results_env_only/Results Adam",
    "GP 5 gen-only" = "Results/Results_gen_only/Results Adam"
  )

  for (m in names(add_paths)) {
    inf <- c()
    for (id in IDS) {
      fp <- paste0(
        add_paths[[m]], "/Results_hyper_ADAM_0.01_1000_0.5_",
        tar, 5, leakage, setup, "_noisy_sd_", id, ".rda"
      )
      res <- load_res(fp)
      if (is.null(res) || length(res) < 26) next
      inf <- c(inf, res[25])
    }
    rows_all[[length(rows_all) + 1]] <- data.frame(
      setup = setup, leakage = leakage, tar = tar,
      method = m,
      time = mmean(inf),
      phase = "Inference",
      stringsAsFactors = FALSE
    )
  }

  # ---- BGLR3 + kerLMM : total only ----
  bglr_tot <- ker_tot <- c()
  for (id in IDS) {
    fp_b <- paste0(
      "Results/Results BGLR/Results_hyper_ADAM_0.1_100_3_",
      tar, 5, leakage, setup, "_", id, ".rda"
    )
    rb <- load_res(fp_b)
    if (!is.null(rb) && length(rb) >= 26) bglr_tot <- c(bglr_tot, rb[26])

    fp_k <- paste0(
      "Results/Results BGLR/Results_hyper_ADAM_0.1_100_kerlmm_",
      tar, 5, leakage, setup, "_", id, ".rda"
    )
    rk <- load_res(fp_k)
    if (!is.null(rk) && length(rk) >= 26) ker_tot <- c(ker_tot, rk[26])
  }

  rows_all[[length(rows_all) + 1]] <- data.frame(
    setup = setup, leakage = leakage, tar = tar,
    method = "BGLR 3",
    time = mmean(bglr_tot),
    phase = "Inference",
    stringsAsFactors = FALSE
  )
  rows_all[[length(rows_all) + 1]] <- data.frame(
    setup = setup, leakage = leakage, tar = tar,
    method = "kerLMM",
    time = mmean(ker_tot),
    phase = "Inference",
    stringsAsFactors = FALSE
  )
}

df_all <- bind_rows(rows_all) %>%
  mutate(
    method  = factor(method, levels = method_levels_all),
    tar     = factor(tar, levels = c("yield", "prot")),
    leakage = factor(leakage, levels = c("no")),
    setup   = factor(setup, levels = c("1")),
    phase   = factor(phase, levels = c("Inference"))
  )

write.csv(df_all, "Results/TimePlots/time_inferenceColor_ALL_setup1_leakno.csv", row.names = FALSE)

p_all <- ggplot(df_all, aes(x = method, y = time, fill = phase)) +
  geom_col(width = 0.78, na.rm = TRUE) +
  facet_wrap(~ tar, nrow = 1, scales = "free_y") +
  labs(
    title = "Computation time by model (setup 1, leakage no)",
    x = NULL, y = "Time (seconds)", fill = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 28, hjust = 1, vjust = 1),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    panel.spacing = unit(1.1, "lines")
  )

ggsave("Results/TimePlots/time_inferenceGP_setup1_no.png",       p_all, width = 14, height = 5, dpi = 260)
ggsave("Results/TimePlots/time_inferenceGP_setup1_no_white.png", p_all, width = 14, height = 5, dpi = 260, bg = "white")
ggsave("Results/TimePlots/time_inferenceGP_setup1_no.pdf",       p_all, width = 14, height = 5)

print(p_all)

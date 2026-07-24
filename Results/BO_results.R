
library(dplyr)
library(tidyr)
library(ggplot2)

cols <- c(
  "Random" = "#b7122b",
  "Greedy" = "#a9bbc6",
  "EI"     = "#2c67a5"
)
#take_1000=0
#absolute=0
setupp  <- 1
leakage <- "no"
lr      <- 0.001
max_iter<- 1000
targets <- c("yield","prot")

df_list <- vector("list", length(targets))
names(df_list) <- targets

for (tar in targets) {

  Res <- matrix(0, nrow = 30, ncol = 4)

  for (id in 1:30) {
    file <- paste0(
      "Results/Results BO/Results_hyper_ADAM_",
      lr, "_", max_iter, "_", tar, keri, leakage, setupp, "_", train_prop, "_", id, ".rda"
    )

    if (file.exists(file)) {
      e <- new.env(parent = emptyenv())
      load(file, envir = e)
      if (exists("res", envir = e, inherits = FALSE)) {
        Res <- Res + get("res", envir = e)
      }
    }
  }
  if(take_1000){
  for (id in 31:1000) {
    file <- paste0(
      "Results/Results BO/Results_hyper_ADAM_",
      lr, "_", max_iter, "_", tar, keri, leakage, setupp, "_", train_prop, "_", id, ".rda"
    )

    if (file.exists(file)) {
      e <- new.env(parent = emptyenv())
      load(file, envir = e)
      if (exists("res", envir = e, inherits = FALSE)) {
        Res <- rbind(Res, get("res", envir = e))
      }
    }else{cat('file does not exist: ', file, '\n')}
  }
  }
  res_mat <- Res[apply(Res, 1, sum) > 0, , drop = FALSE]
  if (nrow(res_mat) == 0) next

  df <- data.frame(
    Random = res_mat[,3]/ifelse(!absolute,res_mat[,1]/100,1), #- res_mat[,1],
    Greedy = res_mat[,4]/ifelse(!absolute,res_mat[,1]/100,1),# - res_mat[,1],
    EI     = res_mat[,2]/ifelse(!absolute,res_mat[,1]/100,1) #- res_mat[,1]
  )

  df_long <- pivot_longer(df, everything(),
                          names_to = "Method",
                          values_to = "Value")

  df_long$Method <- factor(df_long$Method, levels = c("Random","Greedy","EI"))
  df_long$Target <- tar

  df_list[[tar]] <- df_long
}

df_all <- bind_rows(df_list)

df_all$Target <- factor(df_all$Target, levels = c("yield","prot"),labels = c("Yield","Protein content"))
facet_scales <- ifelse(absolute, "free_y", "fixed")
p <- ggplot(df_all, aes(x = Method, y = Value)) +
  geom_boxplot(fill = "cadetblue",coef = 3) +
  facet_wrap(~Target, ncol = 2, scales = facet_scales) +
  theme_minimal(base_size = 14) +
  labs(x = "", y=paste0("Improvement BSF ",ifelse(absolute,"","(%)")))+
   geom_jitter(width = 0.3, alpha = 0.2, size = 0.03)

ggsave(paste0("Results/BO_plots/BO1_",train_prop,"_jitterpoints",ifelse(take_1000,"_1000",""),ifelse(absolute,"","_percentual"),".png"), p, width = 14, height = 7, dpi = 300, units = "cm")

ggsave(paste0("Results/BO_plots/BO1_",train_prop,"_jitterpoints_white",ifelse(take_1000,"_1000",""),ifelse(absolute,"","_percentual"),".png"), p, width = 14, height = 7, dpi = 300, units = "cm",bg='white')


p <- ggplot(df_all, aes(x = Method, y = Value)) +
  geom_boxplot(fill = "cadetblue",coef = 3) +
  facet_wrap(~Target, ncol = 2, scales = facet_scales) +
  theme_minimal(base_size = 14) +
  labs(x = "", y=paste0("Improvement BSF ",ifelse(absolute,"","(%)")))

ggsave(paste0("Results/BO_plots/BO1_",train_prop,ifelse(take_1000,"_1000",""),ifelse(absolute,"","_percentual"),".png"), p, width = 14, height = 7, dpi = 300, units = "cm")

ggsave(paste0("Results/BO_plots/BO1_",train_prop,"_white",ifelse(take_1000,"_1000",""),ifelse(absolute,"","_percentual"),".png"), p, width = 14, height = 7, dpi = 300, units = "cm",bg='white')

set.seed(1)

df_plot <- df_all %>%
    mutate(Value_plot = Value) %>%
    group_by(Target) %>%
    mutate(
        eps = ifelse(Target == "yield", 0.08, 0.008),
        Value_plot = ifelse(Value == 0, runif(n(), 0, eps), Value)
    ) %>%
    ungroup()

p <- ggplot() +
    
    geom_violin(
        data = subset(df_plot, Method == "Random"),
        aes(x = Method, y = Value_plot, fill = Method),
        trim = 1, scale = "width", #adjust = 1.5, 
      alpha = 0.7, color = NA
    ) +
    
    geom_violin(
        data = subset(df_plot, Method == "Greedy"),
        aes(x = Method, y = Value_plot, fill = Method),
        trim = 1, scale = "width",#adjust = .75,
      alpha = 0.7, color = NA
    ) +
    geom_violin(
        data = subset(df_plot, Method == "EI"),
        aes(x = Method, y = Value_plot, fill = Method),
        trim = 1, scale = "width", #adjust = .75,
      alpha = 0.7, color = NA
    ) +
    geom_boxplot(
        data = df_plot,
        aes(x = Method, y = Value_plot),
        width = 0.16, coef = 3, outlier.shape = NA,
        fill = "white", color = "black"
    ) +
    geom_jitter(width = 0.1, alpha = 1, size = 0.15) +
    facet_wrap(~Target, ncol = 2, scales =  facet_scales) +
    scale_fill_manual(values = cols) +
    theme_minimal(base_size = 14) +
    theme(
        legend.position = "none",
        panel.grid.major.x = element_blank()
    ) +
    labs(x = "", y = paste0("Improvement BSF ",ifelse(absolute,"","(%)")))


ggsave(
  paste0("Results/BO_plots/BO1_", train_prop, "_jitterpoints_violin",ifelse(take_1000,"_1000",""),ifelse(absolute,"","_percentual"),".png"),
  p, width = 14, height = 7, dpi = 300, units = "cm"
)

ggsave(
  paste0("Results/BO_plots/BO1_", train_prop, "_jitterpoints_violin_white",ifelse(take_1000,"_1000",""),ifelse(absolute,"","_percentual"),".png"),
  p, width = 14, height = 7, dpi = 300, units = "cm", bg = "white"
)

p <- ggplot() +
    
    geom_violin(
        data = subset(df_plot, Method == "Random"),
        aes(x = Method, y = Value_plot, fill = Method),
        trim = 1, scale = "width", #adjust = 1.5, 
      alpha = 0.7, color = NA
    ) +
    
    geom_violin(
        data = subset(df_plot, Method == "Greedy"),
        aes(x = Method, y = Value_plot, fill = Method),
        trim = 1, scale = "width", #adjust = .75,
      alpha = 0.7, color = NA
    ) +
    geom_violin(
        data = subset(df_plot, Method == "EI"),
        aes(x = Method, y = Value_plot, fill = Method),
        trim = 1, scale = "width", #adjust = .75, 
      alpha = 0.7, color = NA
    ) +
    geom_boxplot(
        data = df_plot,
        aes(x = Method, y = Value_plot),
        width = 0.16, coef = 3, outlier.shape = NA,
        fill = "white", color = "black"
    ) +
    facet_wrap(~Target, ncol = 2, scales = facet_scales) +
    scale_fill_manual(values = cols) +
    theme_minimal(base_size = 14) +
    theme(
        legend.position = "none",
        panel.grid.major.x = element_blank()
    ) +
    labs(x = "",y= paste0("Improvement BSF ",ifelse(absolute,"","(%)")))

ggsave(
  paste0("Results/BO_plots/BO1_", train_prop, "_violin",ifelse(take_1000,"_1000",""),ifelse(absolute,"","_percentual"),".png"),
  p, width = 14, height = 7, dpi = 300, units = "cm"
)

ggsave(
  paste0("Results/BO_plots/BO1_", train_prop, "_violin_white",ifelse(take_1000,"_1000",""),ifelse(absolute,"","_percentual"),".png"),
  p, width = 14, height = 7, dpi = 300, units = "cm", bg = "white"
)

save(list='df_plot',file=paste0('Results/Results BO/df_plot_',train_prop,'.rda'))


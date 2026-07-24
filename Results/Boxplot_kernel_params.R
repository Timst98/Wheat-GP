lr=0.01;max_iter=1000;batch_size=0.5;leakage='no';setupp=1
algs=c("Results/Results Adam/","Results/Results_additive_only/Results Adam/","Results/Results_product_only/Results Adam/")

for(tar in c('yield','prot')){
   for(leakage in c('no','yes')){ 
     for(setupp in 1:2){
        theta_G=theta_E=ceta=matrix(nrow=30,ncol=4);colnames(theta_G)=colnames(theta_E)=colnames(ceta)=c("~1","+5","x5","~5")
        alpha=beta=matrix(nrow=30,ncol=3);colnames(alpha)=colnames(beta)=c("~1","+5","~5")
        gamma=matrix(nrow=30,ncol=2);colnames(gamma)=c("~1","~4")

        for(id in 1:30){
                i=id        
                keri=1
                filename=paste0(algs[1],"Results_hyper_ADAM_", lr, "_", max_iter, "_", batchsize, "_",tar, keri, leakage, setupp, "_", id, ".rda")
                load(filename)
                theta_G[i,1]=res[14] # opt pars were saved in res[13:17], (theta_E,theta_G,w1,w2,alpha).
                theta_E[i,1]=res[13] 
                ceta[i,1]=res[17] 
        
                w1=res[15];w2=res[16] 
                alpha[i,1]=w1*(1-w2)
                beta[i,1]=w2*w1
                gamma[i,1]=1-w1
                
                keri=5
                filename=paste0(algs[3],"Results_hyper_ADAM_", lr, "_", max_iter, "_", batchsize, "_",tar, keri, leakage, setupp, "_", id, ".rda")
                load(filename)
                theta_G[i,2]=res[14] 
                ceta[i,2]=res[15] 
                theta_E[i,2]=res[13] 
                
                filename=paste0(algs[2],"Results_hyper_ADAM_", lr, "_", max_iter, "_", batchsize, "_",tar, keri, leakage, setupp, "_", id, ".rda")
                load(filename)
                theta_G[i,3]=res[14] 
                ceta[i,3]=res[16] 
                theta_E[i,3]=res[13] 
        
                w1=res[15]
                alpha[i,2]=w1
                beta[i,2]=1-w1
                
                filename=paste0(algs[1],"Results_hyper_ADAM_", lr, "_", max_iter, "_", batchsize, "_",tar, keri, leakage, setupp, "_", id, ".rda")
                load(filename)
                theta_G[i,4]=res[14] 
                ceta[i,4]=res[17] 
                theta_E[i,4]=res[13] 
        
                w1=res[15];w2=res[16] 
                alpha[i,3]=w1*(1-w2)
                beta[i,3]=w2*w1
                gamma[i,2]=1-w1
        }
        
        alpha_old=alpha
        alpha=beta
        beta=alpha_old #since we switched alpha and beta in Eq. 4

        theta_G=theta_G[,c(4,2,3)]
        theta_E=theta_E[,c(4,2,3)]
        ceta=ceta[,c(4,2,3)]
        colnames(theta_G)=colnames(theta_E)=colnames(ceta)=c("~1","+4","x4","~4")[c(4,2,3)]
        alpha=alpha[,c(3,2)]
        beta=beta[,c(3,2)]
        colnames(alpha)=colnames(beta)=c("~1","+4","~4")[c(3,2)]
        
        gamma=as.matrix(gamma[,2])
        colnames(gamma)="~4"

       
        suppressPackageStartupMessages({
          library(dplyr)
          library(tidyr)
          library(ggplot2)
        })
        
        
        suppressPackageStartupMessages(library(patchwork))
        
        mat_long <- function(M, param) {
          stopifnot(is.matrix(M))
          df <- as.data.frame(M)
          df$.rep <- seq_len(nrow(df))
          df %>%
            pivot_longer(-.rep, names_to = "setting", values_to = "value") %>%
            mutate(param = param) %>%
            select(param, setting, value)
        }
        
        
        df_long <- bind_rows(
          mat_long(theta_G, "theta[G]"),
          mat_long(theta_E, "theta[E]"),
          mat_long(alpha,   "alpha"),
          mat_long(beta,    "beta"),
          mat_long(gamma,   "gamma"),
          mat_long(ceta,    "zeta")
        ) %>%
          mutate(
            param = factor(param, levels = c("theta[G]","theta[E]","alpha","beta","gamma","zeta")),
            setting = factor(setting, levels = c("~4","+4","x4"))
          )
        
        pretty_theme <- function() {
          theme_minimal(base_size = 21) +
            theme(
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.grid.major.y = element_line(linewidth = 0.35),
              axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
              axis.title.y = element_text(margin = margin(r = 8)),
              strip.placement = "outside",
              strip.background = element_rect(fill = "grey95", color = NA),
              strip.text.x = element_text(size = 21),
              panel.spacing.x = unit(0.8, "lines"),
              plot.margin = margin(8, 8, 8, 8)
            )
        }
        
        make_plot <- function(df, y_limits) {
          ggplot(df, aes(x = setting, y = value)) +
            geom_boxplot(
              width = 0.55,
              outlier.size = 1.3,
              outlier.stroke = 0.35,
              linewidth = 0.7,
              color = "black"
            ) +
            facet_grid(. ~ param, scales = "free_x", space = "free_x", switch = "x",
                       labeller = label_parsed) +
            coord_cartesian(ylim = y_limits) +
            labs(y = "Value", x = NULL) +
            pretty_theme()
        }
        
        pG <- make_plot(filter(df_long, param == "theta[G]"), y_limits = c(0.8, 4))
        
        
        pRest <- make_plot(filter(df_long, param != "theta[G]"), y_limits = c(0, 1)) +
          theme(axis.title.y = element_blank())  
        
        p <- patchwork::wrap_plots(pG, pRest, nrow = 1, widths = c(1.2, 5.0))
        
        
        png(paste0("Boxplots/Optimized_parameters_",tar,"_",leakage,"_",setupp,".png"), width = 16, height = 4, units = "in", res = 1000)
        print(p)
        dev.off()
       }
     }
}

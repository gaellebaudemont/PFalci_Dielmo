rm(list=ls())
library(ggplot2)
library(ggpubr)
library(base)
library(binom)
library(ggridges)
library(coda)
library(ROCit)


#####################
## DATA MANAGEMENT ##
#####################
# Download Data: Samples and negative controls
dielmo_Ndiop <- read.csv("0.Dielmo_Ndiop_PFalci.csv")
vec_malaria <- names(dielmo_Ndiop)[grepl("MAL_Pf", names(dielmo_Ndiop))]
vec_epidemio <- colnames(dielmo_Ndiop)[1:7]

neg <- read.csv("0.Seroped_nMFI.csv")



# Creation on df for density plot depending on age
cohort_long_u5 <- dielmo_Ndiop[dielmo_Ndiop$age <= 5,]
cohort_long_AgeMid <- dielmo_Ndiop[dielmo_Ndiop$age > 5 & dielmo_Ndiop$age <= 20,]
cohort_long_old <- dielmo_Ndiop[dielmo_Ndiop$age > 20,]

dielmo_Ndiop$Age_group <- "All"
cohort_long_u5$Age_group <- "Children under 5"
cohort_long_AgeMid$Age_group <- "Young adults (5-20)"
cohort_long_old$Age_group <- "Adults (>20)"

cohort_tres_longue <- rbind(dielmo_Ndiop, cohort_long_u5, cohort_long_AgeMid, cohort_long_old)
cohort_tres_longue <- cohort_tres_longue[, c("sample_name", vec_malaria, "Age_group")]        
colnames(cohort_tres_longue)[1] <- "id_sample"
neg <- neg[,c(1,12:21)]
neg$Age_group <- "Negatives"
colnames(neg) <- colnames(cohort_tres_longue)
cohort_tres_longue <- rbind(cohort_tres_longue, neg)
cohort_tres_longue$Age_group <- factor(cohort_tres_longue$Age_group, levels=c("All", "Adults (>20)", "Young adults (5-20)", "Children under 5", "Negatives"))



# Download files with results
files_results <- paste0("results_pooled_", vec_malaria, "_nMFI.csv")
list_results <- lapply(files_results, read.csv)






################
### FOR LOOP ###
################
for(ab in 1:length(vec_malaria)){
  # Select results and data from the considered ag
  dfanalyses <- list_results[[ab]]
    colnames(dfanalyses) <- c("theta", "sigma_N", "mu_N", "sigma_P", "mu_P", "LL")  
    dfanalyses$mu_pos <- dfanalyses$mu_N + dfanalyses$mu_P


    # Select results of interest and distinguish chains
    dfanalyses$chain <- as.character(c(rep(c(1), each = 7500)))
    dfanalyses$iteration <- c(rep(seq(1,7500, 1), 1))
    
    # Plots
    ############################
    ## POSTERIOR DISTRIBUTION ##
    ############################
    plot_chains1 <- list()
    for (param in 1: (dim(dfanalyses)[2] - 2)){
      param_tempo <- colnames(dfanalyses)[param]
      df_tempo <- dfanalyses[,c(param_tempo, "chain")]
      colnames(df_tempo) <- c("parametre", "chain")
      
      plot_chains1[[param]]<- ggplot(df_tempo) + geom_density(aes(x = parametre, group = chain, color = chain), linewidth = 1.5, alpha = 0.2) +
        scale_color_manual(values = c("#0033FF"))+
        geom_vline(aes(xintercept = median(parametre)), color = "black", linetype = "dashed", size = 1.2)+
        geom_text(aes(x =  median(parametre)+0.1*median(parametre), label =  round(median(parametre), 3)), 
                  y = 0.1)+
        theme_classic()+ ggtitle(param_tempo)+
        theme(plot.title = element_text(size=22, face = "bold"),
              axis.text=element_text(size=12),
              axis.title=element_text(size=16, face = "bold"),
              axis.title.x=element_blank())
      
    }
    
    
    plot_chains2 <- list()
    for (param in 1: (dim(dfanalyses)[2] - 2)){
      param_tempo <- colnames(dfanalyses)[param]
      df_tempo <- dfanalyses[,c(param_tempo, "chain", "iteration")]
      colnames(df_tempo) <- c("parametre", "chain", "iteration")
      
      plot_chains2[[param]]<- ggplot(df_tempo) + geom_line(aes(x = iteration, y = parametre, group = chain, color = chain), size = 1.5, alpha = 0.2) +
        scale_color_manual(values = c("#0033FF"))+
        theme_classic()+ 
        theme(plot.title = element_text(size=22, face = "bold"),
              axis.text=element_text(size=12),
              axis.title=element_text(size=16, face = "bold"),
              axis.title.x=element_blank())
      
    }
    
    plot_chains <- ggarrange(plot_chains1[[1]], plot_chains1[[2]], plot_chains1[[3]], plot_chains1[[4]], plot_chains1[[5]], plot_chains1[[6]], plot_chains1[[7]],
                             plot_chains2[[1]], plot_chains2[[2]], plot_chains2[[3]], plot_chains2[[4]], plot_chains2[[5]], plot_chains2[[6]], plot_chains2[[7]],
                             ncol = 7, nrow = 2, common.legend = TRUE)
    
    ggsave(filename = paste0(vec_malaria[ab], "_chains.jpeg"), plot_chains, width=22, height=8)
    
    
      
}






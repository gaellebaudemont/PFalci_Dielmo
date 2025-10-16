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



# Df creation to compare AUC and Se between ag
AUC_df <- as.data.frame(matrix(NA, nrow = 10, ncol = 4))
colnames(AUC_df) <- c("Ab", "AUC", "ICmin", "ICmax")

Se_df <- as.data.frame(matrix(NA, nrow = 10, ncol = 4))
colnames(Se_df) <- c("Ab", "Se", "ICmin", "ICmax")

cutoff_df <- as.data.frame(matrix(NA, nrow = 10, ncol = 7))
colnames(cutoff_df) <- c("Ab", "C_SP99", "C_SP977", "C_SP95", "C_SP90", "C_msd", "C_youden")











################
### FOR LOOP ###
################
for(ab in 1:length(vec_malaria)){
  # Select results and data from the considered ag
    parameters <- list_results[[ab]]
    colnames(parameters) <- c("theta", "sigma_N", "mu_N", "sigma_P", "mu_P", "LL")  
    parameters$mu_pos <- parameters$mu_N + parameters$mu_P

    
    dielmo_tempo <- dielmo_Ndiop[, c(vec_malaria[ab], "age", "year")] 
    colnames(dielmo_tempo)[1] <- 'AB'
    
    neg_tempo <- log(neg[, c(vec_malaria[ab])]) 
    
  
  
  # Median Se and SP
    nMFI <- c( seq(min(log(dielmo_tempo$AB), na.rm = TRUE), max(log(dielmo_tempo$AB)+1, na.rm = TRUE), 0.05))
    estimated_density <- as.data.frame(nMFI)
    
    estimated_density$specificity <-pnorm(estimated_density$nMFI, mean = median(parameters$mu_N), sd = median(parameters$sigma_N))
    estimated_density$FPR <- 1 - estimated_density$specificity
    estimated_density$sensitivity <- 1 - (pnorm(estimated_density$nMFI, mean =median(parameters$mu_pos), sd = median(parameters$sigma_P)))
    
    estimated_density <- estimated_density[order(estimated_density$nMFI, decreasing = TRUE),]
    
    
  # AUC computation for median Se and Sp
    AUC =round(sum(estimated_density$sensitivity[1:length(nMFI)]*diff(c(0, estimated_density$FPR[1:length(nMFI)]))),2)
    
  
    
    
  # Cutoff, Se and Sp for mean + 2sd for the observed negatives 
    cutoff_neg <- round(mean(neg_tempo) + 2*sd(neg_tempo), 2)
 
    estimated_density$index_cutoff <- abs(estimated_density$nMFI - cutoff_neg)
    Se_cutoff <- round(estimated_density$sensitivity[estimated_density$index_cutoff == min(estimated_density$index_cutoff)],2)
    Sp_cutoff <- round(estimated_density$specificity[estimated_density$index_cutoff == min(estimated_density$index_cutoff)],2)
    
    
  # Youden Index computation and affiliated cutoff
    estimated_density$Youden_Index <- estimated_density$sensitivity + estimated_density$specificity - 1
    
    cutoff_ROC <- round(estimated_density$nMFI[estimated_density$Youden_Index==max(estimated_density$Youden_Index)],2)
    Se_youden <- round(estimated_density$sensitivity[estimated_density$Youden_Index==max(estimated_density$Youden_Index)],2)
    Sp_youden <- round(estimated_density$specificity[estimated_density$Youden_Index==max(estimated_density$Youden_Index)],2)
    
    
  # Cutoff and Se for Sp = 99%, SP = 97.7%, SP = 95% and for SP = 90%
    estimated_density$index_SP99 <- abs(estimated_density$specificity - 0.99)
    estimated_density$index_SP97.7 <- abs(estimated_density$specificity - 0.977)
    estimated_density$index_SP95 <- abs(estimated_density$specificity - 0.95)
    estimated_density$index_SP90 <- abs(estimated_density$specificity - 0.90)
    
    Se_Sp99 <- round(estimated_density$sensitivity[estimated_density$index_SP99 == min(estimated_density$index_SP99)],2)
    Se_Sp97.7<- round(estimated_density$sensitivity[estimated_density$index_SP97.7 == min(estimated_density$index_SP97.7)],2)
    Se_Sp95 <- round(estimated_density$sensitivity[estimated_density$index_SP95 == min(estimated_density$index_SP95)],2)
    Se_Sp90 <- round(estimated_density$sensitivity[estimated_density$index_SP90 == min(estimated_density$index_SP90)],2)
    
    cutoff_Sp99 <- round(estimated_density$nMFI[estimated_density$index_SP99 == min(estimated_density$index_SP99)],2)
    cutoff_Sp97.7 <- round(estimated_density$nMFI[estimated_density$index_SP97.7 == min(estimated_density$index_SP97.7)],2)
    cutoff_Sp95 <- round(estimated_density$nMFI[estimated_density$index_SP95 == min(estimated_density$index_SP95)],2)
    cutoff_Sp90 <- round(estimated_density$nMFI[estimated_density$index_SP90 == min(estimated_density$index_SP90)],2)
    

    
    
    # IC Se, Sp and AUC
    N_sample = dim(parameters)[1]/10   # Select and keep 10% of posterior distribution (could be more)
    sam_seq = round(seq(from=1, to=nrow(parameters), length=N_sample))
    parameters_IC <- parameters[sam_seq,]
    
    df_Se_IC <- as.data.frame(matrix(NA, nrow = N_sample, ncol = length(nMFI)))
    df_Sp_IC <- as.data.frame(matrix(NA, nrow = N_sample, ncol = length(nMFI)))
    df_FPR_IC <- as.data.frame(matrix(NA, nrow = N_sample, ncol = length(nMFI)))
    for (i in 1:N_sample){
      df_Sp_IC[i,] <- pnorm(estimated_density$nMFI, mean = parameters_IC[i,"mu_N"], sd = parameters_IC[i,"sigma_N"])
      df_Se_IC[i,] <- 1 - pnorm(estimated_density$nMFI, mean = parameters_IC[i,"mu_pos"], sd = parameters_IC[i,"sigma_P"])
    }
    
    
    df_Se_IC <- as.data.frame(t(df_Se_IC))
    df_Se_IC$nMFI <- nMFI
    
    df_Sp_IC <- as.data.frame(t(df_Sp_IC))
    df_Sp_IC$nMFI <- nMFI
    
    
    df_IC <- merge(df_Sp_IC, df_Se_IC, by = "nMFI") # df with each coloumn being the Se or Sp of one iteration for each log nMFI value
    
    
    # Compute CI for AUC and Se
    AUC_IC <- vector()
    for (j in 2:(N_sample+1)){
      tempo <- df_IC[,c(j, j+N_sample)]
      AUC_IC[j] =round(sum(tempo[,2]*diff(c(0, 1 - tempo[,1]))),2)
    }
  
    AUC_ICmin <- quantile(AUC_IC, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)[[1]]
    AUC_ICmax <- quantile(AUC_IC, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)[[3]]
      
    
    Se_IC <- vector()
    for (j in 2:(N_sample+1)){
      df_IC$index_SP99 <- abs(df_IC[,j] - 0.99)
      
      Se_IC[j] <- df_IC[,j+N_sample][df_IC$index_SP99 == min(df_IC$index_SP99)]
    }
    
    Se_ICmin <- quantile(Se_IC, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)[[1]]
    Se_ICmax <- quantile(Se_IC, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)[[3]]
    
    
    
    
    
    
    
      ####################
      ## PLOT ROC CURVE ##
      ####################
      
      ROC_plot <- vector('list') # list to add the curve of each iterations
      
    # empty plot
      ROC_plot[[1]] <- ggplot(data = estimated_density) + 
        ggtitle(paste0("A.ROC curve ", vec_malaria[ab]))+
        annotate("rect", xmin = 0.12, xmax = 0.6, ymin = 0.01, ymax = 0.07, color = 'black', fill = "white") +
        annotate("text", x = 0.35, y = 0.043, label = paste0("AUC = ", AUC, " [", AUC_ICmin, ";", AUC_ICmax, "]"), size = 6)+
        xlab("1 - Specificity")+ylab("Sensitivity")+theme_classic()+
        theme(plot.title = element_text(size=18, face = "bold"),
              axis.text=element_text(size=12),
              axis.title = element_text(size=14, face = "bold"),
              legend.position = "none")
      
      # add curve for each iteraion
        for (i in 2:(N_sample+1)){
          tempo <- df_IC[,c(i, i+N_sample)]
          colnames(tempo) <- c("x", "y")
          ROC_plot[[i]] <- ROC_plot[[i - 1]] + geom_line(data = tempo, aes(x = 1 - x, y = y), linewidth = 1, color = "grey", alpha = 0.1)
        }
      
      # add median curve
      ROC_plot_tempo <- ROC_plot[[N_sample+1]] + geom_line(data = estimated_density, aes(x = FPR, y = sensitivity), size = 2)
      
      # add point for cutoffs and "legend"
      ROC_plot_final <- ROC_plot_tempo + annotate("point", x = 1 - Sp_youden, y = Se_youden, size = 3, color = "red")+
        annotate("rect", xmin = 0.75, xmax = 0.99, ymin = 0.83, ymax = 0.98, color = 'red', fill = "white") +
        annotate("text", x = 0.87, y = 0.91, label = paste0("Youden Index: \n Cutoff = ", cutoff_ROC, ",\n  Se = ", Se_youden, ",\n  Sp = ", Sp_youden), size = 4, color = "red") +
        annotate("point", x = 1 - 0.99, y = Se_Sp99, size = 3, color = "#0033FF")+
        annotate("rect", xmin = 0.75, xmax = 0.99, ymin = 0.01, ymax = 0.14, color = '#0033FF', fill = "white") +
        annotate("text", x = 0.87, y = 0.08, label = paste0("Cutoff = ", cutoff_Sp99, ",\n  Se = ", Se_Sp99, ",\n  Sp = ", 0.99), size = 4, color = "#0033FF")+
        annotate("point", x = 1 - 0.977, y = Se_Sp97.7, size = 3, color = "#0066FF")+
        annotate("rect", xmin = 0.75, xmax = 0.99, ymin = 0.16, ymax = 0.29, color = '#0066FF', fill = "white") +
        annotate("text", x = 0.87, y = 0.23, label = paste0("Cutoff = ", cutoff_Sp97.7, ",\n  Se = ", Se_Sp97.7, ",\n  Sp = ", 0.977), size = 4, color = "#0066FF")+
        annotate("point", x = 1 - 0.95, y = Se_Sp95, size = 3, color = "#0099FF")+
        annotate("rect", xmin = 0.75, xmax = 0.99, ymin = 0.31, ymax = 0.44, color = '#0099FF', fill = "white") +
        annotate("text", x = 0.87, y = 0.38, label = paste0("Cutoff = ", cutoff_Sp95, ",\n  Se = ", Se_Sp95, ",\n  Sp = ", 0.95), size = 4, color = "#0099FF")+
        annotate("point", x = 1 - 0.90, y = Se_Sp90, size = 3, color = "#0099CC")+
        annotate("rect", xmin = 0.75, xmax = 0.99, ymin = 0.46, ymax = 0.59, color = '#0099CC', fill = "white") +
        annotate("text", x = 0.87, y = 0.53, label = paste0("Cutoff = ", cutoff_Sp90, ",\n  Se = ", Se_Sp90, ",\n  Sp = ", 0.90), size = 4, color = "#0099CC")+
        annotate("point", x = 1 - Sp_cutoff, y = Se_cutoff, size = 3, color = "#CC6633")+
        annotate("rect", xmin = 0.75, xmax = 0.99, ymin = 0.66, ymax = 0.81, color = '#CC6633', fill = "white") +
        annotate("text", x = 0.87, y = 0.74, label = paste0("mean + 2sd: \n Cutoff = ", cutoff_neg, ",\n  Se = ", Se_cutoff, ",\n  Sp = ", Sp_cutoff), size = 4, color = "#CC6633")
      
        
      
      # fill up the df for AUC and Se comparison
      AUC_df[ab,] <- c(vec_malaria[ab], AUC, AUC_ICmin, AUC_ICmax)
      Se_df[ab,] <- c(vec_malaria[ab], Se_Sp99, Se_ICmin, Se_ICmax)
      cutoff_df[ab,] <- c(vec_malaria[ab], 
                          estimated_density$nMFI[estimated_density$index_SP99 == min(estimated_density$index_SP99)],
                          estimated_density$nMFI[estimated_density$index_SP97.7 == min(estimated_density$index_SP97.7)],
                          estimated_density$nMFI[estimated_density$index_SP95 == min(estimated_density$index_SP95)],
                          estimated_density$nMFI[estimated_density$index_SP90 == min(estimated_density$index_SP90)],
                          cutoff_neg,
                          estimated_density$nMFI[estimated_density$Youden_Index==max(estimated_density$Youden_Index)])
      
      
      
      
      
      #####################
      ### PLOT SEROPREV ###
      #####################
      age_bins     <- seq(from=0, to=95, by=5)
      age_bins_mid <- seq(from=2.5, to=92.5, by=5) 
      N_bins <- length(age_bins) - 1 
      
      # binarisation
      dielmo_tempo$ab_binaire_youden <- ifelse(log(dielmo_tempo$AB) >= cutoff_ROC, 1, 0)
      dielmo_tempo$ab_binaire_neg <- ifelse(log(dielmo_tempo$AB) >= cutoff_neg, 1, 0)
      dielmo_tempo$ab_binaire_Sp99 <- ifelse(log(dielmo_tempo$AB) >= cutoff_Sp99, 1, 0)
      dielmo_tempo$ab_binaire_Sp97.7 <- ifelse(log(dielmo_tempo$AB) >= cutoff_Sp97.7, 1, 0)
      dielmo_tempo$ab_binaire_Sp95 <- ifelse(log(dielmo_tempo$AB) >= cutoff_Sp95, 1, 0)
      dielmo_tempo$ab_binaire_Sp90 <- ifelse(log(dielmo_tempo$AB) >= cutoff_Sp90, 1, 0)

      
      
      # Cutoff Youden
      SP_bins_youden <- matrix(NA, nrow=N_bins, ncol=4)
      colnames(SP_bins_youden) <- c("med", "low_95", "high_95", "N")
      
      for(i in 1:N_bins){
        index <- which( dielmo_tempo$age >age_bins[i] & dielmo_tempo$age<=age_bins[i+1] ) 
        temp_AB  <- dielmo_tempo[index,'ab_binaire_youden']                                         
        
        SP_bins_youden[i,] <- c(as.numeric(as.vector(
          binom.confint( sum(temp_AB), length(temp_AB), method="wilson")[1,4:6])), length(temp_AB))}
      
      SP_bins_youden <- as.data.frame(SP_bins_youden)
      SP_bins_youden$age <- age_bins_mid
      SP_bins_youden$Cutoff <- "Youden Index"
      
    
      
      # Cutoff GMM
      SP_bins_neg <- matrix(NA, nrow=N_bins, ncol=4)
      colnames(SP_bins_neg) <- c("med", "low_95", "high_95", "N")
      
      for(i in 1:N_bins){
        index <- which( dielmo_tempo$age >age_bins[i] & dielmo_tempo$age<=age_bins[i+1] ) 
        temp_AB  <- dielmo_tempo[index,'ab_binaire_neg']                                           
        
        SP_bins_neg[i,] <- c(as.numeric(as.vector(
          binom.confint( sum(temp_AB), length(temp_AB), method="wilson")[1,4:6])), length(temp_AB))}
      
      SP_bins_neg <- as.data.frame(SP_bins_neg)
      SP_bins_neg$age <- age_bins_mid
      SP_bins_neg$Cutoff <- "mean+2sd"
      
      
      
      # Cutoff Sp = 99
      SP_bins_Sp99 <- matrix(NA, nrow=N_bins, ncol=4)
      colnames(SP_bins_Sp99) <- c("med", "low_95", "high_95", "N")
      
      for(i in 1:N_bins){
        index <- which( dielmo_tempo$age >age_bins[i] & dielmo_tempo$age<=age_bins[i+1] ) 
        temp_AB  <- dielmo_tempo[index,'ab_binaire_Sp99']                                           
        
        SP_bins_Sp99[i,] <- c(as.numeric(as.vector(
          binom.confint( sum(temp_AB), length(temp_AB), method="wilson")[1,4:6])), length(temp_AB))}
      
      SP_bins_Sp99 <- as.data.frame(SP_bins_Sp99)
      SP_bins_Sp99$age <- age_bins_mid
      SP_bins_Sp99$Cutoff <- "Sp99"
      
      
      
      # Cutoff Sp = 99
      SP_bins_Sp99 <- matrix(NA, nrow=N_bins, ncol=4)
      colnames(SP_bins_Sp99) <- c("med", "low_95", "high_95", "N")
      
      for(i in 1:N_bins){
        index <- which( dielmo_tempo$age >age_bins[i] & dielmo_tempo$age<=age_bins[i+1] ) 
        temp_AB  <- dielmo_tempo[index,'ab_binaire_Sp99']                                          
        
        SP_bins_Sp99[i,] <- c(as.numeric(as.vector(
          binom.confint( sum(temp_AB), length(temp_AB), method="wilson")[1,4:6])), length(temp_AB))}
      
      SP_bins_Sp99 <- as.data.frame(SP_bins_Sp99)
      SP_bins_Sp99$age <- age_bins_mid
      SP_bins_Sp99$Cutoff <- "Sp99"
      
      
      
      
      # Cutoff Sp = 97.7
      SP_bins_Sp97.7 <- matrix(NA, nrow=N_bins, ncol=4)
      colnames(SP_bins_Sp97.7) <- c("med", "low_95", "high_95", "N")
      
      for(i in 1:N_bins){
        index <- which( dielmo_tempo$age >age_bins[i] & dielmo_tempo$age<=age_bins[i+1] ) 
        temp_AB  <- dielmo_tempo[index,'ab_binaire_Sp97.7']                                      
        
        SP_bins_Sp97.7[i,] <- c(as.numeric(as.vector(
          binom.confint( sum(temp_AB), length(temp_AB), method="wilson")[1,4:6])), length(temp_AB))}
      
      SP_bins_Sp97.7 <- as.data.frame(SP_bins_Sp97.7)
      SP_bins_Sp97.7$age <- age_bins_mid
      SP_bins_Sp97.7$Cutoff <- "Sp97.7"
      
      
      
      # Cutoff Sp = 95
      SP_bins_Sp95 <- matrix(NA, nrow=N_bins, ncol=4)
      colnames(SP_bins_Sp95) <- c("med", "low_95", "high_95", "N")
      
      for(i in 1:N_bins){
        index <- which( dielmo_tempo$age >age_bins[i] & dielmo_tempo$age<=age_bins[i+1] ) 
        temp_AB  <- dielmo_tempo[index,'ab_binaire_Sp95']                                           
        
        SP_bins_Sp95[i,] <- c(as.numeric(as.vector(
          binom.confint( sum(temp_AB), length(temp_AB), method="wilson")[1,4:6])), length(temp_AB))}
      
      SP_bins_Sp95 <- as.data.frame(SP_bins_Sp95)
      SP_bins_Sp95$age <- age_bins_mid
      SP_bins_Sp95$Cutoff <- "Sp95"
      
      
      
      
      # Cutoff Sp = 90
      SP_bins_Sp90 <- matrix(NA, nrow=N_bins, ncol=4)
      colnames(SP_bins_Sp90) <- c("med", "low_95", "high_95", "N")
      
      for(i in 1:N_bins){
        index <- which( dielmo_tempo$age >age_bins[i] & dielmo_tempo$age<=age_bins[i+1] ) 
        temp_AB  <- dielmo_tempo[index,'ab_binaire_Sp90']                                         
        
        SP_bins_Sp90[i,] <- c(as.numeric(as.vector(
          binom.confint( sum(temp_AB), length(temp_AB), method="wilson")[1,4:6])), length(temp_AB))}
      
      SP_bins_Sp90 <- as.data.frame(SP_bins_Sp90)
      SP_bins_Sp90$age <- age_bins_mid
      SP_bins_Sp90$Cutoff <- "Sp90"
      
      
      
      # Gather all that in one database for plotting
      SP_bins <- rbind(SP_bins_youden, SP_bins_neg, SP_bins_Sp99, SP_bins_Sp97.7, SP_bins_Sp95, SP_bins_Sp90)
      SP_bins$Cutoff <- factor(SP_bins$Cutoff, levels=c("Youden Index", "mean+2sd", "Sp99", "Sp97.7", "Sp95", "Sp90"))
      
      
      plot_seroprev <- ggplot(SP_bins, aes(x = age, y = med, group_by(Cutoff))) +
        geom_point(data = SP_bins, aes(color = Cutoff), size = 2.5,shape = 19, position=position_dodge(2))+
        geom_errorbar(data = SP_bins, aes(ymin=low_95, ymax=high_95, color = Cutoff), width=0,
                      position=position_dodge(2)) +theme_classic()+
        scale_color_manual(values = c("red", "#CC6633", "#0033FF", "#0066FF", "#0099FF", "#0099CC")) +
        ggtitle("D. Seroprevalence based on estimated cutoff") +
        xlab("Age (Years)") + ylab("Seroprevalence")+
        theme(plot.title = element_text(size=22, face = "bold"),
              axis.text=element_text(size=12),
              axis.title=element_text(size=16, face = "bold"))+
        ylim(0,1)
      
      
      
      
      
      ###########################
      ### PLOT AGE VS LOG(AB) ###
      ###########################
       plot_distrib <- ggplot(data = dielmo_tempo, aes(x = age, y = log(AB))) + geom_point(alpha = 0.2, color = "#990000")+
        geom_hline(mapping = aes(yintercept = cutoff_ROC), color = "red", linetype = "dashed", linewidth = 1.2)+
        geom_hline(mapping = aes(yintercept = cutoff_Sp99), color = "#0033FF", linetype = "dashed", size = 1.2)+
        geom_hline(mapping = aes(yintercept = cutoff_Sp97.7), color = "#0066FF", linetype = "dashed", size = 1.2)+
        geom_hline(mapping = aes(yintercept = cutoff_Sp95), color = "#0099FF", linetype = "dashed", size = 1.2)+
        geom_hline(mapping = aes(yintercept = cutoff_Sp90), color = "#0099CC", linetype = "dashed", size = 1.2)+
        geom_hline(mapping = aes(yintercept = cutoff_neg), color = "#CC6633", linetype = "dashed", size = 1.2)+
        ggtitle("C.")+
        xlab("Age")+ylab("log(nMFI)")+theme_classic()+
        theme(plot.title = element_text(size=18, face = "bold"),
              axis.text=element_text(size=12),
              axis.title = element_text(size=14, face = "bold"))#+ ylim(0,10)
      
      
      

        
      
      ###################################
      ### PLOT DISTRIBUTION / DENSITY ###
      ###################################
      # estimated density 
      line_neg_estimated <- c((1 -median(parameters$theta))*(1/(median(parameters$sigma_N)*sqrt(2*pi)))*exp(-(1/2)*((nMFI - median(parameters$mu_N))^2)/(median(parameters$sigma_N))^2))
      line_pos_estimated <- c(median(parameters$theta)*(1/(median(parameters$sigma_P)*sqrt(2*pi)))*exp(-(1/2)*((nMFI - median(parameters$mu_pos))^2)/(median(parameters$sigma_P))^2))
      
      estimated_density <- as.data.frame(cbind(nMFI, line_neg_estimated, line_pos_estimated))
      estimated_density$Age_group <- "Estimated distribution"
      
      
      data_density <- cohort_tres_longue[, c(vec_malaria[ab], "Age_group")] 
      colnames(data_density)[1] <- 'V1'
      
      dens_plot <- ggplot(data = data_density, aes(x = log(V1), y = Age_group, group = Age_group)) + geom_density_ridges(aes(color = Age_group, fill = Age_group), stat="binline", scale = 2, alpha = 0.5)+
        ggtitle("B.")+
        scale_y_discrete(limits = c("Negatives", "Children under 5", "Young adults (5-20)", "Adults (>20)", "All"))+
        scale_fill_manual(values = c("#990000", "#CC3333", "#FF6666", "#FF9999", "#66CC66")) +
        scale_color_manual(values = c("#990000", "#CC3333", "#FF6666", "#FF9999", "#66CC66")) +
        geom_vline(mapping = aes(xintercept = median(parameters$mu_N)), color = "black", linewidth = 1)+
        geom_vline(mapping = aes(xintercept = median(parameters$mu_pos)), color = "black", linewidth = 1)+
        geom_vline(mapping = aes(xintercept = cutoff_ROC), color = "red", linetype = "dashed", linewidth = 0.7)+
        geom_vline(mapping = aes(xintercept = cutoff_Sp99), color = "#0033FF", linetype = "dashed", linewidth = 0.7)+
        geom_vline(mapping = aes(xintercept = cutoff_Sp97.7), color = "#0066FF", linetype = "dashed", linewidth = 0.7)+
        geom_vline(mapping = aes(xintercept = cutoff_Sp95), color = "#0099FF", linetype = "dashed", linewidth = 0.7)+
        geom_vline(mapping = aes(xintercept = cutoff_Sp90), color = "#0099CC", linetype = "dashed", linewidth = 0.7)+
        geom_vline(mapping = aes(xintercept = cutoff_neg), color = "#CC6633", linetype = "dashed", linewidth = 0.7)+
        geom_line(data = estimated_density, mapping = aes(y = line_neg_estimated, x = nMFI))+
        geom_line(data = estimated_density, mapping = aes(y = line_pos_estimated, x = nMFI))+
        xlab("log(nMFI)")+ylab("density")+theme_classic()+
        theme(plot.title = element_text(size=18, face = "bold"),
              axis.text=element_text(size=12),
              axis.title = element_text(size=14, face = "bold"),
              legend.position = "none")#+ xlim(0,13)
      
      
      
      
      
            
      plots_tempo <- ggarrange(plot_distrib, plot_seroprev, ncol = 1) 
      
      plot_final <- ggarrange(ROC_plot_final, dens_plot, plots_tempo, nrow = 1)

      #ggsave(paste0("plot_AUROC_GMM_99_", vec_malaria[ab], ".jpg"), plot_final, width = 22, height = 13)
      
}



write.csv(cutoff_df, "cutoff_df.csv")


# plot AUC ranked
AUC_df <- AUC_df[order(as.numeric(AUC_df$AUC), decreasing = FALSE),]
# To order Ab names 
AUC_df$Ab <- factor(AUC_df$Ab, levels=c(unique(AUC_df$Ab)))


plot_aucs <- ggplot(AUC_df, aes(x = as.numeric(AUC), y = Ab)) + geom_point(size =2) + xlim(0,1)+
  geom_errorbar(aes(xmin=as.numeric(ICmin), xmax=as.numeric(ICmax)), width=0)+
  ggtitle("AUC comparison")+
  xlab("AUC")+ylab("Antigen")+theme_classic()+
  theme(plot.title = element_text(size=18, face = "bold"),
        axis.text=element_text(size=12),
        axis.title = element_text(size=14, face = "bold"),
        axis.title.y = element_blank())
plot_aucs





# plot Se ranked (for SP = 99%)

Se_df <- Se_df[order(as.numeric(Se_df$Se), decreasing = FALSE),]
# To order Ab names 
Se_df$Ab <- factor(Se_df$Ab, levels=c(unique(Se_df$Ab)))


plot_aucs <- ggplot(Se_df, aes(x = as.numeric(Se), y = Ab)) + geom_point(size =2) + xlim(0,1)+
  geom_errorbar(aes(xmin=as.numeric(ICmin), xmax=as.numeric(ICmax)), width=0)+
  ggtitle("Sensitivy comparison for a Specificity of 99%")+
  xlab("Se")+ylab("Antigen")+theme_classic()+
  theme(plot.title = element_text(size=18, face = "bold"),
        axis.text=element_text(size=12),
        axis.title = element_text(size=14, face = "bold"),
        axis.title.y = element_blank())
plot_aucs




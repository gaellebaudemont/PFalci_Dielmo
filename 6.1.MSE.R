# Download packages
rm(list = ls())
library(base)
library(ggplot2)
library(ggpubr)
library(stringr)

# Download validation data
Senegalese <- read.csv(file = "Senegalese_data_NbCases.csv") 
Senegalese <- Senegalese[,-1]
Senegalese$age <- c(seq(26, 4, - 1))


# Download results from modelisation
my_files_3Ab_tempo = list.files(path = "./3Ab_1tc", pattern = "*_dielmo.csv")
my_files_3Ab <- c(paste0("./3Ab_1tc/", my_files_3Ab_tempo))
list_results_3Ab <- lapply(my_files_3Ab, read.csv)
my_files_3Ab_tempo <- str_sub(my_files_3Ab_tempo, end=-16)

my_files_2Ab_tempo = list.files(path = "./2Ab_1tc", pattern = "*_dielmo.csv")
my_files_2Ab <- c(paste0("./2Ab_1tc/", my_files_2Ab_tempo))
list_results_2Ab <- lapply(my_files_2Ab, read.csv)
my_files_2Ab_tempo <- str_sub(my_files_2Ab_tempo, end=-16)

my_files_1Ab_1tc_tempo = list.files(path = "./1Ab_1tc", pattern = "*_dielmo.csv")
my_files_1Ab_1tc <- c(paste0("./1Ab_1tc/", my_files_1Ab_1tc_tempo))
list_results_1Ab_1tc <- lapply(my_files_1Ab_1tc, read.csv)
my_files_1Ab_1tc_tempo <- str_sub(my_files_1Ab_1tc_tempo, end=-16)



# MSE computation
for (combinaison in 1:length(my_files_3Ab)){
  # Ab combination
  Ab1 <- strsplit(my_files_3Ab_tempo[[combinaison]], split = "MAL_")[[1]][2]
  Ab2 <- strsplit(my_files_3Ab_tempo[[combinaison]], split = "MAL_")[[1]][3]
  Ab3 <- strsplit(my_files_3Ab_tempo[[combinaison]], split = "MAL_")[[1]][4]
  
  # parameters vectors : med and IC95
  tempo_med <- c(median(list_results_3Ab[[combinaison]]$lambda), median(list_results_3Ab[[combinaison]]$time_c), median(list_results_3Ab[[combinaison]]$delta))
  tempo_2.5 <- c(quantile(list_results_3Ab[[combinaison]]$lambda, c(0.025)), quantile(list_results_3Ab[[combinaison]]$time_c, c(0.975)), quantile(list_results_3Ab[[combinaison]]$delta, c(0.025)))
  tempo_97.5 <- c(quantile(list_results_3Ab[[combinaison]]$lambda, c(0.975)), quantile(list_results_3Ab[[combinaison]]$time_c, c(0.025)), quantile(list_results_3Ab[[combinaison]]$delta, c(0.975)))
  
  
  # find S to minimise the Standard error
  S<- c(seq(1, 300, by = 0.5))
  MSE <- c(rep(NA, 599))
  Sfinal <- as.data.frame(cbind(S, MSE))
  
  for(j in 1:length(S)){
    tempo <- Senegalese
    
    tempo$Error <- ifelse(tempo$age > tempo_med[2], abs(tempo$nb_cases - S[j]*tempo_med[1])^2, abs(tempo$nb_cases - S[j]*tempo_med[1]*tempo_med[3])^2)
    
    Sfinal[j,"MSE"] <- mean(tempo$Error)
  }
  S_1tc <- Sfinal$S[Sfinal$MSE== min(Sfinal$MSE)]
  
  # MSE computation
  MSE_tempo <- c(rep(NA, 550))
  for (k in 1:length(list_results_3Ab[[combinaison]]$lambda)){
    tempo <- Senegalese
    
    tempo$Error <- ifelse(tempo$age > list_results_3Ab[[combinaison]][k,"time_c"],
                          abs(tempo$nb_cases - S_1tc*list_results_3Ab[[combinaison]][k,"lambda"])^2,
                          abs(tempo$nb_cases - S_1tc*list_results_3Ab[[combinaison]][k,"lambda"]*list_results_3Ab[[combinaison]][k,"delta"])^2)
    
    MSE_tempo[k] <- mean(tempo$Error) 
  }

  MSE_tempo <- as.data.frame(MSE_tempo)
  MSE_tempo$Ab1 <- Ab1
  MSE_tempo$Ab2 <- Ab2
  MSE_tempo$Ab3 <- Ab3
  MSE_tempo$nb_ab <- 3
  
  # save results
  write.csv(MSE_tempo, paste0("MSE_", Ab1, "_", Ab2, "_", Ab3, "_1tc_.csv"))
  

  
  # Plot
  year_sim<- c(seq(1990, 2016, by = 0.1))
  ajusted_lambda <- c(rep(NA, length(year_sim))) 
  ajusted_lambda2.5 <- c(rep(NA, length(year_sim))) 
  ajusted_lambda97.5 <- c(rep(NA, length(year_sim))) 
  
  df_model <- as.data.frame(cbind(year_sim, ajusted_lambda, ajusted_lambda2.5, ajusted_lambda97.5))
  df_model$age <- 2016 - df_model$year_sim
  
  df_model$ajusted_lambda <- ifelse(df_model$age > tempo_med[2],
                                             S_1tc*tempo_med[1],
                                             S_1tc*tempo_med[1]*tempo_med[3])
  
  df_model$ajusted_lambda2.5 <- ifelse(df_model$age > tempo_2.5[2],
                                    S_1tc*tempo_2.5[1],
                                    S_1tc*tempo_2.5[1]*tempo_2.5[3])
  
  df_model$ajusted_lambda97.5 <- ifelse(df_model$age > tempo_97.5[2],
                                    S_1tc*tempo_97.5[1],
                                    S_1tc*tempo_97.5[1]*tempo_97.5[3])

  
  plot_validation <- ggplot(Senegalese, mapping = aes(y = nb_cases, x = year)) + geom_bar(stat = 'identity', width = 0.4, fill="#006699", alpha = 0.8) +
    geom_line(df_model, mapping = aes(y=ajusted_lambda, x = year_sim), color = "#990000", alpha = 1, linewidth = 2, inherit.aes = FALSE)+
    geom_ribbon(df_model, mapping = aes(ymin = ajusted_lambda2.5, ymax = ajusted_lambda97.5, x = year_sim), fill = "#990000", alpha = 0.25, inherit.aes = FALSE)+
    geom_vline(aes(xintercept = 1995), color = "black",  alpha = 1, linewidth = 1, linetype = "dashed")+
    geom_vline(aes(xintercept = 2003), color = "black",  alpha = 1, linewidth = 1, linetype = "dashed")+
    geom_vline(aes(xintercept = 2006), color = "black",  alpha = 1, linewidth = 1, linetype = "dashed")+
    geom_vline(aes(xintercept = 2008), color = "black",  alpha = 1, linewidth = 1, linetype = "dashed")+
    geom_text(mapping=aes(x=1995.2, y=30, label="1"), size=6, vjust=-0.4, hjust=0)+
    geom_text(mapping=aes(x=2003.2, y=30, label="2"), size=6, vjust=-0.4, hjust=0)+
    geom_text(mapping=aes(x=2006.2, y=30, label="3"), size=6, vjust=-0.4, hjust=0)+
    geom_text(mapping=aes(x=2008.2, y=30, label="4"), size=6, vjust=-0.4, hjust=0)+
    theme_light()+ 
    scale_x_continuous(breaks=seq(1990, 2016, 5)) +
    scale_y_continuous(name = "Number of cases per person per year",
                       sec.axis = sec_axis( trans=~./S_1tc, name="Force of infection")) +
    xlab("Year") + ylab("Number of cases per person per year")+  
    labs(title = paste0(str_sub(Ab1, end=-2), "  ", str_sub(Ab2, end=-2), "  ", Ab3))+
    theme(plot.title = element_text(size=22),
          plot.subtitle = element_text(size=16),
          axis.text=element_text(size=12),
          axis.title = element_text(size=16, face = "bold"), 
          plot.caption = element_text(hjust=0))
  
  ggsave(filename = paste0(str_sub(Ab1, end=-2), "_", str_sub(Ab2, end=-2), "_", Ab3, "_validation.jpeg"), plot_validation, width=10, height=8)
  
  
  
}







for (combinaison in 1:length(my_files_2Ab)){
  # Ab combination
  Ab1 <- strsplit(my_files_2Ab_tempo[[combinaison]], split = "MAL_")[[1]][2]
  Ab2 <- strsplit(my_files_2Ab_tempo[[combinaison]], split = "MAL_")[[1]][3]

  
  # parameters vectors : med and IC95
  tempo_med <- c(median(list_results_2Ab[[combinaison]]$lambda), median(list_results_2Ab[[combinaison]]$time_c), median(list_results_2Ab[[combinaison]]$delta))
  tempo_2.5 <- c(quantile(list_results_2Ab[[combinaison]]$lambda, c(0.025)), quantile(list_results_2Ab[[combinaison]]$time_c, c(0.975)), quantile(list_results_2Ab[[combinaison]]$delta, c(0.025)))
  tempo_97.5 <- c(quantile(list_results_2Ab[[combinaison]]$lambda, c(0.975)), quantile(list_results_2Ab[[combinaison]]$time_c, c(0.025)), quantile(list_results_2Ab[[combinaison]]$delta, c(0.975)))
  
  
  # find S to minimise the Standard error
  S<- c(seq(1, 300, by = 0.5))
  MSE <- c(rep(NA, 599))
  Sfinal <- as.data.frame(cbind(S, MSE))
  
  for(j in 1:length(S)){
    tempo <- Senegalese
    
    tempo$Error <- ifelse(tempo$age > tempo_med[2], abs(tempo$nb_cases - S[j]*tempo_med[1])^2, abs(tempo$nb_cases - S[j]*tempo_med[1]*tempo_med[3])^2)
    
    Sfinal[j,"MSE"] <- mean(tempo$Error)
  }
  S_1tc <- Sfinal$S[Sfinal$MSE== min(Sfinal$MSE)]
  
  # MSE computation
  MSE_tempo <- c(rep(NA, 550))
  for (k in 1:length(list_results_2Ab[[combinaison]]$lambda)){
    tempo <- Senegalese
    
    tempo$Error <- ifelse(tempo$age > list_results_2Ab[[combinaison]][k,"time_c"], abs(tempo$nb_cases - S_1tc*list_results_2Ab[[combinaison]][k,"lambda"])^2, abs(tempo$nb_cases - S_1tc*list_results_2Ab[[combinaison]][k,"lambda"]*list_results_2Ab[[combinaison]][k,"delta"])^2)
    
    MSE_tempo[k] <- mean(tempo$Error) 
  }
  
  MSE_tempo <- as.data.frame(MSE_tempo)
  MSE_tempo$Ab1 <- Ab1
  MSE_tempo$Ab2 <- Ab2
  MSE_tempo$Ab3 <- NA
  MSE_tempo$nb_ab <- 2
  
  # save results
  write.csv(MSE_tempo, paste0("MSE_", Ab1, "_", Ab2, "_1tc_.csv"))
  
  
  
  
  # Plot
  year_sim<- c(seq(1990, 2016, by = 0.1))
  ajusted_lambda <- c(rep(NA, length(year_sim))) 
  ajusted_lambda2.5 <- c(rep(NA, length(year_sim))) 
  ajusted_lambda97.5 <- c(rep(NA, length(year_sim))) 
  
  df_model <- as.data.frame(cbind(year_sim, ajusted_lambda, ajusted_lambda2.5, ajusted_lambda97.5))
  df_model$age <- 2016 - df_model$year_sim
  
  df_model$ajusted_lambda <- ifelse(df_model$age > tempo_med[2],
                                    S_1tc*tempo_med[1],
                                    S_1tc*tempo_med[1]*tempo_med[3])
  
  df_model$ajusted_lambda2.5 <- ifelse(df_model$age > tempo_2.5[2],
                                       S_1tc*tempo_2.5[1],
                                       S_1tc*tempo_2.5[1]*tempo_2.5[3])
  
  df_model$ajusted_lambda97.5 <- ifelse(df_model$age > tempo_97.5[2],
                                        S_1tc*tempo_97.5[1],
                                        S_1tc*tempo_97.5[1]*tempo_97.5[3])
  
  
  plot_validation <- ggplot(Senegalese, mapping = aes(y = nb_cases, x = year)) + geom_bar(stat = 'identity', width = 0.4, fill="#006699", alpha = 0.8) +
    geom_line(df_model, mapping = aes(y=ajusted_lambda, x = year_sim), color = "#FF6666", alpha = 1, linewidth = 2, inherit.aes = FALSE)+
    geom_ribbon(df_model, mapping = aes(ymin = ajusted_lambda2.5, ymax = ajusted_lambda97.5, x = year_sim), fill = "#FF6666", alpha = 0.25, inherit.aes = FALSE)+
    geom_vline(aes(xintercept = 1995), color = "black",  alpha = 1, linewidth = 1, linetype = "dashed")+
    geom_vline(aes(xintercept = 2003), color = "black",  alpha = 1, linewidth = 1, linetype = "dashed")+
    geom_vline(aes(xintercept = 2006), color = "black",  alpha = 1, linewidth = 1, linetype = "dashed")+
    geom_vline(aes(xintercept = 2008), color = "black",  alpha = 1, linewidth = 1, linetype = "dashed")+
    geom_text(mapping=aes(x=1995.2, y=30, label="1"), size=6, vjust=-0.4, hjust=0)+
    geom_text(mapping=aes(x=2003.2, y=30, label="2"), size=6, vjust=-0.4, hjust=0)+
    geom_text(mapping=aes(x=2006.2, y=30, label="3"), size=6, vjust=-0.4, hjust=0)+
    geom_text(mapping=aes(x=2008.2, y=30, label="4"), size=6, vjust=-0.4, hjust=0)+
    theme_light()+ 
    scale_x_continuous(breaks=seq(1990, 2016, 5)) +
    scale_y_continuous(name = "Number of cases per person per year",
                       sec.axis = sec_axis( trans=~./S_1tc, name="Force of infection")) +
    #ggtitle("Figure 2. Estimated seroconversion rate validated on the ?? dataset") +
    xlab("Year") + ylab("Number of cases per person per year")+  
    labs(title = paste0(str_sub(Ab1, end=-2), "  ", Ab2))+
    theme(plot.title = element_text(size=22),
          plot.subtitle = element_text(size=16),
          axis.text=element_text(size=12),
          axis.title = element_text(size=16, face = "bold"), 
          plot.caption = element_text(hjust=0))
  
  ggsave(filename = paste0(str_sub(Ab1, end=-2), "_", Ab2, "_validation.jpeg"), plot_validation, width=10, height=8)
  
}







for (combinaison in 1:length(my_files_1Ab_1tc)){
  # Ab combination
  Ab1 <- strsplit(my_files_1Ab_1tc_tempo[[combinaison]], split = "MAL_")[[1]][2]
 
  
  # parameters vectors : med and IC95
  tempo_med <- c(median(list_results_1Ab_1tc[[combinaison]]$lambda), median(list_results_1Ab_1tc[[combinaison]]$tc), median(list_results_1Ab_1tc[[combinaison]]$delta))
  tempo_2.5 <- c(quantile(list_results_1Ab_1tc[[combinaison]]$lambda, c(0.025)), quantile(list_results_1Ab_1tc[[combinaison]]$tc, c(0.975)), quantile(list_results_1Ab_1tc[[combinaison]]$delta, c(0.025)))
  tempo_97.5 <- c(quantile(list_results_1Ab_1tc[[combinaison]]$lambda, c(0.975)), quantile(list_results_1Ab_1tc[[combinaison]]$tc, c(0.025)), quantile(list_results_1Ab_1tc[[combinaison]]$delta, c(0.975)))
  
  
  # find S to minimise the Standard error
  S<- c(seq(1, 300, by = 0.5))
  MSE <- c(rep(NA, 599))
  Sfinal <- as.data.frame(cbind(S, MSE))
  
  for(j in 1:length(S)){
    tempo <- Senegalese
    
    tempo$Error <- ifelse(tempo$age > tempo_med[2], abs(tempo$nb_cases - S[j]*tempo_med[1])^2, abs(tempo$nb_cases - S[j]*tempo_med[1]*tempo_med[3])^2)
    
    Sfinal[j,"MSE"] <- mean(tempo$Error)
  }
  S_1tc <- Sfinal$S[Sfinal$MSE== min(Sfinal$MSE)]
  
  # MSE computation
  MSE_tempo <- c(rep(NA, 550))
  for (k in 1:length(list_results_1Ab_1tc[[combinaison]]$lambda)){
    tempo <- Senegalese
    
    tempo$Error <- ifelse(tempo$age > list_results_1Ab_1tc[[combinaison]][k,"tc"], abs(tempo$nb_cases - S_1tc*list_results_1Ab_1tc[[combinaison]][k,"lambda"])^2, abs(tempo$nb_cases - S_1tc*list_results_1Ab_1tc[[combinaison]][k,"lambda"]*list_results_1Ab_1tc[[combinaison]][k,"delta"])^2)
    
    MSE_tempo[k] <- mean(tempo$Error) 
  }
  
  MSE_tempo <- as.data.frame(MSE_tempo)
  MSE_tempo$Ab1 <- Ab1
  MSE_tempo$Ab2 <- NA
  MSE_tempo$Ab3 <- NA
  MSE_tempo$nb_ab <- 1
  
  # save results
  write.csv(MSE_tempo, paste0("MSE_", Ab1, "_1tc_.csv"))
  
  
  # Plot
  year_sim<- c(seq(1990, 2016, by = 0.1))
  ajusted_lambda <- c(rep(NA, length(year_sim))) 
  ajusted_lambda2.5 <- c(rep(NA, length(year_sim))) 
  ajusted_lambda97.5 <- c(rep(NA, length(year_sim))) 
  
  df_model <- as.data.frame(cbind(year_sim, ajusted_lambda, ajusted_lambda2.5, ajusted_lambda97.5))
  df_model$age <- 2016 - df_model$year_sim
  
  df_model$ajusted_lambda <- ifelse(df_model$age > tempo_med[2],
                                    S_1tc*tempo_med[1],
                                    S_1tc*tempo_med[1]*tempo_med[3])
  
  df_model$ajusted_lambda2.5 <- ifelse(df_model$age > tempo_2.5[2],
                                       S_1tc*tempo_2.5[1],
                                       S_1tc*tempo_2.5[1]*tempo_2.5[3])
  
  df_model$ajusted_lambda97.5 <- ifelse(df_model$age > tempo_97.5[2],
                                        S_1tc*tempo_97.5[1],
                                        S_1tc*tempo_97.5[1]*tempo_97.5[3])
  
  
  plot_validation <- ggplot(Senegalese, mapping = aes(y = nb_cases, x = year)) + geom_bar(stat = 'identity', width = 0.4, fill="#006699", alpha = 0.8) +
    geom_line(df_model, mapping = aes(y=ajusted_lambda, x = year_sim), color = "#FFCCCC", alpha = 1, linewidth = 2, inherit.aes = FALSE)+
    geom_ribbon(df_model, mapping = aes(ymin = ajusted_lambda2.5, ymax = ajusted_lambda97.5, x = year_sim), fill = "#FFCCCC", alpha = 0.25, inherit.aes = FALSE)+
    geom_vline(aes(xintercept = 1995), color = "black",  alpha = 1, linewidth = 1, linetype = "dashed")+
    geom_vline(aes(xintercept = 2003), color = "black",  alpha = 1, linewidth = 1, linetype = "dashed")+
    geom_vline(aes(xintercept = 2006), color = "black",  alpha = 1, linewidth = 1, linetype = "dashed")+
    geom_vline(aes(xintercept = 2008), color = "black",  alpha = 1, linewidth = 1, linetype = "dashed")+
    geom_text(mapping=aes(x=1995.2, y=30, label="1"), size=6, vjust=-0.4, hjust=0)+
    geom_text(mapping=aes(x=2003.2, y=30, label="2"), size=6, vjust=-0.4, hjust=0)+
    geom_text(mapping=aes(x=2006.2, y=30, label="3"), size=6, vjust=-0.4, hjust=0)+
    geom_text(mapping=aes(x=2008.2, y=30, label="4"), size=6, vjust=-0.4, hjust=0)+
    theme_light()+ 
    scale_x_continuous(breaks=seq(1990, 2016, 5)) +
    scale_y_continuous(name = "Number of cases per person per year",
                       sec.axis = sec_axis( trans=~./S_1tc, name="Force of infection")) +
    #ggtitle("Figure 2. Estimated seroconversion rate validated on the ?? dataset") +
    xlab("Year") + ylab("Number of cases per person per year")+  
    labs(title = paste0(Ab1, end=-2))+
    theme(plot.title = element_text(size=22),
          plot.subtitle = element_text(size=16),
          axis.text=element_text(size=12),
          axis.title = element_text(size=16, face = "bold"), 
          plot.caption = element_text(hjust=0))
  
  ggsave(filename = paste0(Ab1, "_validation_1tc.jpeg"), plot_validation, width=10, height=8)
}





# Download packages
rm(list = ls())
library(base)
library(ggplot2)
library(ggpubr)
library(stringr)



# Download results from modelisation
my_files_3Ab_tempo = list.files(path = "./3Ab_1tc", pattern = "*_dielmo.csv")
my_files_3Ab <- c(paste0("./3Ab_1tc/", my_files_3Ab_tempo))
list_results_3Ab <- lapply(my_files_3Ab, read.csv)
my_files_3Ab_tempo <- str_sub(my_files_3Ab_tempo, end=-16)

vec_title <- c("PfMSP1", "PfAMA1", "PfEtramp4", "PfGlurpR2", "PfMSP2-Dd2", "PfSEA1", "PfHSP40", "PfSBP1", "PfMSP2-CH150", "PfCSP")
vec_malaria <- c("MAL_PfMSP1", "MAL_PfAMA1", "MAL_PfEtramp4_Ag2", "MAL_PfGlurpR2", "MAL_PfMSP2_Dd2", "MAL_PfSEA1", "MAL_PfHSP40_Ag1", "MAL_PfSBP1", "MAL_PfMSP2_CH150", "MAL_PfCSP")

for (combinaison in 1:length(my_files_3Ab_tempo)){
  Ab1_tempo <- str_sub(strsplit(my_files_3Ab_tempo[[combinaison]], split = "MAL_")[[1]][2], end=-2)
  Ab2_tempo <- str_sub(strsplit(my_files_3Ab_tempo[[combinaison]], split = "MAL_")[[1]][3], end=-2)
  Ab3_tempo <- strsplit(my_files_3Ab_tempo[[combinaison]], split = "MAL_")[[1]][4]
  
  Ab1 <- vec_title[grep(Ab1_tempo, vec_malaria)]
  Ab2 <- vec_title[grep(Ab2_tempo, vec_malaria)]
  Ab3 <- vec_title[grep(Ab3_tempo, vec_malaria)]
  
  
  
  # Download results files of all models including combinaison of those 3 antigens
  my_files_1Ab = c(list.files(path = "./1Ab_1tc", pattern = Ab1_tempo), list.files(path = "./1Ab_1tc", pattern = Ab2_tempo), list.files(path = "./1Ab_1tc", pattern = Ab3_tempo))
  my_files_1Ab <- c(paste0("./1Ab_1tc/", my_files_1Ab))
  list_results_1Ab <- lapply(my_files_1Ab, read.csv)
  
  my_files_2Ab = c(intersect(list.files(path = "./2Ab_1tc", pattern = Ab1_tempo), list.files(path = "./2Ab_1tc", pattern = Ab2_tempo)),
                   intersect(list.files(path = "./2Ab_1tc", pattern = Ab1_tempo), list.files(path = "./2Ab_1tc", pattern = Ab3_tempo)),
                   intersect(list.files(path = "./2Ab_1tc", pattern = Ab2_tempo), list.files(path = "./2Ab_1tc", pattern = Ab3_tempo)))
  my_files_2Ab <- c(paste0("./2Ab_1tc/", my_files_2Ab))
  list_results_2Ab <- lapply(my_files_2Ab, read.csv)
  
  df_plot <- as.data.frame(matrix(NA, ncol = 6, nrow = 21))
  colnames(df_plot) <- c("Model", "Nb_ag", "Parameter", "mean", "IClow", "IChigh")
  
  
  df_plot$Model <- rep(c(paste0(Ab1),
                     paste0(Ab2),
                     paste0(Ab3),
                     paste0(Ab1, " ", Ab2),
                     paste0(Ab1, " ", Ab3),
                     paste0(Ab2, " ", Ab3),
                     paste0(Ab1, " ", Ab2, " ", Ab3)), each = 3)
  
  df_plot$Parameter <- c(rep(c("lambda", "tc", "delta"), 6),c("lambda", "tc", "delta"))
  df_plot$Nb_ag <- c(rep(1, 9), rep(2, 9), rep(3, 3)) 
  
  
  missing_results1 <- c(length(list.files(path = "./1Ab_1tc", pattern = Ab1_tempo)) == 0,
                               length(list.files(path = "./1Ab_1tc", pattern = Ab2_tempo)) == 0,
                                      length(list.files(path = "./1Ab_1tc", pattern = Ab3_tempo)) == 0)
  present_results1 <- which(missing_results1 == FALSE)
  
  if(length(list.files(path = "./1Ab_1tc", pattern = Ab1_tempo)) != 0){
    df_plot$mean[df_plot$Model == Ab1] <- c(mean(list_results_1Ab[[1]]$lambda), mean(list_results_1Ab[[1]]$tc),mean(list_results_1Ab[[1]]$delta))
    df_plot$IClow[df_plot$Model == Ab1] <- c(quantile(list_results_1Ab[[1]]$lambda, c(0.025)), quantile(list_results_1Ab[[1]]$tc, c(0.025)), quantile(list_results_1Ab[[1]]$delta, c(0.025)))
    df_plot$IChigh[df_plot$Model == Ab1] <- c(quantile(list_results_1Ab[[1]]$lambda, c(0.975)), quantile(list_results_1Ab[[1]]$tc, c(0.975)), quantile(list_results_1Ab[[1]]$delta, c(0.975)))
  }
         
  if(length(list.files(path = "./1Ab_1tc", pattern = Ab2_tempo)) != 0){
    df_plot$mean[df_plot$Model == Ab2] <- c(mean(list_results_1Ab[[which(present_results1 ==2)]]$lambda), mean(list_results_1Ab[[which(present_results1 ==2)]]$tc),mean(list_results_1Ab[[which(present_results1 ==2)]]$delta))
    df_plot$IClow[df_plot$Model == Ab2] <- c(quantile(list_results_1Ab[[which(present_results1 ==2)]]$lambda, c(0.025)), quantile(list_results_1Ab[[which(present_results1 ==2)]]$tc, c(0.025)), quantile(list_results_1Ab[[which(present_results1 ==2)]]$delta, c(0.025)))
    df_plot$IChigh[df_plot$Model == Ab2] <- c(quantile(list_results_1Ab[[which(present_results1 ==2)]]$lambda, c(0.975)), quantile(list_results_1Ab[[which(present_results1 ==2)]]$tc, c(0.975)), quantile(list_results_1Ab[[which(present_results1 ==2)]]$delta, c(0.975)))
  }
  
  if(length(list.files(path = "./1Ab_1tc", pattern = Ab3_tempo)) != 0){
    df_plot$mean[df_plot$Model == Ab3] <- c(mean(list_results_1Ab[[which(present_results1 ==3)]]$lambda), mean(list_results_1Ab[[which(present_results1 ==3)]]$tc),mean(list_results_1Ab[[which(present_results1 ==3)]]$delta))
    df_plot$IClow[df_plot$Model == Ab3] <- c(quantile(list_results_1Ab[[which(present_results1 ==3)]]$lambda, c(0.025)), quantile(list_results_1Ab[[which(present_results1 ==3)]]$tc, c(0.025)), quantile(list_results_1Ab[[which(present_results1 ==3)]]$delta, c(0.025)))
    df_plot$IChigh[df_plot$Model == Ab3] <- c(quantile(list_results_1Ab[[which(present_results1 ==3)]]$lambda, c(0.975)), quantile(list_results_1Ab[[which(present_results1 ==3)]]$tc, c(0.975)), quantile(list_results_1Ab[[which(present_results1 ==3)]]$delta, c(0.975)))
  }
  
  
  missing_results <- c(length(intersect(list.files(path = "./2Ab_1tc", pattern = Ab1_tempo), list.files(path = "./2Ab_1tc", pattern = Ab2_tempo))) == 0,
                       length(intersect(list.files(path = "./2Ab_1tc", pattern = Ab1_tempo), list.files(path = "./2Ab_1tc", pattern = Ab3_tempo))) == 0,
                       length(intersect(list.files(path = "./2Ab_1tc", pattern = Ab2_tempo), list.files(path = "./2Ab_1tc", pattern = Ab3_tempo))) == 0)
  present_results <- which(missing_results == FALSE)
  
  if(length(intersect(list.files(path = "./2Ab_1tc", pattern = Ab1_tempo), list.files(path = "./2Ab_1tc", pattern = Ab2_tempo))) != 0){
    df_plot$mean[df_plot$Model == paste0(Ab1, " ", Ab2)] <- c(mean(list_results_2Ab[[1]]$lambda), mean(list_results_2Ab[[1]]$time_c),mean(list_results_2Ab[[1]]$delta))
    df_plot$IClow[df_plot$Model == paste0(Ab1, " ", Ab2)] <- c(quantile(list_results_2Ab[[1]]$lambda, c(0.025)), quantile(list_results_2Ab[[1]]$time_c, c(0.025)), quantile(list_results_2Ab[[1]]$delta, c(0.025)))
    df_plot$IChigh[df_plot$Model == paste0(Ab1, " ", Ab2)] <- c(quantile(list_results_2Ab[[1]]$lambda, c(0.975)), quantile(list_results_2Ab[[1]]$time_c, c(0.975)), quantile(list_results_2Ab[[1]]$delta, c(0.975)))
  }
  
  if(length(intersect(list.files(path = "./2Ab_1tc", pattern = Ab1_tempo), list.files(path = "./2Ab_1tc", pattern = Ab3_tempo))) != 0){
    df_plot$mean[df_plot$Model == paste0(Ab1, " ", Ab3)] <- c(mean(list_results_2Ab[[which(present_results ==2)]]$lambda), mean(list_results_2Ab[[which(present_results ==2)]]$time_c),mean(list_results_2Ab[[which(present_results ==2)]]$delta))
    df_plot$IClow[df_plot$Model == paste0(Ab1, " ", Ab3)] <- c(quantile(list_results_2Ab[[which(present_results ==2)]]$lambda, c(0.025)), quantile(list_results_2Ab[[which(present_results ==2)]]$time_c, c(0.025)), quantile(list_results_2Ab[[which(present_results ==2)]]$delta, c(0.025)))
    df_plot$IChigh[df_plot$Model == paste0(Ab1, " ", Ab3)] <- c(quantile(list_results_2Ab[[which(present_results ==2)]]$lambda, c(0.975)), quantile(list_results_2Ab[[which(present_results ==2)]]$time_c, c(0.975)), quantile(list_results_2Ab[[which(present_results ==2)]]$delta, c(0.975)))
  }
  
  if(length(intersect(list.files(path = "./2Ab_1tc", pattern = Ab2_tempo), list.files(path = "./2Ab_1tc", pattern = Ab3_tempo))) != 0){
    df_plot$mean[df_plot$Model == paste0(Ab2, " ", Ab3)] <- c(mean(list_results_2Ab[[which(present_results ==3)]]$lambda), mean(list_results_2Ab[[which(present_results ==3)]]$time_c),mean(list_results_2Ab[[which(present_results ==3)]]$delta))
    df_plot$IClow[df_plot$Model == paste0(Ab2, " ", Ab3)] <- c(quantile(list_results_2Ab[[which(present_results ==3)]]$lambda, c(0.025)), quantile(list_results_2Ab[[which(present_results ==3)]]$time_c, c(0.025)), quantile(list_results_2Ab[[which(present_results ==3)]]$delta, c(0.025)))
    df_plot$IChigh[df_plot$Model == paste0(Ab2, " ", Ab3)] <- c(quantile(list_results_2Ab[[which(present_results ==3)]]$lambda, c(0.975)), quantile(list_results_2Ab[[which(present_results ==3)]]$time_c, c(0.975)), quantile(list_results_2Ab[[which(present_results ==3)]]$delta, c(0.975)))
  }
  
  
  df_plot$mean[df_plot$Model == paste0(Ab1, " ", Ab2, " ", Ab3)] <- c(mean(list_results_3Ab[[combinaison]]$lambda), mean(list_results_3Ab[[combinaison]]$time_c),mean(list_results_3Ab[[combinaison]]$delta))
  df_plot$IClow[df_plot$Model == paste0(Ab1, " ", Ab2, " ", Ab3)] <- c(quantile(list_results_3Ab[[combinaison]]$lambda, c(0.025)), quantile(list_results_3Ab[[combinaison]]$time_c, c(0.025)), quantile(list_results_3Ab[[combinaison]]$delta, c(0.025)))
  df_plot$IChigh[df_plot$Model == paste0(Ab1, " ", Ab2, " ", Ab3)] <- c(quantile(list_results_3Ab[[combinaison]]$lambda, c(0.975)), quantile(list_results_3Ab[[combinaison]]$time_c, c(0.975)), quantile(list_results_3Ab[[combinaison]]$delta, c(0.975)))
  
  
  
  plot_comparaison_postdistrib <- ggplot(df_plot, aes(x = mean, y = Model)) + geom_point(aes(color = as.character(Nb_ag)), size = 3)+
    geom_pointrange(aes(xmin = as.numeric(IClow), xmax = as.numeric(IChigh), color = as.character(Nb_ag)), linewidth = 1.5)+
    scale_color_manual("Number of antigens", values=c("#FF9999", "#FF6666", "#990000"))+
    facet_grid(Nb_ag ~ Parameter, scales = "free")+
    theme_light()+
    ggtitle("Posterior distribution comparison") +
    xlab("") + ylab("")+
    theme(plot.title = element_text(size=22),
          axis.line.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.y=element_blank(),
          axis.title.x=element_blank(),
          axis.text=element_text(size=12),
          axis.title=element_text(size=14),
          strip.text.x = element_text(size = 16),
          strip.text.y = element_text(size = 16))

  ggsave(filename =paste0("Supplementary_", combinaison, "_postdistrib.jpeg"), plot_comparaison_postdistrib, width=12, height=8)

  
}




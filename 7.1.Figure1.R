rm(list=ls())

library("bayesplot")
library("RColorBrewer")
library("dplyr")
library("binom")
library(ggplot2)
library(patchwork)
library(ggpubr)



#----------#
#-- Data --#
#----------#
# Download Data: Samples and negative controls
dielmo_Ndiop <- read.csv("0.Dielmo_Ndiop_PFalci.csv")
vec_malaria <- names(dielmo_Ndiop)[grepl("MAL_Pf", names(dielmo_Ndiop))]
vec_epidemio <- colnames(dielmo_Ndiop)[1:7]

neg <- read.csv("0.Seroped_nMFI.csv")





cutoffs <- read.csv("cutoff_df.csv")




dielmo_Ndiop$population <- NA
dielmo_Ndiop$population[dielmo_Ndiop$village == "Ndiop" & dielmo_Ndiop$year == "2016"] <- "Ndiop 2016"
dielmo_Ndiop$population[dielmo_Ndiop$village == "Ndiop" & dielmo_Ndiop$year == "2018"] <- "Ndiop 2018"
dielmo_Ndiop$population[dielmo_Ndiop$village == "Dielmo" & dielmo_Ndiop$year == "2016"] <- "Dielmo 2016"
dielmo_Ndiop$population[dielmo_Ndiop$village == "Dielmo" & dielmo_Ndiop$year == "2018"] <- "Dielmo 2018"

dielmo_Ndiop <- dielmo_Ndiop %>%
  mutate(population = factor(population, levels = c("Negatives", "Ndiop 2016", "Ndiop 2018", "Dielmo 2016", "Dielmo 2018"))) %>%  # Optional, for legend order
  arrange(population)

vec_title <- c("PfMSP1", "PfAMA1", "PfEtramp4", "PfGlurpR2", "PfMSP2-Dd2", "PfSEA1", "PfHSP40", "PfSBP1", "PfMSP2-CH150", "PfCSP")


plot_AgeVSRAU <- vector('list')
for (ag in 1: length(vec_malaria)){
  tempo <- dielmo_Ndiop[c(vec_malaria[ag], "population", "year", "age")]
  colnames(tempo) <- c("Ab", "Population", "Year", "Age")
  tempo <- tempo[!is.na(tempo$Ab),]
  tempo$Ab <- ifelse(tempo$Ab !=0, tempo$Ab, 0.001)
  
  
  tempo_neg <- as.data.frame(neg[c(vec_malaria_neg_nMFI[ag], "id_sample")])
  colnames(tempo_neg) <- c("Ab", "id")
  tempo_neg <- tempo_neg[!is.na(tempo_neg$Ab),]
  tempo_neg$Population <- "Negatives"
  tempo_neg$Population <- factor(tempo_neg$Population, levels = c("Negatives", "Ndiop 2016", "Ndiop 2018", "Dielmo 2016", "Dielmo 2018"))
  
  plot_AgeVSRAU[[ag]] <-local({
  temp_value <- exp(cutoffs$cutoff[ag])


  p1 <- ggplot(data = tempo, aes(x = Age, y = Ab, by = Population)) + geom_point(aes(color = Population), size = 2.5, position="jitter", alpha = c(0.5), show.legend = TRUE)+
    scale_y_continuous(trans="log", limits = c(0.0005, 16), breaks = c(0, 0.01, 0.1,1, 10))+
    scale_fill_manual(name = 'Population', values = c("Negatives" = "#66CC66", "Ndiop 2016" = "#003366", "Ndiop 2018" = "#0099CC", "Dielmo 2016" = "#990000", "Dielmo 2018" = "#FF9999"), drop = FALSE) +
    scale_color_manual(name = 'Population', values = c("Negatives" = "#66CC66", "Ndiop 2016" = "#003366", "Ndiop 2018" = "#0099CC", "Dielmo 2016" = "#990000", "Dielmo 2018" = "#FF9999"), drop = FALSE) +
    geom_hline(aes(yintercept = temp_value), linewidth = 1.5, linetype = "dashed")+
    # scale_fill_manual(values = c("#990000", "#FF9999", "#0099CC")) +
    # scale_color_manual(values = c("#990000", "#FF9999", "#0099CC")) +
    xlab("Age")+ylab("nMFI")+theme_light()+
    theme(plot.title = element_text(size=20, face = "bold"),
          axis.text=element_text(size=18),
          axis.title = element_text(size=20, face = "bold"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          legend.text = element_text(size=20),
          legend.ticks = element_line(linewidth = 16),
          legend.title = element_text(size=22))+
    guides(color = guide_legend(override.aes = list(size = 12)))

  
  p2 <- ggplot(tempo_neg, aes(x = 0, y = Ab)) +
    geom_point(aes(color = Population), size = 2.5, position="jitter", alpha = 0.5, show.legend = FALSE) +
    geom_hline(aes(yintercept = temp_value), linewidth = 1.5, linetype = "dashed")+
    scale_fill_manual(name = 'Population', values = c("Negatives" = "#66CC66", "Ndiop 2016" = "#003366", "Ndiop 2018" = "#0099CC", "Dielmo 2016" = "#990000", "Dielmo 2018" = "#FF9999"), drop = FALSE) +
    scale_color_manual(name = 'Population', values = c("Negatives" = "#66CC66", "Ndiop 2016" = "#003366", "Ndiop 2018" = "#0099CC", "Dielmo 2016" = "#990000", "Dielmo 2018" = "#FF9999"), drop = FALSE) +
    scale_y_continuous(trans="log", limits = c(0.0005, 16), breaks = c(0,0.01, 0.1,1, 10))+
    scale_x_continuous(limits = c(-1, 1), breaks = 0, labels = "Negatives") +
    ylab("nMFI")+theme_light()+
    labs(x = NULL)+
    theme(plot.title = element_text(size=20, face = "bold"),
          axis.text=element_text(size=18),
          axis.title = element_text(size=20, face = "bold"),
          legend.text = element_text(size=20),
          legend.ticks = element_line(linewidth = 16),
          legend.title = element_text(size=22))+
    guides(color = guide_legend(override.aes = list(size = 12)))

  
 p2 + p1 + plot_layout(ncol = 2, widths = c(1, 4), guides = "collect")+ 
    plot_annotation(title = paste0(vec_title[ag]))& 
    theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.11, vjust = 0))
})
}

plot_distrib_CHik<-ggarrange(plot_AgeVSRAU[[1]], plot_AgeVSRAU[[2]], plot_AgeVSRAU[[3]],plot_AgeVSRAU[[4]],
                             plot_AgeVSRAU[[5]], plot_AgeVSRAU[[6]], plot_AgeVSRAU[[7]], plot_AgeVSRAU[[8]], 
                             plot_AgeVSRAU[[9]], plot_AgeVSRAU[[10]], common.legend = TRUE, legend="bottom") 

plot_distrib_CHik
ggsave(filename = paste0("Figure1.jpeg"), plot_distrib_CHik, width=25, height=15)


rm(list = ls())
library(ggplot2)
library(ggpubr)
library(base)
library("deSolve")
library(binom)



#----------#
#-- Data --#
#----------#
# Download Data: Samples and negative controls
dielmo_Ndiop <- read.csv("0.Dielmo_Ndiop_PFalci.csv")
vec_malaria <- names(dielmo_Ndiop)[grepl("MAL_Pf", names(dielmo_Ndiop))]
vec_epidemio <- colnames(dielmo_Ndiop)[1:7]

dielmo <- dielmo_Ndiop[dielmo_Ndiop$village == "Ndiop",]




# Binarization of data
cutoffs_dielmo <- read.csv("cutoff_df.csv")
for (i in 1:length(vec_malaria)){
  dielmo_Ndiop[,vec_malaria[i]] <- ifelse(log(dielmo_Ndiop[, vec_malaria[i]]) >= cutoffs_dielmo$cutoff[cutoffs_dielmo$Ab == vec_malaria[i]], 1, 0)
}


# Split in two cohorts
dielmo_2016 <- dielmo_Ndiop[dielmo_Ndiop$year == "2016",]
dielmo_2018 <- dielmo_Ndiop[dielmo_Ndiop$year == "2018",]



# download results
results_1AB_1tc = list.files(pattern = "*ndiop.csv")
list_results_1AB_1tc <- lapply(results_1AB_1tc, read.csv)



for (Ab in 1:length(vec_malaria)){
dfanalyses <- list_results_1AB_1tc[[which(grepl(vec_malaria[Ab], results_1AB_1tc))]]
dfanalyses$chain <- as.character(c(rep(c(1,2,3,4), each = 3000)))
dfanalyses$iteration <- c(rep(seq(1,3000, 1), 4))

# Plots
############################
## POSTERIOR DISTRIBUTION ##
############################
plot_chains1 <- list()
for (param in 1: (dim(dfanalyses)[2] - 2)){
  param_tempo <- colnames(dfanalyses)[param]
  df_tempo <- dfanalyses[,c(param_tempo, "chain")]
  colnames(df_tempo) <- c("parametre", "chain")
  
  plot_chains1[[param]]<- ggplot(df_tempo) + geom_density(aes(x = parametre, group = chain, color = chain), size = 1.5, alpha = 0.2) +
  scale_color_manual(values = c("#0033FF", "#0066FF", "#0099FF", "#0099CC"))+
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
    scale_color_manual(values = c("#0033FF", "#0066FF", "#0099FF", "#0099CC"))+
    theme_classic()+ 
    theme(plot.title = element_text(size=22, face = "bold"),
          axis.text=element_text(size=12),
          axis.title=element_text(size=16, face = "bold"),
          axis.title.x=element_blank())
  
}

plot_chains <- ggarrange(plot_chains1[[1]], plot_chains1[[2]], plot_chains1[[3]], plot_chains1[[4]], 
                  plot_chains2[[1]], plot_chains2[[2]], plot_chains2[[3]], plot_chains2[[4]], 
                  ncol = 4, nrow = 2, common.legend = TRUE)

 ggsave(filename = paste0(vec_malaria[Ab], "_1tc_chains_ndiop_Priors.jpeg"), plot_chains, width=12, height=8)






N_par <- 4   # number of parameters
N_mcmc <- dim(dfanalyses)[1]



##############################
## ESTIMATED SEROPREVALENCE ##
##############################
# ODE systems
equa_noChange <- function(Time, State, Pars){
  with(as.list(c(State,Pars)), {
    SPA <- State[1]
    
    dSPA <- lambda*(1-SPA) - rhoA*SPA
    
    return(list(c(dSPA)))
  })
}

equa_Change1 <- function(Time, State, Pars){
  with(as.list(c(State,Pars)), {
    SPA <- State[1]
    
    dSPA <- delta*lambda*(1-SPA) - rhoA*SPA
    
    return(list(c(dSPA)))
  })
}


# Initial probabilities
yini <- c(SPA = 0)




# Median estimated seroprevalence
age_seq <- seq(from=0, to=95, by=0.1)                                            # age range 
par_median <- apply(X=dfanalyses[,1:N_par], MARGIN=2, FUN=median)                # median of each parameters posterior distribution
tc <- par_median[4]

age_change_2016 <- age_seq[age_seq <= tc]                                        # age before drop of transmission
age_nochange_2016 <- age_seq[age_seq> tc] - tc                                   # age after drop of transmission


Est_median_2016 <- as.data.frame(ode(yini, age_change_2016, equa_Change1, par_median))   
Est_median_2016_bis <- as.data.frame(ode(as.numeric(Est_median_2016[length(age_change_2016),2]), age_nochange_2016, equa_noChange, par_median))
colnames(Est_median_2016_bis) <- colnames(Est_median_2016)
Est_median_2016 <- rbind(Est_median_2016, Est_median_2016_bis)

for (i in 1:length(Est_median_2016$time)){
  Est_median_2016[i,"A"] <- Est_median_2016[i,2]
}


age_change_2018 <- age_seq[age_seq <= (tc+2)]                                        
age_nochange_2018 <- age_seq[age_seq> (tc+2)] - (tc+2)                                   


Est_median_2018 <- as.data.frame(ode(yini, age_change_2018, equa_Change1, par_median))   
Est_median_2018_bis <- as.data.frame(ode(as.numeric(Est_median_2018[length(age_change_2018),2]), age_nochange_2018, equa_noChange, par_median))
colnames(Est_median_2018_bis) <- colnames(Est_median_2018)
Est_median_2018 <- rbind(Est_median_2018, Est_median_2018_bis)

for (i in 1:length(Est_median_2018$time)){
  Est_median_2018[i,"A"] <- Est_median_2018[i,2]
}





# Confidence Intervalle
N_sam = 700                                                         # Number of samples in the posterior distribution
sam_seq = round(seq(from=1, to=nrow(dfanalyses), length=N_sam))     # Sample N_sam in posterior distribution

# 2016 cohort____________________________________________________
Est_IC_2016 = as.data.frame(matrix(NA, nrow=N_sam, ncol=length(age_seq)))
for(k in 1:N_sam)
{
  param_tempo <- dfanalyses[sam_seq[k],1:N_par]

  tc <- as.numeric(param_tempo[4])
  age_change1 <- age_seq[age_seq <= tc]
  age_nochange <- c(0.000001, age_seq[age_seq> tc] - tc, 100-tc)
  
  
  tempo2 <- as.data.frame(ode(yini, age_change1, equa_Change1, param_tempo))   
  tempo <- as.data.frame(ode(as.numeric(tempo2[length(age_change1),2]), age_nochange, equa_noChange, param_tempo))
  colnames(tempo) <- colnames(tempo2)
  tempo <- rbind(tempo2, tempo[-c(1,length(tempo$time)),])
  
  
  for (l in 1:length(tempo$time)){
      tempo[l,"A"] <- tempo[l,2]
  }


  Est_IC_2016[k,] <- t(tempo$A)
}


IC_2016 = matrix(NA, nrow=3, ncol=length(age_seq))
for(j in 1:length(age_seq))
{
  IC_2016[,j] = quantile( Est_IC_2016[,j], prob=c(0.025, 0.5, 0.975), na.rm = TRUE )
}


Est_2016<- as.data.frame(t(IC_2016))
Est_2016$med <- Est_median_2016$A
Est_2016$time <- age_seq
Est_2016$cohort <- "2016"





# 2018 cohort____________________________________________________
Est_IC_2018 = as.data.frame(matrix(NA, nrow=N_sam, ncol=length(age_seq)))
for(k in 1:N_sam)
{
  param_tempo <- dfanalyses[sam_seq[k],1:N_par]
  
  tc <- as.numeric(param_tempo[4])
  age_change1 <- age_seq[age_seq <= (tc+2)]
  age_nochange <- c(0.000001, age_seq[age_seq> (tc+2)] - (tc+2), 100-(tc+2))
  
  
  tempo2 <- as.data.frame(ode(yini, age_change1, equa_Change1, param_tempo))   
  tempo <- as.data.frame(ode(as.numeric(tempo2[length(age_change1),2]), age_nochange, equa_noChange, param_tempo))
  colnames(tempo) <- colnames(tempo2)
  tempo <- rbind(tempo2, tempo[-c(1,length(tempo$time)),])
  
  
  for (l in 1:length(tempo$time)){
    tempo[l,"A"] <- tempo[l,2]
  }
  
  
  Est_IC_2018[k,] <- t(tempo$A)
}


IC_2018 = matrix(NA, nrow=3, ncol=length(age_seq))
for(j in 1:length(age_seq))
{
  IC_2018[,j] = quantile( Est_IC_2018[,j], prob=c(0.025, 0.5, 0.975), na.rm = TRUE )
}


Est_2018<- as.data.frame(t(IC_2018))
Est_2018$med <- Est_median_2018$A
Est_2018$time <- age_seq
Est_2018$cohort <- "2018"


Est_final <- rbind(Est_2016, Est_2018)
Est_final$Data <- "Estimated" 

#write.csv(Est_final,paste0("Est_seroprev_", vec_malaria[Ab], "_1tc.csv"))



##############
## PLOTTING ##
##############
# Prepare obs data for plotting 
age_bins     <- seq(from=0, to=95, by=5)
age_bins_mid <- seq(from=2.5, to=92.5, by=5) 
N_bins <- length(age_bins) - 1 


SP_bins2016 <- matrix(NA, nrow=N_bins, ncol=3)
colnames(SP_bins2016) <- c("med", "low_95", "high_95")

for(i in 1:N_bins)
{
  index <- which( dielmo_2016$age>age_bins[i] & dielmo_2016$age<=age_bins[i+1] ) 
  temp_AB  <- dielmo_2016[index, vec_malaria[Ab]]                                           
  
  SP_bins2016[i,] <- as.numeric(as.vector(
    binom.confint( sum(temp_AB), length(temp_AB), method="wilson")[1,4:6]
  ))
}

SP_bins2016 <- as.data.frame(SP_bins2016)
SP_bins2016$age <- age_bins_mid
SP_bins2016$cohort <- "2016"



SP_bins2018 <- matrix(NA, nrow=N_bins, ncol=3)
colnames(SP_bins2018) <- c("med", "low_95", "high_95")

for(i in 1:N_bins)
{
  index <- which( dielmo_2018$age>age_bins[i] & dielmo_2018$age<=age_bins[i+1] ) 
  temp_AB  <- dielmo_2018[index, vec_malaria[Ab]]                                           
  
  SP_bins2018[i,] <- as.numeric(as.vector(
    binom.confint( sum(temp_AB), length(temp_AB), method="wilson")[1,4:6]
  ))
}

SP_bins2018 <- as.data.frame(SP_bins2018)
SP_bins2018$age <- age_bins_mid
SP_bins2018$cohort <- "2018"


SP_bins <- rbind(SP_bins2016, SP_bins2018)



plot_seroprev_vs_age <- ggplot(SP_bins, aes(x = age, y = med, group = cohort)) +
  geom_point(data = SP_bins, aes(color = cohort), size = 2.5, shape = 19, position=position_dodge(1))+
  geom_errorbar(data = SP_bins, aes(ymin=low_95, ymax=high_95,  group = cohort, color = cohort), width=0, alpha = 0.5,
                position=position_dodge(1)) +
  geom_ribbon(data = Est_final, aes(x = time, ymin = V1, ymax = V3, group = cohort, color = cohort, fill = cohort), alpha = 0.2)+
  geom_line(data = Est_final, aes(x = time, y = med,  group = cohort, color = cohort), alpha = 1, linewidth = 2)+
  scale_fill_manual(values = c("#990000", "#FF9999")) +
  scale_color_manual(values = c("#990000", "#FF9999")) +
  theme_classic()+
  ggtitle(paste0(vec_malaria[Ab])) +
  xlab("Age (Years)") + ylab("Seroprevalence")+
  theme(plot.title = element_text(size=22, face = "bold"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=16, face = "bold"))+
  ylim(0,1)
plot_seroprev_vs_age


ggsave(filename = paste0(vec_malaria[Ab], "_1tc_VPC_ndiop_Priors.jpeg"), plot_seroprev_vs_age, width=12, height=8)



}











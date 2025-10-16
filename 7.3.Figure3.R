rm(list = ls())
library(ggplot2)
library(ggpubr)
library(base)
library(forcats)
library(stringr)
library("deSolve")
library(binom)



#######################
### PLOT MSE RANKED ###
#######################
# Download results from modelisation
my_files = list.files(pattern = "MSE_Pf*")
list_results <- lapply(my_files, read.csv)
my_files_tempo <- str_sub(my_files, start =5, end=-6)



# Create global database
df_final <- as.data.frame(matrix(ncol = 6, nrow = 26))
colnames(df_final) <- c("Combinaison", "Model", "mean", "min", "max", "nb_ab")

for (combinaison in 1:length(my_files)){
  
  df_final[combinaison,] <- c(str_sub(my_files_tempo[combinaison], end=-5),
                              str_sub(my_files_tempo[combinaison], start = -3),
                              mean(list_results[[combinaison]]$MSE_tempo), 
                              quantile(list_results[[combinaison]]$MSE_tempo, c(0.025)), 
                              quantile(list_results[[combinaison]]$MSE_tempo, c(0.975)), 
                              list_results[[combinaison]]$nb_ab)
  
  
}

df_final$Combinaison <- gsub("_", " ", df_final$Combinaison)
df_final$Combinaison <- gsub("Ag1", "", df_final$Combinaison)
df_final$Combinaison <- gsub("Ag2", "", df_final$Combinaison)
df_final$Combinaison <- gsub("PfMSP2 Dd2", "PfMSP2-Dd2", df_final$Combinaison)
df_final <- df_final[order(as.numeric(df_final$mean), decreasing = FALSE),]
df_final$Combinaison <- factor(df_final$Combinaison, levels=rev(c(unique(df_final$Combinaison))))

#######################################################################################################
##  Kruskal-Wallis rank sum test  
kruskal.test(as.numeric(mean)~nb_ab, data=df_final)
#######################################################################################################

plot_final <- ggplot(data = df_final, aes(y=Combinaison, x = as.numeric(mean), color = nb_ab)) + geom_point(size = 5, position = position_dodge(width = 0.2))+ 
  geom_pointrange(aes(xmin = as.numeric(min), xmax = as.numeric(max)), position = position_dodge(width = 0.2), linewidth = 2) +
  scale_x_log10(breaks = c(10, 50, 100, 200), position = "top")+
  theme_light()+
  ggtitle("c. Mean Squared Error") +
  xlab("Mean Squared Error") + ylab("Antibody")+
  theme(plot.title = element_text(size=27, face = "bold"),
        axis.line.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text=element_text(size=20),
        axis.title=element_text(size=14),
        legend.text = element_text(size=24),
        legend.ticks = element_line(linewidth = 22),
        legend.title = element_text(size=24),
        legend.position="bottom",
        legend.key.size = unit(1.5,"line"))+
  scale_color_manual("Number of antigens", values=c("#FFCCCC", "#FF6666", "#990000"))
plot_final






#######################
### PLOT VALIDATION ###
#######################
my_files_2Ab_tempo = list.files(pattern = "MAL_PfAMA1_MAL_PfGlurpR2_1tc_dielmo.csv")
list_results_2Ab <- lapply(my_files_2Ab_tempo, read.csv)
my_files_2Ab_tempo <- str_sub(my_files_2Ab_tempo, end=-16)


# Download validation data
Senegalese <- read.csv(file = "Senegalese_data_NbCases.csv") 
Senegalese <- Senegalese[,-1]
Senegalese$age <- c(seq(26, 4, - 1))

for (combinaison in 1:length(my_files_2Ab_tempo)){
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
  #write.csv(MSE_tempo, paste0("MSE_", Ab1, "_", Ab2, "_1tc_.csv"))
  
  
  
  
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
    geom_text(mapping=aes(x=1995.2, y=25, label="1"), size=8, vjust=-0.4, hjust=0)+
    geom_text(mapping=aes(x=2003.2, y=25, label="2"), size=8, vjust=-0.4, hjust=0)+
    geom_text(mapping=aes(x=2006.2, y=25, label="3"), size=8, vjust=-0.4, hjust=0)+
    geom_text(mapping=aes(x=2008.2, y=25, label="4"), size=8, vjust=-0.4, hjust=0)+
    theme_light()+ 
    scale_x_continuous(breaks=seq(1990, 2016, 5)) +
    scale_y_continuous(name = "Number of cases per person per year",
                       sec.axis = sec_axis( trans=~./S_1tc, name="Sero-conversion rate"), limits = c(0,30)) +
    xlab("Year") + ylab("Number of cases per person per year")+  
    labs(title = "b. Estimated sero-conversion rate", 
         caption = "   1: Introduction of Chloroquine
   2: Introduction of AQ + SP
   3: Introduction of AQ + AS
   4: Introduction of LLINs")+
    theme(plot.title = element_text(size=27, face = "bold"),
          plot.subtitle = element_text(size=20),
          axis.text=element_text(size=20),
          axis.title = element_text(size=25, face = "bold"), 
          plot.caption = element_text(hjust=0, size = 25))
  
  #ggsave(filename = paste0(str_sub(Ab1, end=-2), "_", Ab2, "_validation.jpeg"), plot_validation, width=10, height=8)
  

  
  
  
  
}














################
### PLOT VPC ###
################
# Cohort already in nMFI
dielmo_Ndiop <- read.csv("df_Dielmo_NDiop_2016_2018_nMFI.csv")
vec_malaria <- names(dielmo_Ndiop)[grepl("MAL_Pf", names(dielmo_Ndiop))]
vec_epidemio <- colnames(dielmo_Ndiop)[2:8]
dielmo_Ndiop <- dielmo_Ndiop[, c(vec_epidemio, vec_malaria)]

Dielmo <- dielmo_Ndiop[dielmo_Ndiop$village == "Dielmo",]
Dielmo <- Dielmo[order(Dielmo$age),]



# Binarization of data
cutoffs_dielmo <- read.csv("cutoff_df.csv")
for (i in 1:length(vec_malaria)){
  Dielmo[,vec_malaria[i]] <- ifelse(log(Dielmo[, vec_malaria[i]]) >= cutoffs_dielmo$cutoff[cutoffs_dielmo$Ab == vec_malaria[i]], 1, 0)
}


# Split in two cohorts
Dielmo_2016 <- Dielmo[Dielmo$year == "2016",]
Dielmo_2018 <- Dielmo[Dielmo$year == "2018",]


# download results
dfanalyses <- read.csv("MAL_PfAMA1_MAL_PfGlurpR2_1tc_dielmo.csv")
N_par <- 7

##############################
## ESTIMATED SEROPREVALENCE ##
##############################
# ODE systems
equa_noChange <- function(Time, State, Pars){
  with(as.list(c(State,Pars)), {
    SPA <- State[1]
    SPB <- State[2]
    SPAB <- State[3]
    SP0 <- State[4]
    
    dSP0 <- rhoA*SPA + rhoB*SPB + (-lambda+lambda*(1-gammaA)*(1-gammaB))*SP0
    dSPA <- lambda*gammaA*(1-gammaB)*SP0 + rhoB*SPAB - (rhoA + lambda*gammaB)*SPA
    dSPB <- lambda*(1-gammaA)*gammaB*SP0 + rhoA*SPAB - (rhoB + lambda*gammaA)*SPB 
    dSPAB <- lambda*gammaA*gammaB*SP0 + lambda*gammaB*SPA +lambda*gammaA*SPB - (rhoA+rhoB)*SPAB
    
    return(list(c(dSPA, dSPB, dSPAB, dSP0)))
  })
}

equa_Change <- function(Time, State, Pars){
  with(as.list(c(State,Pars)), {
    SPA <- State[1]
    SPB <- State[2]
    SPAB <- State[3]
    SP0 <- State[4]
    
    dSP0 <- rhoA*SPA + rhoB*SPB + (-(delta*lambda)+(delta*lambda)*(1-gammaA)*(1-gammaB))*SP0
    dSPA <- (delta*lambda)*gammaA*(1-gammaB)*SP0 + rhoB*SPAB - (rhoA + (delta*lambda)*gammaB)*SPA
    dSPB <- (delta*lambda)*(1-gammaA)*gammaB*SP0 + rhoA*SPAB - (rhoB + (delta*lambda)*gammaA)*SPB 
    dSPAB <- (delta*lambda)*gammaA*gammaB*SP0 + (delta*lambda)*gammaB*SPA +(delta*lambda)*gammaA*SPB - (rhoA+rhoB)*SPAB
    
    return(list(c(dSPA, dSPB, dSPAB, dSP0)))
  })
}


# Initial probabilities
yini <- c(SPA = 0, SPB = 0, SPAB = 0, SP0 = 1)





# Median estimated seroprevalence
age_seq <- seq(from=0, to=95, by=0.1)                                            # age range 
par_median <- apply(X=dfanalyses[,1:N_par], MARGIN=2, FUN=median)                # median of each parameters posterior distribution
tc <- par_median[7]


age_change_2016 <- age_seq[age_seq <= tc]                                        # age before drop of transmission
age_nochange_2016 <- age_seq[age_seq> tc] - tc                                   # age after drop of transmission

Est_median_2016 <- as.data.frame(ode(yini, age_change_2016, equa_Change, par_median))   
Est_median_2016_bis <- as.data.frame(ode(as.numeric(Est_median_2016[length(age_change_2016),c(2:5)]), age_nochange_2016, equa_noChange, par_median))
colnames(Est_median_2016_bis) <- colnames(Est_median_2016)
Est_median_2016 <- rbind(Est_median_2016, Est_median_2016_bis)

for (i in 1:length(Est_median_2016$time)){
  Est_median_2016[i,"A"] <- Est_median_2016[i,2] + Est_median_2016[i,4]
  Est_median_2016[i,"B"] <- Est_median_2016[i,3] + Est_median_2016[i,4]
}



age_change_2018 <- age_seq[age_seq <= (tc+2)]                                        
age_nochange_2018 <- age_seq[age_seq> (tc+2)] - (tc+2)                                   

Est_median_2018 <- as.data.frame(ode(yini, age_change_2018, equa_Change, par_median))   
Est_median_2018_bis <- as.data.frame(ode(as.numeric(Est_median_2018[length(age_change_2018),c(2:5)]), age_nochange_2018, equa_noChange, par_median))
colnames(Est_median_2018_bis) <- colnames(Est_median_2018)
Est_median_2018 <- rbind(Est_median_2018, Est_median_2018_bis)

for (i in 1:length(Est_median_2018$time)){
  Est_median_2018[i,"A"] <- Est_median_2018[i,2] + Est_median_2018[i,4]
  Est_median_2018[i,"B"] <- Est_median_2018[i,3] + Est_median_2018[i,4]
}





# Confidence Intervalle
N_sam = 700
sam_seq = round(seq(from=1, to=nrow(dfanalyses), length=N_sam))


Est_IC_2016_A = as.data.frame(matrix(NA, nrow=N_sam, ncol=length(age_seq)))
Est_IC_2016_B = as.data.frame(matrix(NA, nrow=N_sam, ncol=length(age_seq)))
Est_IC_2016_AB = as.data.frame(matrix(NA, nrow=N_sam, ncol=length(age_seq)))
for(k in 1:N_sam)
{
  param_tempo <- dfanalyses[sam_seq[k],1:N_par]
  
  tc <- as.numeric(param_tempo[7])
  age_change1 <- age_seq[age_seq <= tc]
  age_nochange <- c(0.000001, age_seq[age_seq> tc] - tc, 100-tc)
  
  
  tempo2 <- as.data.frame(ode(yini, age_change1, equa_Change, param_tempo))   
  tempo <- as.data.frame(ode(as.numeric(tempo2[length(age_change1),c(2:5)]), age_nochange, equa_noChange, param_tempo))
  colnames(tempo) <- colnames(tempo2)
  tempo <- rbind(tempo2, tempo[-c(1,length(tempo$time)),])
  
  
  for (l in 1:length(tempo$time)){
    tempo[l,"A"] <- tempo[l,2] +tempo[l,4]
    tempo[l,"B"] <- tempo[l,3] +tempo[l,4]
  }
  
  
  Est_IC_2016_A[k,] <- t(tempo$A)
  Est_IC_2016_B[k,] <- t(tempo$B)
  Est_IC_2016_AB[k,] <- t(tempo$SPAB)
}


IC_2016_A = matrix(NA, nrow=3, ncol=length(age_seq))
IC_2016_B = matrix(NA, nrow=3, ncol=length(age_seq))
IC_2016_AB = matrix(NA, nrow=3, ncol=length(age_seq))
for(j in 1:length(age_seq))
{
  IC_2016_A[,j] = quantile( Est_IC_2016_A[,j], prob=c(0.025, 0.5, 0.975), na.rm = TRUE )
  IC_2016_B[,j] = quantile( Est_IC_2016_B[,j], prob=c(0.025, 0.5, 0.975), na.rm = TRUE )
  IC_2016_AB[,j] = quantile( Est_IC_2016_AB[,j], prob=c(0.025, 0.5, 0.975), na.rm = TRUE )
}


Est_2016_A<- as.data.frame(t(IC_2016_A))
Est_2016_A$med <- Est_median_2016$A
Est_2016_A$time <- age_seq
Est_2016_A$cohort <- "2016"

Est_2016_B<- as.data.frame(t(IC_2016_B))
Est_2016_B$med <- Est_median_2016$B
Est_2016_B$time <- age_seq
Est_2016_B$cohort <- "2016"

Est_2016_AB<- as.data.frame(t(IC_2016_AB))
Est_2016_AB$med <- Est_median_2016$SPAB
Est_2016_AB$time <- age_seq
Est_2016_AB$cohort <- "2016"





Est_IC_2018_A = as.data.frame(matrix(NA, nrow=N_sam, ncol=length(age_seq)))
Est_IC_2018_B = as.data.frame(matrix(NA, nrow=N_sam, ncol=length(age_seq)))
Est_IC_2018_AB = as.data.frame(matrix(NA, nrow=N_sam, ncol=length(age_seq)))
for(k in 1:N_sam)
{
  param_tempo <- dfanalyses[sam_seq[k],1:N_par]
  
  tc <- as.numeric(param_tempo[7])
  age_change1 <- age_seq[age_seq <= tc]
  age_nochange <- c(0.000001, age_seq[age_seq> tc] - tc, 100-tc)
  
  
  tempo2 <- as.data.frame(ode(yini, age_change1, equa_Change, param_tempo))   
  tempo <- as.data.frame(ode(as.numeric(tempo2[length(age_change1),c(2:5)]), age_nochange, equa_noChange, param_tempo))
  colnames(tempo) <- colnames(tempo2)
  tempo <- rbind(tempo2, tempo[-c(1,length(tempo$time)),])
  
  
  for (l in 1:length(tempo$time)){
    tempo[l,"A"] <- tempo[l,2] +tempo[l,4]
    tempo[l,"B"] <- tempo[l,3] +tempo[l,4]
  }
  
  
  Est_IC_2018_A[k,] <- t(tempo$A)
  Est_IC_2018_B[k,] <- t(tempo$B)
  Est_IC_2018_AB[k,] <- t(tempo$SPAB)
}


IC_2018_A = matrix(NA, nrow=3, ncol=length(age_seq))
IC_2018_B = matrix(NA, nrow=3, ncol=length(age_seq))
IC_2018_AB = matrix(NA, nrow=3, ncol=length(age_seq))
for(j in 1:length(age_seq))
{
  IC_2018_A[,j] = quantile( Est_IC_2018_A[,j], prob=c(0.025, 0.5, 0.975), na.rm = TRUE )
  IC_2018_B[,j] = quantile( Est_IC_2018_B[,j], prob=c(0.025, 0.5, 0.975), na.rm = TRUE )
  IC_2018_AB[,j] = quantile( Est_IC_2018_AB[,j], prob=c(0.025, 0.5, 0.975), na.rm = TRUE )
}


Est_2018_A<- as.data.frame(t(IC_2018_A))
Est_2018_A$med <- Est_median_2018$A
Est_2018_A$time <- age_seq
Est_2018_A$cohort <- "2018"

Est_2018_B<- as.data.frame(t(IC_2018_B))
Est_2018_B$med <- Est_median_2018$B
Est_2018_B$time <- age_seq
Est_2018_B$cohort <- "2018"

Est_2018_AB<- as.data.frame(t(IC_2018_AB))
Est_2018_AB$med <- Est_median_2018$SPAB
Est_2018_AB$time <- age_seq
Est_2018_AB$cohort <- "2018"


Est_final_A <- rbind(Est_2016_A, Est_2018_A)
Est_final_A$Data <- "Estimated" 

Est_final_B <- rbind(Est_2016_B, Est_2018_B)
Est_final_B$Data <- "Estimated" 

Est_final_AB <- rbind(Est_2016_AB, Est_2018_AB)
Est_final_AB$Data <- "Estimated" 

#write.csv(Est_final,paste0("Est_seroprev_", vec_malaria[Ab], "_1tc.csv"))



##############
## PLOTTING ##
##############
# Prepare obs data for plotting 
age_bins     <- seq(from=0, to=95, by=5)
age_bins_mid <- seq(from=2.5, to=92.5, by=5) 
N_bins <- length(age_bins) - 1 
ab1 <- "MAL_PfAMA1"
ab2 <- "MAL_PfGlurpR2"
A1 <- "PfAMA1"
A2 <- "PfGlurpR2"


age_cat <- as.data.frame(c(0.0001,seq(1,65,1), seq(65,85, 5), 100))
colnames(age_cat) <- c("age")
statut_2016 <- as.data.frame(matrix(ncol = 4))
colnames(statut_2016) <- c("SP0", "SPA", "SPB", "SPAB")
for(k in 1: dim(age_cat)[1]){
  tempo <- Dielmo_2016[Dielmo_2016$age>=age_cat$age[k] & Dielmo_2016$age<age_cat$age[k+1],]
  
  age_cat$nb_2016[k] <- dim(tempo)[1]
  
  statut_2016[k,1] <- length(tempo$age[tempo[,ab1] == 0 & tempo[,ab2] == 0])
  statut_2016[k,2] <- length(tempo$age[tempo[,ab1] == 1 & tempo[,ab2] == 0])
  statut_2016[k,3] <- length(tempo$age[tempo[,ab1] == 0 & tempo[,ab2] == 1])
  statut_2016[k,4] <- length(tempo$age[tempo[,ab1] == 1 & tempo[,ab2] == 1])
}

statut_2016$tot <- as.numeric(statut_2016$SP0)+as.numeric(statut_2016$SPA)+as.numeric(statut_2016$SPB)+as.numeric(statut_2016$SPAB)
statut_2016$A <- (as.numeric(statut_2016$SPA)+as.numeric(statut_2016$SPAB))
statut_2016$B <- (as.numeric(statut_2016$SPB)+as.numeric(statut_2016$SPAB))
statut_2016$age <- age_cat$age


SP_bins2016<- matrix(NA, nrow=N_bins, ncol=9)
colnames(SP_bins2016) <- c("medA", "low_95A", "high_95A", "medB", "low_95B", "high_95B", "medAB", "low_95AB", "high_95AB")
for(i in 1:N_bins)
{
  index <- which( statut_2016$age>age_bins[i] & statut_2016$age<=age_bins[i+1] )
  temp_A  <- statut_2016[index,'A']    
  temp_B  <- statut_2016[index,'B']    
  temp_AB  <- statut_2016[index,'SPAB']
  tot_temp <- statut_2016[index,'tot']
  
  SP_bins2016[i,] <- c(as.numeric(as.vector(binom.confint( sum(temp_A), sum(tot_temp), method="wilson")[1,4:6])),
                       as.numeric(as.vector(binom.confint( sum(temp_B), sum(tot_temp), method="wilson")[1,4:6])),   
                       as.numeric(as.vector(binom.confint( sum(temp_AB), sum(tot_temp), method="wilson")[1,4:6])))
}

SP_bins2016 <- as.data.frame(SP_bins2016)
SP_bins2016$age <- age_bins_mid
SP_bins2016$cohort <- "2016"



statut_2018 <- as.data.frame(matrix(ncol = 4))
colnames(statut_2018) <- c("SP0", "SPA", "SPB", "SPAB")
for(k in 1: dim(age_cat)[1]){
  tempo <- Dielmo_2018[Dielmo_2018$age>=age_cat$age[k] & Dielmo_2018$age<age_cat$age[k+1],]
  
  age_cat$nb_2018[k] <- dim(tempo)[1]
  
  statut_2018[k,1] <- length(tempo$age[tempo[,ab1] == 0 & tempo[,ab2] == 0])
  statut_2018[k,2] <- length(tempo$age[tempo[,ab1] == 1 & tempo[,ab2] == 0])
  statut_2018[k,3] <- length(tempo$age[tempo[,ab1] == 0 & tempo[,ab2] == 1])
  statut_2018[k,4] <- length(tempo$age[tempo[,ab1] == 1 & tempo[,ab2] == 1])
}

statut_2018$tot <- as.numeric(statut_2018$SP0)+as.numeric(statut_2018$SPA)+as.numeric(statut_2018$SPB)+as.numeric(statut_2018$SPAB)
statut_2018$A <- (as.numeric(statut_2018$SPA)+as.numeric(statut_2018$SPAB))
statut_2018$B <- (as.numeric(statut_2018$SPB)+as.numeric(statut_2018$SPAB))
statut_2018$age <- age_cat$age


SP_bins2018<- matrix(NA, nrow=N_bins, ncol=9)
colnames(SP_bins2018) <- c("medA", "low_95A", "high_95A", "medB", "low_95B", "high_95B", "medAB", "low_95AB", "high_95AB")
for(i in 1:N_bins)
{
  index <- which( statut_2018$age>age_bins[i] & statut_2018$age<=age_bins[i+1] )
  temp_A  <- statut_2018[index,'A']    
  temp_B  <- statut_2018[index,'B']    
  temp_AB  <- statut_2018[index,'SPAB']    # change Zika to PGP3
  tot_temp <- statut_2018[index,'tot']
  
  SP_bins2018[i,] <- c(as.numeric(as.vector(binom.confint( sum(temp_A), sum(tot_temp), method="wilson")[1,4:6])),
                       as.numeric(as.vector(binom.confint( sum(temp_B), sum(tot_temp), method="wilson")[1,4:6])),   
                       as.numeric(as.vector(binom.confint( sum(temp_AB), sum(tot_temp), method="wilson")[1,4:6])))
}

SP_bins2018 <- as.data.frame(SP_bins2018)
SP_bins2018$age <- age_bins_mid
SP_bins2018$cohort <- "2018"


SP_bins <- rbind(SP_bins2016, SP_bins2018)



plot_A <- ggplot(SP_bins, aes(x = age, y = medA, group = cohort)) +
  geom_point(data = SP_bins, aes(color = cohort), size = 2.5, shape = 19, position=position_dodge(1))+
  geom_errorbar(data = SP_bins, aes(ymin=low_95A, ymax=high_95A,  group = cohort, color = cohort), width=0, alpha = 0.5,
                position=position_dodge(1)) +
  geom_ribbon(data = Est_final_A, mapping = aes(x = time, ymin = V1, ymax = V3, group = cohort, color = cohort, fill = cohort), alpha = 0.2, inherit.aes = F)+
  geom_line(data = Est_final_A, mapping = aes(x = time, y = med,  group = cohort, color = cohort), alpha = 1, linewidth = 2, inherit.aes = F)+
  scale_fill_manual(values = c("#990000", "#FF9999")) +
  scale_color_manual(values = c("#990000", "#FF9999")) +
  theme_light()+ 
  ggtitle(paste0("a. Estimated seroprevalence")) +
  xlab("Age (Years)") + ylab(paste0("Seroprevalence to ", A1))+
  theme(plot.title = element_text(size=27, face = "bold", hjust = 0.9),
        axis.text=element_text(size=20),
        axis.title=element_text(size=25, face = "bold"),
        legend.text = element_text(size=24),
        legend.ticks = element_line(linewidth = 22),
        legend.title = element_text(size=24),
        legend.key.size = unit(1.5,"line"))+
  ylim(0,1)
plot_A


plot_B <- ggplot(SP_bins, aes(x = age, y = medB, group = cohort)) +
  geom_point(data = SP_bins, aes(color = cohort), size = 2.5, shape = 19, position=position_dodge(1))+
  geom_errorbar(data = SP_bins, aes(ymin=low_95B, ymax=high_95B,  group = cohort, color = cohort), width=0, alpha = 0.5,
                position=position_dodge(1)) +
  geom_ribbon(data = Est_final_B, mapping = aes(x = time, ymin = V1, ymax = V3, group = cohort, color = cohort, fill = cohort), alpha = 0.2, inherit.aes = F)+
  geom_line(data = Est_final_B, mapping = aes(x = time, y = med,  group = cohort, color = cohort), alpha = 1, linewidth = 2, inherit.aes = F)+
  scale_fill_manual(values = c("#990000", "#FF9999")) +
  scale_color_manual(values = c("#990000", "#FF9999")) +
  theme_light()+ 
  ggtitle(paste0(" ")) +
  xlab("Age (Years)") + ylab(paste0("Seroprevalence to ", A2))+
  theme(plot.title = element_text(size=27, face = "bold"),
        axis.text=element_text(size=20),
        axis.title=element_text(size=25, face = "bold"),
        legend.text = element_text(size=24),
        legend.ticks = element_line(linewidth = 22),
        legend.title = element_text(size=24),
        legend.key.size = unit(1.5,"line"))+
  ylim(0,1)
plot_B




plot_SP <- ggarrange(plot_A, plot_B, ncol = 2, common.legend = TRUE, legend = "right") 


plot_SP_validation <- ggarrange(plot_SP, "", plot_validation, ncol = 1, heights = c(1, 0.1, 1))


plot_total <- ggarrange(plot_SP_validation, "", plot_final, ncol = 3, widths = c(1, 0.2, 1))
plot_total



ggsave(filename = paste0("Figure3.jpeg"), plot_total, width=27, height=19)



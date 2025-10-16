rm(list = ls())
library(ggplot2)
library(ggpubr)
library(base)
library("deSolve")
library(binom)
library(stringr)



#----------#
#-- Data --#
#----------#
# Download Data 
dielmo_Ndiop <- read.csv("0.Dielmo_Ndiop_PFalci.csv")
vec_malaria <- names(dielmo_Ndiop)[grepl("MAL_Pf", names(dielmo_Ndiop))]
vec_epidemio <- colnames(dielmo_Ndiop)[1:7]

Dielmo <- dielmo_Ndiop[dielmo_Ndiop$village == "Dielmo",]
Dielmo <- Dielmo[order(Dielmo$age),]



# Binarization of data
cutoffs_dielmo <- read.csv("cutoff_df.csv")
for (i in 1:length(vec_malaria)){
  Dielmo[,vec_malaria[i]] <- ifelse(log(Dielmo[, vec_malaria[i]]) >= cutoffs_dielmo$cutoff[cutoffs_dielmo$Ab == vec_malaria[i]], 1, 0)
}


# Split in two cohorts
dielmo_2016 <- Dielmo[Dielmo$year == "2016",]
dielmo_2018 <- Dielmo[Dielmo$year == "2018",]


# download results
results_3AB_1tc = list.files(pattern = "1tc_dielmo.csv")
list_results_3AB_1tc <- lapply(results_3AB_1tc, read.csv)



for (combinaison in 1:length(list_results_3AB_1tc)){
  # Antigens names and place 
  Ab1_str <- str_split(results_3AB_1tc[[combinaison]], "_")[[1]][2]
  Ab2_str <- str_split(results_3AB_1tc[[combinaison]], "_")[[1]][4]
  Ab3_str <- str_split(results_3AB_1tc[[combinaison]], "_")[[1]][6]
  if(Ab3_str == "PfMSP2"){Ab3_str <- "Dd2"}
  if(Ab2_str == "PfMSP2"){Ab2_str <- "Dd2"}
  if(Ab1_str == "PfMSP2"){Ab1_str <- "Dd2"}
  if(Ab2_str == "MAL"){Ab2_str <- str_split(results_3AB_1tc[[combinaison]], "_")[[1]][5]}
  if(Ab3_str == "MAL"){Ab3_str <- str_split(results_3AB_1tc[[combinaison]], "_")[[1]][7]}
  if(Ab3_str == "Ag1"){Ab3_str <- str_split(results_3AB_1tc[[combinaison]], "_")[[1]][8]}
  
  Ab1 <- which(grepl(Ab1_str, vec_malaria))
  Ab2 <-  which(grepl(Ab2_str, vec_malaria))
  Ab3 <-  which(grepl(Ab3_str, vec_malaria))
  
  
  
  # Select results of interest and distinguish chains
  dfanalyses <- list_results_3AB_1tc[[combinaison]]
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
  
  plot_chains1[[param]]<- ggplot(df_tempo) + geom_density(aes(x = parametre, group = chain, color = chain), linewidth = 1.5, alpha = 0.2) +
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

plot_chains <- ggarrange(plot_chains1[[1]], plot_chains1[[2]], plot_chains1[[3]], plot_chains1[[4]], plot_chains1[[5]], plot_chains1[[6]], plot_chains1[[7]],plot_chains1[[8]],plot_chains1[[9]],
                  plot_chains2[[1]], plot_chains2[[2]], plot_chains2[[3]], plot_chains2[[4]], plot_chains2[[5]], plot_chains2[[6]], plot_chains2[[7]],plot_chains2[[8]],  plot_chains2[[9]],
                  ncol = 9, nrow = 2, common.legend = TRUE)

 ggsave(filename = paste0(vec_malaria[Ab1], "_", vec_malaria[Ab2], "_", vec_malaria[Ab3], "_1tc_chains_dielmo.jpeg"), plot_chains, width=22, height=8)






N_par <- 9   
N_mcmc <- dim(dfanalyses)[1]




##############################
## ESTIMATED SEROPREVALENCE ##
##############################
# ODE systems
equa_noChange <- function(Time, State, Pars){
  with(as.list(c(State,Pars)), {
    SP0 <- State[1]
    SPA <- State[2]
    SPB <- State[3]
    SPC <- State[4]
    SPAB <- State[5]
    SPAC <- State[6]
    SPBC <- State[7]
    SPABC <- State[8]
    
    dSP0 <- rhoA*SPA + rhoB*SPB + rhoC*SPC + (-lambda+lambda*(1 - gammaA)*(1 - gammaB)*(1 - gammaC))*SP0
    dSPA <- lambda*gammaA*(1 - gammaB)*(1 - gammaC)*SP0 + rhoB*SPAB + rhoC*SPAC + (-lambda + lambda*(1 - gammaB)*(1 - gammaC) - rhoA)*SPA
    dSPB <- lambda*(1 - gammaA)*gammaB*(1 - gammaC)*SP0 + rhoA*SPAB + rhoC*SPBC + (-lambda + lambda*(1 - gammaA)*(1 - gammaC) - rhoB)*SPB
    dSPC <- lambda*(1 - gammaA)*(1 - gammaB)*gammaC*SP0 + rhoA*SPAC + rhoB*SPBC + (-lambda + lambda*(1 - gammaA)*(1 - gammaB) - rhoC)*SPC
    dSPAB <- lambda*gammaA*gammaB*(1 - gammaC)*SP0 + lambda*gammaA*(1 - gammaC)*SPB + lambda*gammaB*(1 - gammaC)*SPA + rhoC*SPABC + (-lambda + lambda*(1 - gammaC) - rhoB - rhoA)*SPAB
    dSPAC <- lambda*gammaA*(1 - gammaB)*gammaC*SP0 + lambda*gammaA*(1 - gammaB)*SPC + lambda*(1 - gammaB)*gammaC*SPA + rhoB*SPABC + (-lambda + lambda*(1 - gammaB) - rhoA - rhoC)*SPAC
    dSPBC <- lambda*(1 - gammaA)*gammaB*gammaC*SP0 + lambda*(1 - gammaA)*gammaC*SPB + lambda*(1 - gammaA)*gammaB*SPC + rhoA*SPABC + (-lambda + lambda*(1 - gammaA) - rhoB - rhoC)*SPBC
    dSPABC <-  lambda*(gammaA*SPBC + gammaB*SPAC + gammaC*SPAB + gammaA*gammaB*SPC + gammaA*gammaC*SPB + gammaB*gammaC*SPA + gammaA*gammaB*gammaC*SP0) + (-rhoA - rhoB - rhoC)*SPABC
    
    return(list(c(dSP0, dSPA, dSPB, dSPC, dSPAB, dSPAC, dSPBC, dSPABC)))
  })
}



equa_Change <- function(Time, State, Pars){
  with(as.list(c(State,Pars)), {
    SP0 <- State[1]
    SPA <- State[2]
    SPB <- State[3]
    SPC <- State[4]
    SPAB <- State[5]
    SPAC <- State[6]
    SPBC <- State[7]
    SPABC <- State[8]
    
    dSP0 <- rhoA*SPA + rhoB*SPB + rhoC*SPC + (-delta*lambda+delta*lambda*(1 - gammaA)*(1 - gammaB)*(1 - gammaC))*SP0
    dSPA <- delta*lambda*gammaA*(1 - gammaB)*(1 - gammaC)*SP0 + rhoB*SPAB + rhoC*SPAC + (-delta*lambda + delta*lambda*(1 - gammaB)*(1 - gammaC) - rhoA)*SPA
    dSPB <- delta*lambda*(1 - gammaA)*gammaB*(1 - gammaC)*SP0 + rhoA*SPAB + rhoC*SPBC + (-delta*lambda + delta*lambda*(1 - gammaA)*(1 - gammaC) - rhoB)*SPB
    dSPC <- delta*lambda*(1 - gammaA)*(1 - gammaB)*gammaC*SP0 + rhoA*SPAC + rhoB*SPBC + (-delta*lambda + delta*lambda*(1 - gammaA)*(1 - gammaB) - rhoC)*SPC
    dSPAB <- delta*lambda*gammaA*gammaB*(1 - gammaC)*SP0 + delta*lambda*gammaA*(1 - gammaC)*SPB + delta*lambda*gammaB*(1 - gammaC)*SPA + rhoC*SPABC + (-delta*lambda + delta*lambda*(1 - gammaC) - rhoB - rhoA)*SPAB
    dSPAC <- delta*lambda*gammaA*(1 - gammaB)*gammaC*SP0 + delta*lambda*gammaA*(1 - gammaB)*SPC + delta*lambda*(1 - gammaB)*gammaC*SPA + rhoB*SPABC + (-delta*lambda + delta*lambda*(1 - gammaB) - rhoA - rhoC)*SPAC
    dSPBC <- delta*lambda*(1 - gammaA)*gammaB*gammaC*SP0 + delta*lambda*(1 - gammaA)*gammaC*SPB + delta*lambda*(1 - gammaA)*gammaB*SPC + rhoA*SPABC + (-delta*lambda + delta*lambda*(1 - gammaA) - rhoB - rhoC)*SPBC
    dSPABC <-  delta*lambda*(gammaA*SPBC + gammaB*SPAC + gammaC*SPAB + gammaA*gammaB*SPC + gammaA*gammaC*SPB + gammaB*gammaC*SPA + gammaA*gammaB*gammaC*SP0) + (-rhoA - rhoB - rhoC)*SPABC
    
    return(list(c(dSP0, dSPA, dSPB, dSPC, dSPAB, dSPAC, dSPBC, dSPABC)))
  })
}



# Initial probabilities
yini <- c(SP0 = 1, SPA = 0, SPB = 0, SPC = 0, SPAB = 0, SPAC = 0, SPBC = 0, SPABC = 0)






# Median estimated seroprevalence
age_seq <- seq(from=0, to=95, by=0.1)                                            # age range 
par_median <- apply(X=dfanalyses[,1:N_par], MARGIN=2, FUN=median)                # median of each parameters posterior distribution
tc <- par_median[9]


age_change_2016 <- age_seq[age_seq <= tc]                                        # age before drop of transmission
age_nochange_2016 <- age_seq[age_seq> tc] - tc                                   # age after drop of transmission

Est_median_2016 <- as.data.frame(ode(yini, age_change_2016, equa_Change, par_median))   
Est_median_2016_bis <- as.data.frame(ode(as.numeric(Est_median_2016[length(age_change_2016),c(2:9)]), age_nochange_2016, equa_noChange, par_median))
colnames(Est_median_2016_bis) <- colnames(Est_median_2016)
Est_median_2016 <- rbind(Est_median_2016, Est_median_2016_bis)

for (i in 1:length(Est_median_2016$time)){
  Est_median_2016[i,"A"] <- Est_median_2016[i,3] + Est_median_2016[i,6] + Est_median_2016[i,7] + Est_median_2016[i,9]
  Est_median_2016[i,"B"] <- Est_median_2016[i,4] + Est_median_2016[i,6] + Est_median_2016[i,8] + Est_median_2016[i,9]
  Est_median_2016[i,"C"] <- Est_median_2016[i,5] + Est_median_2016[i,7] + Est_median_2016[i,8] + Est_median_2016[i,9]
}



age_change_2018 <- age_seq[age_seq <= (tc+2)]                                        
age_nochange_2018 <- age_seq[age_seq> (tc+2)] - (tc+2)                                   

Est_median_2018 <- as.data.frame(ode(yini, age_change_2018, equa_Change, par_median))   
Est_median_2018_bis <- as.data.frame(ode(as.numeric(Est_median_2018[length(age_change_2018),c(2:9)]), age_nochange_2018, equa_noChange, par_median))
colnames(Est_median_2018_bis) <- colnames(Est_median_2018)
Est_median_2018 <- rbind(Est_median_2018, Est_median_2018_bis)

for (i in 1:length(Est_median_2018$time)){
  Est_median_2018[i,"A"] <- Est_median_2018[i,3] + Est_median_2018[i,6] + Est_median_2018[i,7] + Est_median_2018[i,9]
  Est_median_2018[i,"B"] <- Est_median_2018[i,4] + Est_median_2018[i,6] + Est_median_2018[i,8] + Est_median_2018[i,9]
  Est_median_2018[i,"C"] <- Est_median_2018[i,5] + Est_median_2018[i,7] + Est_median_2018[i,8] + Est_median_2018[i,9]
}





# Confidence Intervalle
N_sam = 700
sam_seq = round(seq(from=1, to=nrow(dfanalyses), length=N_sam))


Est_IC_2016_A = as.data.frame(matrix(NA, nrow=N_sam, ncol=length(age_seq)))
Est_IC_2016_B = as.data.frame(matrix(NA, nrow=N_sam, ncol=length(age_seq)))
Est_IC_2016_C = as.data.frame(matrix(NA, nrow=N_sam, ncol=length(age_seq)))
Est_IC_2016_ABC = as.data.frame(matrix(NA, nrow=N_sam, ncol=length(age_seq)))
for(k in 1:N_sam)
{
  param_tempo <- dfanalyses[sam_seq[k],1:N_par]

  tc <- as.numeric(param_tempo[9])
  age_change1 <- age_seq[age_seq <= tc]
  age_nochange <- c(0.000001, age_seq[age_seq> tc] - tc, 100-tc)
  
  
  tempo2 <- as.data.frame(ode(yini, age_change1, equa_Change, param_tempo))   
  tempo <- as.data.frame(ode(as.numeric(tempo2[length(age_change1),c(2:9)]), age_nochange, equa_noChange, param_tempo))
  colnames(tempo) <- colnames(tempo2)
  tempo <- rbind(tempo2, tempo[-c(1,length(tempo$time)),])
  
  
  for (l in 1:length(tempo$time)){
      tempo[l,"A"] <- tempo[l,3] + tempo[l,6] + tempo[l,7] + tempo[l,9]
      tempo[l,"B"] <- tempo[l,4] + tempo[l,6] + tempo[l,8] + tempo[l,9]
      tempo[l,"C"] <- tempo[l,5] + tempo[l,7] + tempo[l,8] + tempo[l,9]
  }


  Est_IC_2016_A[k,] <- t(tempo$A)
  Est_IC_2016_B[k,] <- t(tempo$B)
  Est_IC_2016_C[k,] <- t(tempo$C)
  Est_IC_2016_ABC[k,] <- t(tempo$SPABC)
}


IC_2016_A = matrix(NA, nrow=3, ncol=length(age_seq))
IC_2016_B = matrix(NA, nrow=3, ncol=length(age_seq))
IC_2016_C = matrix(NA, nrow=3, ncol=length(age_seq))
IC_2016_ABC = matrix(NA, nrow=3, ncol=length(age_seq))
for(j in 1:length(age_seq))
{
  IC_2016_A[,j] = quantile( Est_IC_2016_A[,j], prob=c(0.025, 0.5, 0.975), na.rm = TRUE )
  IC_2016_B[,j] = quantile( Est_IC_2016_B[,j], prob=c(0.025, 0.5, 0.975), na.rm = TRUE )
  IC_2016_C[,j] = quantile( Est_IC_2016_C[,j], prob=c(0.025, 0.5, 0.975), na.rm = TRUE )
  IC_2016_ABC[,j] = quantile( Est_IC_2016_ABC[,j], prob=c(0.025, 0.5, 0.975), na.rm = TRUE )
}


Est_2016_A<- as.data.frame(t(IC_2016_A))
Est_2016_A$med <- Est_median_2016$A
Est_2016_A$time <- age_seq
Est_2016_A$cohort <- "2016"

Est_2016_B<- as.data.frame(t(IC_2016_B))
Est_2016_B$med <- Est_median_2016$B
Est_2016_B$time <- age_seq
Est_2016_B$cohort <- "2016"

Est_2016_C<- as.data.frame(t(IC_2016_C))
Est_2016_C$med <- Est_median_2016$C
Est_2016_C$time <- age_seq
Est_2016_C$cohort <- "2016"

Est_2016_ABC<- as.data.frame(t(IC_2016_ABC))
Est_2016_ABC$med <- Est_median_2016$SPABC
Est_2016_ABC$time <- age_seq
Est_2016_ABC$cohort <- "2016"





Est_IC_2018_A = as.data.frame(matrix(NA, nrow=N_sam, ncol=length(age_seq)))
Est_IC_2018_B = as.data.frame(matrix(NA, nrow=N_sam, ncol=length(age_seq)))
Est_IC_2018_C = as.data.frame(matrix(NA, nrow=N_sam, ncol=length(age_seq)))
Est_IC_2018_ABC = as.data.frame(matrix(NA, nrow=N_sam, ncol=length(age_seq)))
for(k in 1:N_sam)
{
  param_tempo <- dfanalyses[sam_seq[k],1:N_par]
  
  tc <- as.numeric(param_tempo[9])
  age_change1 <- age_seq[age_seq <= tc]
  age_nochange <- c(0.000001, age_seq[age_seq> tc] - tc, 100-tc)
  
  
  tempo2 <- as.data.frame(ode(yini, age_change1, equa_Change, param_tempo))   
  tempo <- as.data.frame(ode(as.numeric(tempo2[length(age_change1),c(2:9)]), age_nochange, equa_noChange, param_tempo))
  colnames(tempo) <- colnames(tempo2)
  tempo <- rbind(tempo2, tempo[-c(1,length(tempo$time)),])
  
  
  for (l in 1:length(tempo$time)){
    tempo[l,"A"] <- tempo[l,3] + tempo[l,6] + tempo[l,7] + tempo[l,9]
    tempo[l,"B"] <- tempo[l,4] + tempo[l,6] + tempo[l,8] + tempo[l,9]
    tempo[l,"C"] <- tempo[l,5] + tempo[l,7] + tempo[l,8] + tempo[l,9]
  }
  
  
  Est_IC_2018_A[k,] <- t(tempo$A)
  Est_IC_2018_B[k,] <- t(tempo$B)
  Est_IC_2018_C[k,] <- t(tempo$C)
  Est_IC_2018_ABC[k,] <- t(tempo$SPABC)
}


IC_2018_A = matrix(NA, nrow=3, ncol=length(age_seq))
IC_2018_B = matrix(NA, nrow=3, ncol=length(age_seq))
IC_2018_C = matrix(NA, nrow=3, ncol=length(age_seq))
IC_2018_ABC = matrix(NA, nrow=3, ncol=length(age_seq))
for(j in 1:length(age_seq))
{
  IC_2018_A[,j] = quantile( Est_IC_2018_A[,j], prob=c(0.025, 0.5, 0.975), na.rm = TRUE )
  IC_2018_B[,j] = quantile( Est_IC_2018_B[,j], prob=c(0.025, 0.5, 0.975), na.rm = TRUE )
  IC_2018_C[,j] = quantile( Est_IC_2018_C[,j], prob=c(0.025, 0.5, 0.975), na.rm = TRUE )
  IC_2018_ABC[,j] = quantile( Est_IC_2018_ABC[,j], prob=c(0.025, 0.5, 0.975), na.rm = TRUE )
}


Est_2018_A<- as.data.frame(t(IC_2018_A))
Est_2018_A$med <- Est_median_2018$A
Est_2018_A$time <- age_seq
Est_2018_A$cohort <- "2018"

Est_2018_B<- as.data.frame(t(IC_2018_B))
Est_2018_B$med <- Est_median_2018$B
Est_2018_B$time <- age_seq
Est_2018_B$cohort <- "2018"

Est_2018_C<- as.data.frame(t(IC_2018_C))
Est_2018_C$med <- Est_median_2018$C
Est_2018_C$time <- age_seq
Est_2018_C$cohort <- "2018"

Est_2018_ABC<- as.data.frame(t(IC_2018_ABC))
Est_2018_ABC$med <- Est_median_2018$SPABC
Est_2018_ABC$time <- age_seq
Est_2018_ABC$cohort <- "2018"



Est_final_A <- rbind(Est_2016_A, Est_2018_A)
Est_final_A$Data <- "Estimated" 

Est_final_B <- rbind(Est_2016_B, Est_2018_B)
Est_final_B$Data <- "Estimated" 

Est_final_C <- rbind(Est_2016_C, Est_2018_C)
Est_final_C$Data <- "Estimated" 

Est_final_ABC <- rbind(Est_2016_ABC, Est_2018_ABC)
Est_final_ABC$Data <- "Estimated" 

#write.csv(Est_final,paste0("Est_seroprev_", vec_malaria[Ab], "_1tc.csv"))



##############
## PLOTTING ##
##############
# Prepare obs data for plotting 
age_bins     <- seq(from=0, to=95, by=5)
age_bins_mid <- seq(from=2.5, to=92.5, by=5) 
N_bins <- length(age_bins) - 1 


age_cat <- as.data.frame(c(0.0001,seq(1,65,1), seq(65,85, 5), 100))
colnames(age_cat) <- c("age")
statut_2016 <- as.data.frame(matrix(ncol = 8))
colnames(statut_2016) <- c("SP0", "SPA", "SPB", "SPC", "SPAB", "SPAC", "SPBC", "SPABC")
for(j in 1:dim(age_cat)[1]){
  tempo <- dielmo_2016[dielmo_2016$age>=age_cat$age[j] & dielmo_2016$age<age_cat$age[j+1],]
  
  age_cat$nb_2016[j] <- dim(tempo)[1]
  
  statut_2016[j,1] <- length(tempo$age[tempo[,Ab1+7] == 0 & tempo[,Ab2+7] == 0 & tempo[,Ab3+7] == 0])
  statut_2016[j,2] <- length(tempo$age[tempo[,Ab1+7] == 1 & tempo[,Ab2+7] == 0 & tempo[,Ab3+7] == 0])
  statut_2016[j,3] <- length(tempo$age[tempo[,Ab1+7] == 0 & tempo[,Ab2+7] == 1 & tempo[,Ab3+7] == 0])
  statut_2016[j,4] <- length(tempo$age[tempo[,Ab1+7] == 0 & tempo[,Ab2+7] == 0 & tempo[,Ab3+7] == 1])
  statut_2016[j,5] <- length(tempo$age[tempo[,Ab1+7] == 1 & tempo[,Ab2+7] == 1 & tempo[,Ab3+7] == 0])
  statut_2016[j,6] <- length(tempo$age[tempo[,Ab1+7] == 1 & tempo[,Ab2+7] == 0 & tempo[,Ab3+7] == 1])
  statut_2016[j,7] <- length(tempo$age[tempo[,Ab1+7] == 0 & tempo[,Ab2+7] == 1 & tempo[,Ab3+7] == 1])
  statut_2016[j,8] <- length(tempo$age[tempo[,Ab1+7] == 1 & tempo[,Ab2+7] == 1 & tempo[,Ab3+7] == 1])
}

statut_2016$tot <- as.numeric(statut_2016$SP0)+as.numeric(statut_2016$SPA)+as.numeric(statut_2016$SPB)+as.numeric(statut_2016$SPC)+as.numeric(statut_2016$SPAB)+as.numeric(statut_2016$SPAC)+as.numeric(statut_2016$SPBC)+as.numeric(statut_2016$SPABC)
statut_2016$A <- as.numeric(statut_2016$SPA)+as.numeric(statut_2016$SPAB)+as.numeric(statut_2016$SPAC)+as.numeric(statut_2016$SPABC)
statut_2016$B <- as.numeric(statut_2016$SPB)+as.numeric(statut_2016$SPAB)+as.numeric(statut_2016$SPBC)+as.numeric(statut_2016$SPABC)
statut_2016$C <- as.numeric(statut_2016$SPC)+as.numeric(statut_2016$SPAC)+as.numeric(statut_2016$SPBC)+as.numeric(statut_2016$SPABC)
statut_2016$age <- age_cat$age


SP_bins2016<- matrix(NA, nrow=N_bins, ncol=12)
colnames(SP_bins2016) <- c("medA", "low_95A", "high_95A", "medB", "low_95B", "high_95B","medC", "low_95C", "high_95C", "medABC", "low_95ABC", "high_95ABC")
for(i in 1:N_bins)
{
  index <- which( statut_2016$age>age_bins[i] & statut_2016$age<=age_bins[i+1] )
  temp_A  <- statut_2016[index,'A']    
  temp_B  <- statut_2016[index,'B']  
  temp_C  <- statut_2016[index,'C']  
  temp_ABC  <- statut_2016[index,'SPABC']    # change Zika to PGP3
  tot_temp <- statut_2016[index,'tot']
  
  SP_bins2016[i,] <- c(as.numeric(as.vector(binom.confint( sum(temp_A), sum(tot_temp), method="wilson")[1,4:6])),
                         as.numeric(as.vector(binom.confint( sum(temp_B), sum(tot_temp), method="wilson")[1,4:6])), 
                         as.numeric(as.vector(binom.confint( sum(temp_C), sum(tot_temp), method="wilson")[1,4:6])), 
                         as.numeric(as.vector(binom.confint( sum(temp_ABC), sum(tot_temp), method="wilson")[1,4:6])))
}

SP_bins2016 <- as.data.frame(SP_bins2016)
SP_bins2016$age <- age_bins_mid
SP_bins2016$cohort <- "2016"



statut_2018 <- as.data.frame(matrix(ncol = 8))
colnames(statut_2018) <- c("SP0", "SPA", "SPB", "SPC", "SPAB", "SPAC", "SPBC", "SPABC")
for(j in 1:dim(age_cat)[1]){
  tempo <- dielmo_2018[dielmo_2018$age>=age_cat$age[j] & dielmo_2018$age<age_cat$age[j+1],]
  
  age_cat$nb_2018[j] <- dim(tempo)[1]
  
  statut_2018[j,1] <- length(tempo$age[tempo[,Ab1+7] == 0 & tempo[,Ab2+7] == 0 & tempo[,Ab3+7] == 0])
  statut_2018[j,2] <- length(tempo$age[tempo[,Ab1+7] == 1 & tempo[,Ab2+7] == 0 & tempo[,Ab3+7] == 0])
  statut_2018[j,3] <- length(tempo$age[tempo[,Ab1+7] == 0 & tempo[,Ab2+7] == 1 & tempo[,Ab3+7] == 0])
  statut_2018[j,4] <- length(tempo$age[tempo[,Ab1+7] == 0 & tempo[,Ab2+7] == 0 & tempo[,Ab3+7] == 1])
  statut_2018[j,5] <- length(tempo$age[tempo[,Ab1+7] == 1 & tempo[,Ab2+7] == 1 & tempo[,Ab3+7] == 0])
  statut_2018[j,6] <- length(tempo$age[tempo[,Ab1+7] == 1 & tempo[,Ab2+7] == 0 & tempo[,Ab3+7] == 1])
  statut_2018[j,7] <- length(tempo$age[tempo[,Ab1+7] == 0 & tempo[,Ab2+7] == 1 & tempo[,Ab3+7] == 1])
  statut_2018[j,8] <- length(tempo$age[tempo[,Ab1+7] == 1 & tempo[,Ab2+7] == 1 & tempo[,Ab3+7] == 1])
}

statut_2018$tot <- as.numeric(statut_2018$SP0)+as.numeric(statut_2018$SPA)+as.numeric(statut_2018$SPB)+as.numeric(statut_2018$SPC)+as.numeric(statut_2018$SPAB)+as.numeric(statut_2018$SPAC)+as.numeric(statut_2018$SPBC)+as.numeric(statut_2018$SPABC)
statut_2018$A <- as.numeric(statut_2018$SPA)+as.numeric(statut_2018$SPAB)+as.numeric(statut_2018$SPAC)+as.numeric(statut_2018$SPABC)
statut_2018$B <- as.numeric(statut_2018$SPB)+as.numeric(statut_2018$SPAB)+as.numeric(statut_2018$SPBC)+as.numeric(statut_2018$SPABC)
statut_2018$C <- as.numeric(statut_2018$SPC)+as.numeric(statut_2018$SPAC)+as.numeric(statut_2018$SPBC)+as.numeric(statut_2018$SPABC)
statut_2018$age <- age_cat$age


SP_bins2018<- matrix(NA, nrow=N_bins, ncol=12)
colnames(SP_bins2018) <- c("medA", "low_95A", "high_95A", "medB", "low_95B", "high_95B","medC", "low_95C", "high_95C", "medABC", "low_95ABC", "high_95ABC")
for(i in 1:N_bins)
{
  index <- which( statut_2018$age>age_bins[i] & statut_2018$age<=age_bins[i+1] )
  temp_A  <- statut_2018[index,'A']    
  temp_B  <- statut_2018[index,'B']  
  temp_C  <- statut_2018[index,'C']  
  temp_ABC  <- statut_2018[index,'SPABC']    # change Zika to PGP3
  tot_temp <- statut_2018[index,'tot']
  
  SP_bins2018[i,] <- c(as.numeric(as.vector(binom.confint( sum(temp_A), sum(tot_temp), method="wilson")[1,4:6])),
                       as.numeric(as.vector(binom.confint( sum(temp_B), sum(tot_temp), method="wilson")[1,4:6])), 
                       as.numeric(as.vector(binom.confint( sum(temp_C), sum(tot_temp), method="wilson")[1,4:6])), 
                       as.numeric(as.vector(binom.confint( sum(temp_ABC), sum(tot_temp), method="wilson")[1,4:6])))
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
  theme_classic()+
  ggtitle(paste0(Ab1_str)) +
  xlab("Age (Years)") + ylab("Seroprevalence")+
  theme(plot.title = element_text(size=22, face = "bold"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=16, face = "bold"))+
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
  theme_classic()+
  ggtitle(paste0(Ab2_str)) +
  xlab("Age (Years)") + ylab("Seroprevalence")+
  theme(plot.title = element_text(size=22, face = "bold"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=16, face = "bold"))+
  ylim(0,1)
plot_B


plot_C <- ggplot(SP_bins, aes(x = age, y = medC, group = cohort)) +
  geom_point(data = SP_bins, aes(color = cohort), size = 2.5, shape = 19, position=position_dodge(1))+
  geom_errorbar(data = SP_bins, aes(ymin=low_95C, ymax=high_95C,  group = cohort, color = cohort), width=0, alpha = 0.5,
                position=position_dodge(1)) +
  geom_ribbon(data = Est_final_C, mapping = aes(x = time, ymin = V1, ymax = V3, group = cohort, color = cohort, fill = cohort), alpha = 0.2, inherit.aes = F)+
  geom_line(data = Est_final_C, mapping = aes(x = time, y = med,  group = cohort, color = cohort), alpha = 1, linewidth = 2, inherit.aes = F)+
  scale_fill_manual(values = c("#990000", "#FF9999")) +
  scale_color_manual(values = c("#990000", "#FF9999")) +
  theme_classic()+
  ggtitle(paste0(Ab3_str)) +
  xlab("Age (Years)") + ylab("Seroprevalence")+
  theme(plot.title = element_text(size=22, face = "bold"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=16, face = "bold"))+
  ylim(0,1)
plot_C


plot_ABC <- ggplot(SP_bins, aes(x = age, y = medABC, group = cohort)) +
  geom_point(data = SP_bins, aes(color = cohort), size = 2.5, shape = 19, position=position_dodge(1))+
  geom_errorbar(data = SP_bins, aes(ymin=low_95ABC, ymax=high_95ABC,  group = cohort, color = cohort), width=0, alpha = 0.5,
                position=position_dodge(1)) +
  geom_ribbon(data = Est_final_ABC, mapping = aes(x = time, ymin = V1, ymax = V3, group = cohort, color = cohort, fill = cohort), alpha = 0.2, inherit.aes = F)+
  geom_line(data = Est_final_ABC, mapping = aes(x = time, y = med,  group = cohort, color = cohort), alpha = 1, linewidth = 2, inherit.aes = F)+
  scale_fill_manual(values = c("#990000", "#FF9999")) +
  scale_color_manual(values = c("#990000", "#FF9999")) +
  theme_classic()+
  ggtitle(paste0(Ab1_str, " ", Ab2_str, " ", Ab3_str)) +
  xlab("Age (Years)") + ylab("Seroprevalence")+
  theme(plot.title = element_text(size=22, face = "bold"),
        axis.text=element_text(size=12),
        axis.title=element_text(size=16, face = "bold"))+
  ylim(0,1)
plot_ABC

plot_SP_1 <- ggarrange(plot_A, plot_B, plot_C, ncol = 1, legend = "none") 
plot_SP <- ggarrange(plot_SP_1, plot_ABC, ncol = 2, common.legend = TRUE, legend = "right") 

ggsave(filename = paste0(vec_malaria[Ab1], "_",vec_malaria[Ab2], "_",vec_malaria[Ab3], "_1tc_VPC_dielmo.jpeg"), plot_SP, width=15, height=10)



}









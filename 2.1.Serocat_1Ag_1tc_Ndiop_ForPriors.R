rm(list=ls())

library("rstan")
library("bayesplot")
library("RColorBrewer")
library("dplyr")
library("binom")
library("parallelly")
library("deSolve")



args <- commandArgs(trailingOnly=TRUE)


######################################################################
# parameters to be used
Ab1 =  as.numeric(args[1])

options(mc.cores = parallel::detectCores())

#----------#
#-- Data --#
#----------#
# Download Data: Samples and negative controls
dielmo_Ndiop <- read.csv("0.Dielmo_Ndiop_PFalci.csv")
vec_malaria <- names(dielmo_Ndiop)[grepl("MAL_Pf", names(dielmo_Ndiop))]
vec_epidemio <- colnames(dielmo_Ndiop)[1:7]

NDiop <- dielmo_Ndiop[dielmo_Ndiop$village == "Ndiop",]




# Binarization of data
cutoffs_dielmo <- read.csv("0.cutoff_df.csv")
for (i in 1:length(vec_malaria)){
  NDiop[,vec_malaria[i]] <- ifelse(log(NDiop[, vec_malaria[i]]) >= cutoffs_dielmo$cutoff[cutoffs_dielmo$Ab == vec_malaria[i]], 1, 0)
}


# Split in two cohorts
ndiop_2016 <- NDiop[NDiop$year == "2016",]
ndiop_2018 <- NDiop[NDiop$year == "2018",]


# Bin age
age_cat <- as.data.frame(c(0.0001,seq(1,65,1), seq(65,85, 5), 100))
colnames(age_cat) <- c("age")
statut_2016 <- as.data.frame(matrix(ncol = 1))
for(j in 1: dim(age_cat)[1]){
  tempo <- ndiop_2016[ndiop_2016$age>=age_cat$age[j] & ndiop_2016$age<age_cat$age[j+1],]
  
  age_cat$nb_2016[j] <- dim(tempo)[1]
  
  statut_2016[j,1] <- length(tempo$age[tempo[,Ab1+7] == 1])
}
statut_2016$age_cat <- age_cat$age
statut_2016$nb <- age_cat$nb_2016
statut_2016 <- statut_2016[statut_2016$nb !=0,]



statut_2018 <- as.data.frame(matrix(ncol = 1))
for(j in 1: dim(age_cat)[1]){
  tempo <- ndiop_2018[ndiop_2018$age>=age_cat$age[j] & ndiop_2018$age<age_cat$age[j+1],]
  
  age_cat$nb_2018[j] <- dim(tempo)[1]
  
  statut_2018[j,1] <- length(tempo$age[tempo[,Ab1+7] == 1])
}

statut_2018$age_cat <- age_cat$age
statut_2018$nb <- age_cat$nb_2018
statut_2018 <- statut_2018[statut_2018$nb !=0,]





#-----------------------#
#-- data for analysis --#
#-----------------------#
listmodel <- list(N_2016 = dim(statut_2016)[1],          # Number of ages / age bin to look at
                  N_2018 = dim(statut_2018)[1],
                  age_2016 = statut_2016$age_cat,        # vector of ages
                  age_2018 = statut_2018$age_cat,
                  n_2016 = statut_2016$nb,               # vector of number of participant by age
                  n_2018 = statut_2018$nb,
                  Statut_2016 = statut_2016[,-c(2,3)],   # matrix of number of participant in each serostatus
                  Statut_2018 = statut_2018[,-c(2,3)],
                  initial_state = c(0,1),                # for ODE: everybody is negative at age 0
                  age0 = 0) 













#-------------------#
#-- Model fitting --#
#-------------------#
# Stan program
model <- stan_model(file="2.1.1Ag_1tc.stan") 


# Initial parameters: need to implement because stan sample in -2:2 otherwise
init_fun <- function() {
  list(param = c(0.05, 0.1, 0.2), time_c = 15)
}



# run model
fit <- sampling(model,                                                           # Stan program
                data = listmodel,                                                # named list of data
                chains = 4,                                                      # number of Markov chains
                warmup = 1000,                                                   # number of warmup iterations per chain 
                iter = 4000,                                                     # total number of iterations per chain 
                refresh = 500,                                                   # show progress every 'refresh' iterations
                init = init_fun,                                                               
                control=list(adapt_delta=0.99,
                             max_treedepth=12),
                verbose = TRUE)


# data management of results
results_model <- extract(fit)[c(1:3)]
tempo <- matrix(unlist(results_model), ncol = 5)
colnames(tempo) <- c("rhoA", "lambda", "delta", "tc", "ll")


# save results as csv
write.csv(tempo, paste0(vec_malaria[Ab1],"_1tc_ndiop.csv"), row.names=FALSE)                                   


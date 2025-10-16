rm(list=ls())
library("rstan")




# Download Data: Samples and negative controls
dielmo_ndiop <- read.csv("0.Dielmo_Ndiop_PFalci.csv")
vec_malaria <- names(dielmo_ndiop)[grepl("MAL_Pf", names(dielmo_ndiop))]
vec_epidemio <- colnames(dielmo_ndiop)[1:7]

neg <- read.csv("0.Seroped_nMFI.csv")



# Loop for each antigen
for (Ag in 1:length(1:vec_malaria)){
  
    #-------------------------#
    #-- data for Management --#
    #-------------------------#
    print("Creating data for model")  
    list_data_model <- list(N_unkown = length(dielmo_ndiop$sample_name),                # Number of unknown samples
                            N_neg =  length(neg$id_sample),                             # Number of negatives samples
                            MFI_unkown = log(dielmo_ndiop[,7+Ag]),                      # Vector with nMFI 
                            MFI_neg = log(neg[,11+Ag]),
                            MFI_max = log(max(dielmo_ndiop[,7+Ag], neg[,11+Ag])),       # Maximal value of nMFI (for priors constraints)
                            MFI_min = log(min(dielmo_ndiop[,7+Ag], neg[,11+Ag])))
    
  
    
    
    
    #-------------------#
    #-- Model fitting --#
    #-------------------#
    # Stan program
    model <- stan_model(file="1.1.GMM_with_confirmed_Neg.stan")
    
    
    # Defining initial parameters (if not stan samples in [-2;2])
    init_fun <- function() {
      list(theta = 0.5, sigmaN = 1, muN = log(min(dielmo_ndiop[,7+Ag], neg[,11+Ag]))+0.1, sigmaP = 1, muP = 2)
    }
    
    
    # Fit model
    print("Fitting model")  
    fit <- sampling(model,                            # Stan program
                      data = list_data_model,         # list of data
                      chains = 1,                     # number of Markov chains
                      warmup = 2500,                  # number of warmup iterations per chain (2500)
                      iter = 10000,                   # total number of iterations per chain (10000)
                      refresh = 200,                  # show progress every 'refresh' iterations
                      init = init_fun,                # Initial values set at random
                      verbose = TRUE                  # Gives more explicit error messages
      )
    
    
    print("Model is done")  
    
    # Data Management of the results
    results_model <- extract(fit)[c(1:6)] # gets the estimated parameters and Log Likelihood for each iterations 
    tempo <- matrix(unlist(results_model), ncol = 6) 
    colnames(tempo) <- c("theta", "sigmaN", "muN", "sigmaP", "muP", "LL")
      
    write.csv(tempo, paste0("results_pooled_",vec_malaria[Ag],"_nMFI.csv"), row.names=FALSE) # saves results in your folder

}



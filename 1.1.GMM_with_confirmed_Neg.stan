
data {
  int<lower=0> N_unkown;       // Number of unkown
  int<lower=0> N_neg;          // Number of negatives
  real MFI_unkown[N_unkown];   // MFI value for unkown
  real MFI_neg[N_neg];         // MFI value for negatives
  real MFI_max;                // maximum value of MFI
  real MFI_min;                // minimum value of MFI (for priors constraint)
}


// The parameters accepted by the model. 
parameters {
  real<lower=0, upper = 1> theta;             // prevalence
  real<lower=0, upper = 100> sigmaN;          // Std deviation of muN     
  real<lower=MFI_min, upper = MFI_max> muN;   // Mean MFI value for Negatives
  real<lower=0, upper = 100> sigmaP;          // Std deviation of muP
  real<lower=0, upper = 100> muP;             // difference between muN and muP, has to be positive to have muPos bigger than muN
}


transformed parameters {
    real muPos;
    
    muPos = muN+muP;  // Mean of the positives
}



// The model to be estimated. 
model {
  // Priors definition
  theta ~ uniform(0,1);
  sigmaN ~ uniform(0,100);
  muN ~ uniform(MFI_min,MFI_max);
  sigmaP ~ uniform(0,100);
  muP ~ normal(log(100), 0.2);  

// likelihood for the unknown
  for (l in 1:N_unkown){
    target += log_mix(theta,
                        normal_lpdf(MFI_unkown[l] | muPos, sigmaP),
                        normal_lpdf(MFI_unkown[l] | muN, sigmaN));
  }

// likelihood for the negatives
  for (m in 1:N_neg){
    target += log_mix(0, // Only considering the negatives
                        normal_lpdf(MFI_neg[m] | muPos, sigmaP),
                        normal_lpdf(MFI_neg[m] | muN, sigmaN));
  }
}







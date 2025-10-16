functions{
  vector ode_noChange(real age,       // time
             vector SP,               // system states
             real[] param)            // parameters
             {
               vector[2] dSPdt;
               dSPdt[1] = param[2]*(1 - SP[1]) - param[1]*SP[1];   //ODE
               dSPdt[2] = (1 - SP[1]);

               return dSPdt;
             }
             
   vector ode_Change1(real age,       // time
             vector SP,               // system states
             real[] param)            // parameters
             {
               vector[2] dSPdt;
               dSPdt[1] = param[3]*param[2]*(1 - SP[1]) - param[1]*SP[1]; //ODE
               dSPdt[2] = (1 - SP[1]);

               return dSPdt;
             }
}




data {
  int<lower=0> N_2016;
  int<lower=0> N_2018;
  real age_2016[N_2016];
  real age_2018[N_2018];
  int n_2016[N_2016];
  int n_2018[N_2018];
  int Statut_2016[N_2016];
  int Statut_2018[N_2018];
  vector[2] initial_state;
  real age0;
  real median_rho;
  real sd_rho;
}



// The parameters accepted by the model. 
parameters {
  real<lower=0> param[3]; // 1. rhoA, 2. lambda, 3.delta
  
  real<lower=0, upper = 100> time_c;   // Time from which the decrease happend
}



// The model to be estimated. 
model {
  param[1]~ normal(median_rho, sd_rho) T[0, ];
  param[2]~ exponential(1);
  param[3]~ uniform(0,2);


  time_c ~ uniform(0,50);

  int N_tempo_2016 = 0;
  int N_tempo_2018 = 0;

  real age_c_2016[N_2016];     // So age start back at 0 after the decrease (for ODE)
  real age_c_2018[N_2018];



  // Create new vector with changed ages for ODE after decrease in transmission
  for (i in 1:N_2016){
     if(age_2016[i] <= time_c){
       age_c_2016[i] = age_2016[i];
     }else if (age_2016[i] > time_c){
       age_c_2016[i] = age_2016[i] - time_c;
     }
  }
  
  for (i in 1:N_2018){
     if(age_2018[i] <= (time_c+2)){
       age_c_2018[i] = age_2018[i];
     }else if (age_2018[i] > (time_c+2)){
       age_c_2018[i] = age_2018[i] - (time_c+2);
     }
  }



  for (i in 1:N_2016){if(age_2016[i]<= time_c){N_tempo_2016 = N_tempo_2016+1;}}                                        // count participant aged younger than the time of decrease in transmission
  array[N_2016] vector[2] theta_2016 = ode_rk45(ode_Change1, initial_state, age0, age_2016, param);                    // apply ODE to all participants (will only look at those younger than decrease in LL)
  array[N_2016 - N_tempo_2016] vector[2] theta_bis_2016 = ode_rk45(ode_noChange, theta_2016[N_tempo_2016], age0, age_c_2016[(N_tempo_2016+1):N_2016], param); // aply ODE for participants older than decrease
  
  for (k in 1:N_2018){if(age_2018[k]<= (time_c+2)){N_tempo_2018 = N_tempo_2018+1;}}
  array[N_2018] vector[2] theta_2018 = ode_rk45(ode_Change1, initial_state, age0, age_2018, param);
  array[N_2018 - N_tempo_2018] vector[2] theta_bis_2018 = ode_rk45(ode_noChange, theta_2018[N_tempo_2018], age0, age_c_2018[(N_tempo_2018+1):N_2018], param);
   

  // Likelihood
  Statut_2016[1:N_tempo_2016] ~ binomial(n_2016[1:N_tempo_2016], theta_2016[1:N_tempo_2016,1]);   // Applied to participants younger than decreased 
  Statut_2016[N_tempo_2016+1:N_2016] ~ binomial(n_2016[N_tempo_2016+1:N_2016], theta_bis_2016[1:(N_2016 - N_tempo_2016),1]);  // applied to participants older than decreased

  Statut_2018[1:N_tempo_2018] ~ binomial(n_2018[1:N_tempo_2018], theta_2018[1:N_tempo_2018,1]);
  Statut_2018[N_tempo_2018+1:N_2018] ~ binomial(n_2018[N_tempo_2018+1:N_2018], theta_bis_2018[1:(N_2018 - N_tempo_2018),1]);

}






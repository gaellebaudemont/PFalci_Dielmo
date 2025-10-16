functions{
  vector ode_noChange(real age,       // time
             vector SP,    // system state
             real[] param)    // parameters
             {
               vector[4] dSPdt;
               dSPdt[1] = param[1]*SP[2] + param[2]*SP[3] + (-param[3]+param[3]*(1-param[5])*(1-param[6]))*SP[1];
               dSPdt[2] = param[3]*param[5]*(1-param[6])*SP[1] + param[2]*SP[4] - (param[1] + param[3]*param[6])*SP[2];
               dSPdt[3] = param[3]*(1-param[5])*param[6]*SP[1] + param[1]*SP[4] - (param[2] + param[3]*param[5])*SP[3];
               dSPdt[4] = param[3]*param[5]*param[6]*SP[1] + param[3]*param[6]*SP[2] +param[3]*param[5]*SP[3] - (param[1]+param[2])*SP[4];

               return dSPdt;
             }
             
   vector ode_Change1(real age,       // time
             vector SP,    // system state
             real[] param)    // parameters
             {
               vector[4] dSPdt;
               dSPdt[1] = param[1]*SP[2] + param[2]*SP[3] + (-(param[3]*param[4])+(param[3]*param[4])*(1-param[5])*(1-param[6]))*SP[1];
               dSPdt[2] = (param[3]*param[4])*param[5]*(1-param[6])*SP[1] + param[2]*SP[4] - (param[1] + (param[3]*param[4])*param[6])*SP[2];
               dSPdt[3] = (param[3]*param[4])*(1-param[5])*param[6]*SP[1] + param[1]*SP[4] - (param[2] + (param[3]*param[4])*param[5])*SP[3];
               dSPdt[4] = (param[3]*param[4])*param[5]*param[6]*SP[1] + (param[3]*param[4])*param[6]*SP[2] +(param[3]*param[4])*param[5]*SP[3] - (param[1]+param[2])*SP[4];

               return dSPdt;
             }
}




data {
  int<lower=0> N_2016;
  int<lower=0> N_2018;
  real age_2016[N_2016];
  real age_2018[N_2018];
  int Statut_2016[N_2016,4];
  int Statut_2018[N_2018,4];
  vector[4] initial_state;
  real age0;
  real median_rhoA;
  real median_rhoB;
  real sd_rhoA;
  real sd_rhoB;
}



// The parameters accepted by the model. 
parameters {
   real<lower=0> param[6]; // 1. rhoA, 2. rhoB, 3.lambda, 4.delta, 5.gammaA, 6.gammaB
 
  real<lower=0, upper = 100> time_c;
}




// The model to be estimated. 
model {
  param[1] ~ normal(median_rhoA, sd_rhoA) T[0, ]; // esperence of 0.1
  param[2] ~ normal(median_rhoB, sd_rhoB) T[0, ];
  param[3] ~ exponential(1);
  param[4] ~ uniform(0,2);
  param[5] ~ beta(3,1.3);
  param[6] ~ beta(3,1.3);
  
  time_c ~ uniform(0,50);
  
  
  int N_tempo_2016 = 0;
  int N_tempo_2018 = 0;

  real age_c_2016[N_2016];
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
  array[N_2016] vector[4] theta_2016 = ode_rk45(ode_Change1, initial_state, age0, age_2016, param);                    // apply ODE to all participants (will only look at those younger than decrease in LL)
  array[N_2016 - N_tempo_2016] vector[4] theta_bis_2016 = ode_rk45(ode_noChange, theta_2016[N_tempo_2016], age0, age_c_2016[(N_tempo_2016+1):N_2016], param); // aply ODE for participants older than decrease
  
  for (k in 1:N_2018){if(age_2018[k]<= (time_c+2)){N_tempo_2018 = N_tempo_2018+1;}}
  array[N_2018] vector[4] theta_2018 = ode_rk45(ode_Change1, initial_state, age0, age_2018, param);
  array[N_2018 - N_tempo_2018] vector[4] theta_bis_2018 = ode_rk45(ode_noChange, theta_2018[N_tempo_2018], age0, age_c_2018[(N_tempo_2018+1):N_2018], param);



  // Likelihood
  for(l in 1:N_tempo_2016){Statut_2016[l,] ~ multinomial(theta_2016[l,]);}
  for(j in 1:(N_2016 - N_tempo_2016)){Statut_2016[(j+N_tempo_2016),] ~ multinomial(theta_bis_2016[j,]);}
  
  for(m in 1:N_tempo_2018){Statut_2018[m,] ~ multinomial(theta_2018[m,]);}
  for(n in 1:(N_2018 - N_tempo_2018)){Statut_2018[(n+N_tempo_2018),] ~ multinomial(theta_bis_2018[n,]);}
}






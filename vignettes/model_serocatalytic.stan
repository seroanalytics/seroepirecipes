data {

  int<lower=0> AGE_GROUPS;
  vector[AGE_GROUPS] age_mid;
  int<lower=0> N[AGE_GROUPS];
  int<lower=0> n_seropos[AGE_GROUPS];
  
}

parameters {
  
  // log FOI
  real log_lambda;

}

transformed parameters {
  
  // FOI
  real<lower=0.0> lambda;
  
  lambda = exp(log_lambda);

}

model {
  
  for(a in 1:AGE_GROUPS) n_seropos[a] ~ binomial(N[a], 1.0-exp(-lambda*age_mid[a]));

}

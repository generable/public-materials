functions {
  vector SIR(real t,vector y, array[] real theta) {
    real S     = y[1];
    real I     = y[2];
    real R     = y[3];
    
    real beta  = theta[1];
    real gamma = theta[2];
    
    vector[3] dydt;
    
    dydt[1] = -beta * S * I;
    dydt[2] =  beta * S * I - gamma * I;
    dydt[3] =  gamma * I;
    
    return dydt;
  }
}

data {
  int<lower=1> n_obs;   // number of observation times
  int<lower=1> n_pop;   // total population size
  array[n_obs] int y;   // observed infected counts
  real t0;              // initial time (e.g. 0)
  array[n_obs] real ts; // times at which y was observed
}

parameters {
  array[2] real<lower=0> theta; // {beta, gamma}
  real<lower=0,upper=1> S0;     // initial susceptible fraction
}

transformed parameters {
  vector[3] y_init;             // [S(0), I(0), R(0)]
  array[n_obs] vector[3] y_hat; // solution at each ts
  array[n_obs] real lambda;     // Poisson rates

  // set initial conditions
  y_init[1] = S0;
  y_init[2] = 1 - S0;
  y_init[3] = 0;

  y_hat = ode_rk45(SIR, y_init, t0, ts, theta);

  // convert infected fraction → expected counts
  for (i in 1:n_obs)
    lambda[i] = y_hat[i, 2] * n_pop;
}

model {
  theta ~ lognormal(0, 1);
  S0    ~ beta(1, 1);
  y     ~ poisson(lambda);
}

generated quantities {
  real R_0 = theta[1] / theta[2];
}

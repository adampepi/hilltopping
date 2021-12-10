data {
  int<lower=0> N;// this is the time series length
  int<lower=0> M;//Number of time series
  int<lower=0> states[M]; // vector assigning time series to states
  int<lower=0> S; // number of states (number of true populations)
  int<lower=0> W; // wet vs dry parameters
  int<lower=0> wetdry[M]; // vector assigning wet/dry
  int<lower=0> obsVariances[M];
  int<lower=0> n_obsvar; 
  int<lower=0> n_pos; // number of non-NA values
  int<lower=0> col_indx_pos[n_pos];
  int<lower=0> row_indx_pos[n_pos];
  vector[n_pos] y; // data
  vector[N] CI [M]; // connectivity index
  vector [N] precip; //precipitation covariate
  vector [N] temp; //temp covariate
}
parameters {
  vector[S] x0; // initial states t1
  vector[S] x1; // initial states t2
  real alpha_0;
  real alpha_1;
  real alpha_2;
  real beta_precip;
  real beta_CI;
  real<lower=0> sigma_obs[n_obsvar];
}
transformed parameters {
  vector[M] pred[N];// this is the matrix of time series data?
  vector[S] x[N]; // elements accessed [N,K] -- this is the matrix of unobserved states??
  // process model
  for(s in 1:S) {
    x[1,s] = x0[s]; // initial state, vague prior below
    x[2,s] = x1[s]; // initial state, vague prior below
    for(t in 3:N) {
      x[t,s] = alpha_0 + alpha_1*x[t-1,s] + alpha_2*x[t-2,s] + beta_precip*precip[t-1] + beta_CI*CI[s,t-1]; 
      }
  }
  // map predicted states from process model to time series
  for(m in 1:M) {
    for(t in 1:N) {
      pred[t,m] = x[t,states[m]];
    }
  }
}
model {
  x0 ~ normal(0,10);
  x1 ~ normal(0,10);
   for(i in 1:n_obsvar) {
    sigma_obs[i] ~ cauchy(0,5); // observation variance
  }


  alpha_0 ~ normal(0, 10); // intercepts
  alpha_1 ~ normal(0, 10); // direct density-dependence
  alpha_2 ~ normal(0, 10); //delayed density-dependence
  beta_precip ~ normal(0, 10); //temperature
  beta_CI ~ normal(0, 10); //connectivity
  // likelihood
  for(i in 1:n_pos) {
    
   y[i] ~ normal(pred[col_indx_pos[i], row_indx_pos[i]], sigma_obs[obsVariances[row_indx_pos[i]]]);
  }
  

}
generated quantities {
  vector[n_pos] log_lik;
  // regression example in loo() package
 
 for (n in 1:n_pos) log_lik[n] = normal_lpdf(y[n] | pred[col_indx_pos[n], row_indx_pos[n]], sigma_obs[obsVariances[row_indx_pos[n]]]);

}

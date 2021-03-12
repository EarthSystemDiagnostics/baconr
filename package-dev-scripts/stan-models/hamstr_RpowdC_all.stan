// HAMStR
// 20.12.2020 Andrew Dolman
// use AR1 coefficient at each level

// scale R for each level as R^delta_c

data {
  // age control points
  int<lower=0> N; 
  vector[N] depth;
  vector[N] obs_age;
  vector[N] obs_err;

  // resolution of age-depth model
  int<lower=0> K_fine;  // number of highest resolution sections
  int<lower=0> K_tot;  // total no of gamma parameters
  
  int parent[K_tot]; // index sections to their parent sections
  
  int K_lvls; // number of hierarchical levels
  
  // level start and end indices of alpha parameters
  int lvl_SE[K_lvls, 2];
  
  int finest_idx[K_fine];
  
  // modelled depths 
  vector[K_fine] c_depth_bottom;
  vector[K_fine] c_depth_top;
  
  vector<lower = 0>[K_lvls] delta_c; // width of each section at each resolution level
  
 
  // hyperparameters for the gamma innovations
  
  // prior for the oversall mean accumulation rate
  real<lower = 0> acc_mean_prior;
  
  // shape of the gamma distributions
  real<lower = 0> shape; 
  
  // hyperparameters for prior distribution on memory strength (the AR1 coefficient)
  real<lower = 0> mem_mean;
  real<lower = 0> mem_strength;

  // observation error model parameters 
  
  int<lower=0> nu; // degrees of freedom of t error distribution
  int which_c[N]; // index observations to their fine sections
 
  //int<lower=0, upper=1> scale_R; // scale the AR1 coefficient or not
  int<lower=0, upper=1> inflate_errors; // use error inflation model or not


 // hyperparameters for the error inflation model
  real<lower = 0> infl_shape_shape;  
  real<lower = 0> infl_shape_mean;  
 
  real<lower = 0> infl_sigma_sd;  
  
 // real<lower = 0> infl_mean_shape;  
 // real<lower = 0> infl_mean_mean;  
 
}
transformed data{
  
  // transform mean and strength of memory beta distribution to alpha and beta
  // as used to parameterise beta dist in stan function
  real<lower=0> mem_alpha = mem_strength * mem_mean;
  real<lower=0> mem_beta = mem_strength * (1-mem_mean);
  
  // position of the first highest resolution innovation (alpha)
  int<lower = 1> first_K_fine = K_tot - K_fine+1;
  
 
}
parameters {
  
  // AR1 coeffiecient at 1 depth unit
  real<lower = 0, upper = 1> R;
  
  // the hierarchical gamma innovations in one long vector that will be indexed
  vector<lower = 0>[K_tot] alpha;
  
  // the age at the first modelled depth
  real age0;
  
  // the measurement error inflation factors
  // these have length 0 if inflate_errors == 0 meaning that the parameters are 
  // in scope, so the model runs, but are zero length so nothing is sampled
  real<lower = 0> infl_mean[inflate_errors];
  real<lower = 0> infl_shape_1[inflate_errors];
  vector<lower = 0>[inflate_errors ? N : 0] infl;
  //real<lower = 0> infl_sigma[inflate_errors];
}
transformed parameters{
  
 
  
  // the AR1 coefficient scaled for the thickness of the modelled sediment sections
  vector<lower = 0, upper = 1>[K_lvls] w;
  
  // the AR1 correlated innovations
  vector[K_tot] x;
  
  // the modelled ages
  vector[K_fine+1] c_ages;
  
  // the modelled ages interpolated to the positions of the data
  vector[N] Mod_age;
  
  // the inflated observation errors
  real<lower = 0> infl_shape[inflate_errors];
  vector[N] obs_err_infl;
  
  
  // AR1 coeffiecient at modelled thicknesses
  for (i in 1:K_lvls){
    real d = delta_c[i];
    w[i] = R^d;
  }
  
  if (inflate_errors == 1){
    for (n in 1:N)
    //obs_err_infl[n] = obs_err[n] + infl_sigma[1] * infl[n];
    obs_err_infl[n] = obs_err[n] + infl[n];
    infl_shape[1] = infl_shape_1[1] + 1;
  } else {
    obs_err_infl = obs_err;
  }
  
 
  // need to loop through hierarchical levels

  // scaled AR1 parameter applied to each hierarchical level
  x[1] = alpha[1];
  
  for(l in 2:K_lvls){
    x[lvl_SE[l, 1]] = alpha[lvl_SE[l, 1]];
    for(i in (lvl_SE[l, 1]+1):(lvl_SE[l, 2])){
      x[i] = w[l]*x[i-1] + (1-w[l])*alpha[i];
      }
  }
  
  
  // the cumulative sum of the highest resolution innovations
  c_ages[1] = age0;
  c_ages[2:(K_fine+1)] = age0 + cumulative_sum(x[finest_idx] * delta_c[K_lvls]);
  
  // age model interpolated to the positions of the observations
  Mod_age = c_ages[which_c] + x[which_c] .* (depth - c_depth_top[which_c]);
}

model {
  // the overall mean accumulation rate
  // weak half normal prior
  alpha[1] ~ normal(0, 10*acc_mean_prior);
  
  // the gamma distributed innovations
  // prior parameterised by use set shape and the value of its parent section
  alpha[2:K_tot] ~ gamma(shape, shape ./ x[parent[2:K_tot]]);
  
  // the memory parameters
  R ~ beta(mem_alpha, mem_beta);
  
  // the observation error inflation model
  if (inflate_errors == 1){
    //infl_mean ~ gamma(infl_mean_shape, infl_mean_shape / infl_mean_mean);
    infl_shape ~ gamma(infl_shape_shape, infl_shape_shape / infl_shape_mean);
    
    infl ~ gamma(infl_shape[1], infl_shape[1] / infl_mean[1]);
    //infl ~ gamma(infl_shape[1], infl_shape[1] / 1);
    
    infl_mean ~ normal(0, infl_sigma_sd);
  }
  
  // the Likelihood of the data given the model
  obs_age ~ student_t(nu, Mod_age, obs_err_infl);
}

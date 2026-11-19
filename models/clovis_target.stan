//////////////////////////////////////
// Stan code to model variation in Clovis tip cross-sectional area (TCSA).
// We combine individual Vail measurements with Eren's published summary from the plains (Plains data).
// This is a more developed model from Stan code provided in Otarola-Castillo and Torquato 2018, ARA
// https://www.annualreviews.org/deliver/fulltext/anthro/47/1/an47_otarola-castillo_supmat.zip?itemId=/content/suppdata/10.1146/annurev-anthro-102317-045834/1
// and Otarola-Castillo et al 2023 (https://onlinelibrary.wiley.com/doi/abs/10.1002/9781119592112.ch59; https://github.com/dreoc/the-bayesian-inferential-paradigm-in-archaeology)
//////////////////////////////////////
data {
  // Vail gives us individual TCSA measurements, in square millimeters.
  int<lower=1> N_vail;
  vector<lower=0>[N_vail] vail_tcsa;

  // Plains data are summary data providing sample size, mean, and standard deviation.
  int<lower=2> N_plains;
  real<lower=0> plains_mean_tcsa;
  real<lower=0> plains_sd_tcsa;

  // These prior settings are supplied by the calling R code.
  real mu0;
  real<lower=0> mu0_sd;
  real<lower=0> sigma_scale;
}

transformed data {
  // Put the Vail measurements on the natural-log scale.
  vector[N_vail] vail_log_tcsa = log(vail_tcsa);

  // Approximate plain's data log-scale moments under a lognormal model.
  // We use the model's lognormal assumption to derive these quantities.
  real plains_var_log = log1p(square(plains_sd_tcsa / plains_mean_tcsa));
  real plains_mean_log = log(plains_mean_tcsa) - plains_var_log / 2;
  // Keep the implied total squared deviation for the likelihood below.
  real plains_sse_log = (N_plains - 1) * plains_var_log;
}

parameters {
  // On the log scale, mu is the mean and sigma is the standard deviation.
  real mu;
  real<lower=0> sigma;
}

model {
  //////////////////////////////////////
  // The priors describe plausible values before these observations contribute. 
  // These are weakly informative priors (WIPs), we use WIPs throughout the analyses 
  // and we tested them using sensitivity analyses to gauge minimal impact on 
  // results and thus inference before any data were present.
  // Because sigma is positive, its normal prior is a half-normal.
  // However, unlike packages like rethinking, direct stan code does not provide functions for 
  // half normal distribution functions like half_normal(), so they explicitly state 
  // to use truncated versions of full distributions, see: https://mc-stan.org/docs/reference-manual/statements.html#truncated-distributions
  ////////////////////////////////////// 
  mu ~ normal(mu0, mu0_sd);
  sigma ~ normal(0, sigma_scale)T[0,];
  
  //////////////////////////////////////
  //start the likelihood using target+=
  // Each Vail measurement contributes directly to the likelihood.
  //////////////////////////////////////
  target += normal_lpdf(vail_log_tcsa | mu, sigma);

  //////////////////////////////////////
  // Add the plain's data contribution to the likelihood using the log-scale summaries derived above.
  // Both sources inform the same Clovis mean and variance on the log scale.
  //////////////////////////////////////
  target += -N_plains * log(sigma)
    - plains_sse_log / (2 * square(sigma))
    - N_plains * square(plains_mean_log - mu) / (2 * square(sigma));
}

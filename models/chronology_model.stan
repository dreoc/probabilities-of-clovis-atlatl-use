//////////////////////////////////////
// ATLATL AND BOW CHRONOLOGY
// Fit the atlatl and bow radiocarbon datasets separately.
// Estimate deposition intensity lambda and compare two  chronologies.
// H0 (ABSENT): deposition after Clovis; none during Clovis.
// H1 (PRESENT): deposition after Clovis and during Clovis.
// Account for age-dependent preservation and radiocarbon calibration.
// The deposition intensity has a Gamma(rate_shape, rate_rate) prior.

// Here we also conduct RC calibration within Stan so that the Bayesian model
// accounts for laboratory measurement error and calibration
// curve uncertainty, within stan and we do not need to call outside functions from rcarbon or bchron, etc. 
// So we input the observed RC measurements and lab errors from
// our literature sources and the IntCal20 curve text file downloaded from https://intcal.org/curves.html;
// This, among other things, makes use of stan's speedy C++ compiled code.
//////////////////////////////////////

// FUNCTIONS
functions {  
  //////////////////////////////////////
  // we use the following functions to evaluate possible calendar ages for each RC mesurement 
  // across a set of Intcal20's time intervals, using the function cell_index;
  //////////////////////////////////////
  int cell_index(real t, data vector edges) {  // Find which ordered calendar-age cell contains t.
    int left = 1;  
    int right = num_elements(edges) - 1;  
    while (left < right) {  
      int middle = (left + right + 1) / 2;  
      if (t >= edges[middle]) left = middle;  
      else right = middle - 1;  
    }
    return left;  
  }
    // We then use linear_value to interpolate the IntCal20 curve in order to calculate the expected RC measurement 
    // and calibration uncertainty at that calendar age;
  real linear_value(real t, data vector x, data vector y) {  
    int j = cell_index(t, x);  
    return y[j] + (t-x[j]) * (y[j+1]-y[j]) / (x[j+1]-x[j]); 
  }

//////////////////////////////////////
// In this block, we calculate the observed RC measurement likelihood, given that proposed calendar age 
// and both sources of error = P(RC age | calendar age, curve variance + lab variance), 
// and account for the survival probability of atlatls of that age using 
// Guyette and Perez's wood decay functions and surrovel's archaeological site preservation function,
//////////////////////////////////////
  vector age_kernel(real age, data vector cal_age, data vector cal_mean, 
                    data vector cal_sd, data vector precision,  
                    data vector mean_age, data real h, data real rate_unit,  
                    data real guyette_k, data real perez_k,  // Two wood-survival decay rates.
                    data real site_offset, data real site_exponent) {  // Surovell's site-survival offset and exponent.
    int N = num_elements(precision);  
    real mu = linear_value(age, cal_age, cal_mean);  
    real cs = linear_value(age, cal_age, cal_sd);  
    real vis = (h*exp(-guyette_k*age) + (1-h)*exp(-perez_k*age))  // Mix the two two exponential wood-survival functions.
               * pow(1+age/site_offset, -site_exponent) / rate_unit;  // Adjust by site survival and convert rate units.
    vector[N+1] out;  
    out[1] = vis;  
    for (i in 1:N) {  
      real d = 1 + square(cs)*precision[i];  
      // Radiocarbon likelihood contribution.
      out[i+1] = vis * exp(-0.5 * (log1p(square(cs)*precision[i])  // Combine visibility with marginalized calibration likelihood.
                      + square(mean_age[i]-mu)*precision[i]/d));  
    }
    return out;  
  }
  // Integrate calendar-age preservation and dating contributions.
  vector calendar_interval(real lo, real hi,
                           data vector qx, data vector qw,
                           data vector cal_age, data vector cal_mean,
                           data vector cal_sd, data vector precision,
                           data vector mean_age, data real h, data real rate_unit,
                           data real guyette_k, data real perez_k,
                           data real site_offset, data real site_exponent) {
    vector[num_elements(precision)+1] out = rep_vector(0, num_elements(precision)+1);
    if (hi > lo) {
      real half_width = (hi-lo)/2;
      for (j in 1:num_elements(qx)) {
        real age = (lo+hi)/2 + half_width*qx[j];
        real weight = half_width*qw[j];
        out += weight * age_kernel(age, cal_age, cal_mean, cal_sd,
                                  precision, mean_age, h, rate_unit,
                                  guyette_k, perez_k, site_offset, site_exponent);
      }
    }
    return out;
  }

  vector chron_log_likelihood(matrix parts, real lambda, data real shift) {
    int N = rows(parts)-1;
    vector[2] out;
    out[1] = shift - lambda*parts[1,1];
    out[2] = out[1] - lambda*parts[1,2];
    for (i in 1:N) {
      real postclovis = lambda*parts[i+1,1];
      real clovis = lambda*parts[i+1,2];
      out[1] += log(postclovis);
      out[2] += log(postclovis+clovis);
    }
    return out;
  }
}
//////////////////////////////////////
// DATA
// R brings in the dated observations, calibration grids, prior settings, 
//and preservation parameters.
//////////////////////////////////////
data {
  int<lower=1> N;
  int<lower=N> M;
  array[M] int<lower=1,upper=N> specimen_id;
  vector[M] rc_age;
  vector<lower=0>[M] rc_error;
  int<lower=2> N_curve;
  ordered[N_curve] cal_age;
  vector[N_curve] cal_mean;
  vector<lower=0>[N_curve] cal_sd;
  int<lower=2> E;
  ordered[E] edges;
  int<lower=2> E_clovis;
  ordered[E_clovis] clovis_edges;
  int<lower=2> Q_mat;
  vector[Q_mat] qx;
  vector<lower=0>[Q_mat] qw;
  real<lower=0,upper=1> h;
  real<lower=0> clovis_young;
  real<lower=clovis_young> clovis_old;
  real<lower=0> rate_unit;
  real<lower=0> guyette_k;
  real<lower=0> perez_k;
  real<lower=0> site_offset;
  real<lower=0> site_exponent;
  real<lower=0> rate_shape;
  real<lower=0> rate_rate;
  real<lower=0,upper=1> prior_absent;
}
//////////////////////////////////////
// TRANSFORMED DATA
// Precompute data.
//////////////////////////////////////
transformed data {
  vector[N] precision = rep_vector(0,N);
  vector[N] mean_age = rep_vector(0,N);
  vector[N] residual = rep_vector(0,N);
  vector[N] logvar = rep_vector(0,N);
  array[N] int n_RC_date = rep_array(0,N);
  vector[N] shift;
  vector[(E-1)*Q_mat] cached_age;
  matrix[N+1,(E-1)*Q_mat] cached_kernel;
  matrix[N+1,E] suffix = rep_matrix(0,N+1,E);
  vector[N+1] early = rep_vector(0,N+1);
  matrix[N+1,2] interval_parts;
  // Confirm calendar intervals match the two periods we compare.
  if (abs(edges[1]) > 1e-6 || abs(edges[E] - clovis_young) > 1e-6 ||
      abs(clovis_edges[1] - clovis_young) > 1e-6 ||
      abs(clovis_edges[E_clovis] - clovis_old) > 1e-6)
    reject("Calendar boundaries do not match the post-Clovis and Clovis intervals.");
  // Combine specimen RC dates.
  for (m in 1:M) {
    int i = specimen_id[m];
    real p = inv_square(rc_error[m]);
    precision[i] += p;
    mean_age[i] += p*rc_age[m];
    logvar[i] += log(square(rc_error[m]));
    n_RC_date[i] += 1;
  }
  for (i in 1:N) {
    mean_age[i] /= precision[i];
  }
  for (m in 1:M) residual[specimen_id[m]] +=
    square((rc_age[m]-mean_age[specimen_id[m]])/rc_error[m]);
  for (i in 1:N) shift[i] = -0.5*(n_RC_date[i]*log(2*pi())+logvar[i]+residual[i]);
  // Cache age integrands.
  for (c in 1:(E-1)) {
    real half_width = (edges[c+1]-edges[c])/2;
    for (k in 1:Q_mat) {
      int ix = (c-1)*Q_mat+k;
      cached_age[ix] = (edges[c]+edges[c+1])/2 + half_width*qx[k];
      cached_kernel[,ix] = half_width*qw[k] * age_kernel(cached_age[ix],
        cal_age,cal_mean,cal_sd,precision,mean_age,h,rate_unit,
        guyette_k, perez_k, site_offset, site_exponent);
    }
  }
  for (reverse_c in 1:(E-1)) {
    int c = E-reverse_c;
    suffix[,c] = suffix[,c+1];
    for (k in 1:Q_mat) suffix[,c] += cached_kernel[,(c-1)*Q_mat+k];
  }
  // Clovis period integral.
  for (c in 1:(E_clovis-1)) early += calendar_interval(
    clovis_edges[c],clovis_edges[c+1],qx,qw,
    cal_age,cal_mean,cal_sd,precision,mean_age,h,rate_unit,
    guyette_k, perez_k, site_offset, site_exponent);
  // Combine calendar-age integrals for the two periods.
  interval_parts[,1] = suffix[,1];
  interval_parts[,2] = early;
}
//////////////////////////////////////
// PARAMETERS
//////////////////////////////////////
parameters {
  real<lower=0> lambda;
}
//////////////////////////////////////
// TRANSFORMED PARAMETERS
// compute the likelihood of the chronologies at every posterior draw.
//////////////////////////////////////
transformed parameters {
  vector[2] log_likelihood;
  log_likelihood = chron_log_likelihood(
    interval_parts, lambda, sum(shift));
}
//////////////////////////////////////
// MODEL
// Combine Gamma prior and the likelihood of the chronologies.
//////////////////////////////////////
model {
  target += gamma_lpdf(lambda | rate_shape, rate_rate);

  // Marginalize over the two chronologies.
  target += log_mix(prior_absent, log_likelihood[1], log_likelihood[2]);
}
//////////////////////////////////////
// GENERATED QUANTITIES
// conditional probability of absence for each posterior draw.
//////////////////////////////////////
generated quantities {
  real p_absent = inv_logit(logit(prior_absent)
    + log_likelihood[1] - log_likelihood[2]);
}

//////////////////////////////////////
// Stan code to model joint reference distributions of log-TCSA by delivery class.
// Integrate the joint Gaussian means during fitting and recover conditional draws below.
//////////////////////////////////////

//////////////////////////////////////
// FUNCTIONS BLOCK
//////////////////////////////////////
functions {
  //////////////////////////////////////
  // gp_mean_covariance computes and returns, VCV, the variance-covariance matrix 
  // values among the means of the 15 ref data groups. They all start with their prior.
   // as the parameters are estimated by stan, the variance components are added per group.
   //////////////////////////////////////
  matrix gp_mean_covariance(array[] int class_id, array[] int dataset_id,
                             vector tau_dataset, vector tau_subgroup,
                             real mu0_sd, real tau_class) {
    int n_ref = size(class_id);
    matrix[n_ref, n_ref] VCV;
    for (g in 1:n_ref) {
      for (h in 1:n_ref) {
        // Uncertainty in the mean connects every group.
        VCV[g, h] = square(mu0_sd);
        // Class effects connect groups in the same delivery class.
        if (class_id[g] == class_id[h]) {
          VCV[g, h] += square(tau_class);
          // dataset effects also require matching dataset labels.
          if (dataset_id[g] == dataset_id[h])
            VCV[g, h] += square(tau_dataset[class_id[g]]);
        }
      }
      // Subgroup variation contributes to each group's own variance.
      VCV[g, g] += square(tau_subgroup[class_id[g]]);
    }
    return VCV;
  }
  // Computes likelihood from group means, sample sizes, and within-group variation.
  real ref_group_lik_lpdf(vector ybar, vector sse, array[] int n,
                            vector sigma_group, matrix VCV, real mu0) {
    int n_ref = num_elements(ybar);
    matrix[n_ref, n_ref] vcv_means = VCV;
    real lp = 0;
    for (g in 1:n_ref) {
      // Add the sampling variance of each group mean.
      vcv_means[g, g] += square(sigma_group[g]) / n[g];
      // Contribution from variation among points within the group.
      lp += -0.5 * (n[g] - 1) * log(2 * pi())
            - (n[g] - 1) * log(sigma_group[g]) - 0.5 * log(1.0 * n[g])
            - sse[g] / (2 * square(sigma_group[g]));
    }
    //////////////////////////////////////
    // Combine within-group variation with the joint likelihood of group means.
    // multi_normal_cholesky_lpdf calculates the multivariate normal log-likelihood using stan's 
    // Cholesky decomposition (https://mc-stan.org/docs/functions-reference/distributions_over_unbounded_vectors.html),
    // similar to R's "chol(vcv_means)".
    //////////////////////////////////////
    return lp + multi_normal_cholesky_lpdf(
      ybar | rep_vector(mu0, n_ref), cholesky_decompose(vcv_means));
  }
}
//////////////////////////////////////
// DATA BLOCK
//////////////////////////////////////
data {
  // Delivery classes, reference groups, and dataset labels.
  int<lower=2> n_deliv; // n delivery classes
  int<lower=n_deliv> n_ref; // n reference groups
  int<lower=1> ID; // ID
  array[n_ref] int<lower=1, upper=n_deliv> class_id;
  array[n_ref] int<lower=1, upper=ID> dataset_id;
  // Sample sizes, mean log-TCSA, and sums of squared deviations on that scale.
  array[n_ref] int<lower=2> n;
  vector[n_ref] ybar;
  vector<lower=0>[n_ref] sse;
  // Prior values from settings.csv.
  // mu0 and mu0_sd define the normal prior for the shared mean.
  real mu0;
  real<lower=0> mu0_sd;
  real<lower=0> tau_class_prior_scale;
  real<lower=0> sigma0_prior_scale;
  real<lower=0> tau_dataset_prior_scale;
  real<lower=0> tau_subgroup_prior_scale;
}
//////////////////////////////////////
// PARAMETER BLOCK
//////////////////////////////////////
parameters {
  // Within-group standard deviation for each delivery class.
  vector<lower=0>[n_deliv] sigma0;
  // Variation among datasets and subgroups within each class.
  vector<lower=0>[n_deliv] tau_dataset;
  vector<lower=0>[n_deliv] tau_subgroup;
  // Variation among class means around the shared mean.
  vector<lower=0>[1] tau_class;
}
//////////////////////////////////////
// TRANSFORMED PARAMETER BLOCK
//////////////////////////////////////
transformed parameters {
  // recalculates sigma0 during sampling (https://mc-stan.org/docs/reference-manual/blocks.html),
  // grouping the same classes' residual standard deviation.
  vector<lower=0>[n_ref] sigma_group;
  for (g in 1:n_ref) {
    sigma_group[g] = sigma0[class_id[g]];
  }
}
//////////////////////////////////////
// MODEL BLOCK
//////////////////////////////////////
model {
  real t = tau_class[1];
  matrix[n_ref, n_ref] VCV;
  //////////////////////////////////////
  // Half-normal priors on the residual and hierarchy standard deviations.
  // Unlike packages like rethinking, stan does not provide functions for 
  // half normal distribution functions like half_normal(), so they explicitly state 
  // to use truncated versions of full distributions, see: https://mc-stan.org/docs/reference-manual/statements.html#truncated-distributions
  //////////////////////////////////////
  sigma0 ~ normal(0, sigma0_prior_scale) T[0, ];
  tau_dataset ~ normal(0, tau_dataset_prior_scale) T[0, ];
  tau_subgroup ~ normal(0, tau_subgroup_prior_scale) T[0, ];
  tau_class ~ normal(0, tau_class_prior_scale) T[0, ];
  // Reference likelihood.
  VCV = gp_mean_covariance(class_id, dataset_id, tau_dataset, tau_subgroup,
                          mu0_sd, t);
  target += ref_group_lik_lpdf(ybar | sse, n, sigma_group, VCV, mu0);
}
//////////////////////////////////////
//GENERATED QUANTITIES BLOCK
//////////////////////////////////////
generated quantities {
  // Posterior draws of means, means for new datasets, and replicated summaries.
  vector[1] mu_shared;
  vector[n_deliv] mu_class;
  vector[n_ref] theta_group;
  vector[n_ref] mean_log_rep;
  vector[n_ref] sd_log_rep;
  vector[n_deliv] theta_new_dataset;
  vector<lower=0>[n_deliv] sigma_new_dataset;
  real log_reference;
  {
    real t = tau_class[1];
    real shared_prior;
    vector[n_deliv] mu_prior;
    matrix[n_deliv, ID] dataset_prior;
    vector[n_ref] theta_prior;
    vector[n_ref] sim_ybar;
    matrix[n_ref, n_ref] VCV;
    matrix[n_ref, n_ref] vcv_means;
    vector[n_ref] alpha;
    // Draw the Gaussian hierarchy at the current scale values.
    shared_prior = normal_rng(mu0, mu0_sd);
    VCV = gp_mean_covariance(class_id, dataset_id, tau_dataset, tau_subgroup,
                            mu0_sd, t);
    vcv_means = VCV;
    for (k in 1:n_deliv) {
      mu_prior[k] = normal_rng(shared_prior, t);
      for (s in 1:ID)
        dataset_prior[k, s] = normal_rng(0, tau_dataset[k]);
    }
    // Group means and simulated observations, include sampling variation.
    for (g in 1:n_ref) {
      int k = class_id[g];
      theta_prior[g] = mu_prior[k] + dataset_prior[k, dataset_id[g]]
                        + normal_rng(0, tau_subgroup[k]);
      sim_ybar[g] = normal_rng(theta_prior[g], sigma_group[g] / sqrt(1.0 * n[g]));
      vcv_means[g, g] += square(sigma_group[g]) / n[g];
    }
    // Condition the simulated hierarchy on the observed group means.
    alpha = mdivide_left_spd(vcv_means, ybar - sim_ybar);
    theta_group = theta_prior + VCV * alpha;
    // Recover class means and draw a mean for a new dataset and subgroup.
    for (k in 1:n_deliv) {
      real correction = 0;
      for (g in 1:n_ref) {
        real covariance = square(mu0_sd);
        if (class_id[g] == k)
          covariance += square(t);
        correction += covariance * alpha[g];
      }
      mu_class[k] = mu_prior[k] + correction;
      theta_new_dataset[k] = normal_rng(mu_class[k],
        sqrt(square(tau_dataset[k]) + square(tau_subgroup[k])));
      sigma_new_dataset[k] = sigma0[k];
    }
    // Shared mean across delivery classes.
    mu_shared[1] = shared_prior + square(mu0_sd) * sum(alpha);
    // Replicate group means and standard deviations for model checks.
    for (g in 1:n_ref) {
      mean_log_rep[g] = normal_rng(theta_group[g], sigma_group[g] / sqrt(1.0 * n[g]));
      sd_log_rep[g] = sigma_group[g]
                     * sqrt(chi_square_rng(n[g] - 1) / (n[g] - 1));
    }
    // Log likelihood with Gaussian means integrated out, conditional on sse, n, sigma_group, VCV, and mu0.
    log_reference = ref_group_lik_lpdf(ybar | sse, n, sigma_group, VCV, mu0);
  }
}

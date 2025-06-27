data {
  int<lower=0> n_obs;
  int<lower=0> n_site;
  int<lower=0> n_plot;
  vector[n_obs] ba;
  vector[n_obs] y;
  array[n_obs] int<lower=0, upper=n_site> site;
  array[n_obs] int<lower=0, upper=n_plot> plot;
}
parameters {
  vector<lower=0>[n_site] muba;
  vector[n_plot] delta_ba;
  real<lower=0> sigma_ba;
  real<lower=0> sigma_ba_p;
  vector<lower=0>[n_site] mu;
  vector[n_plot] delta;
  real<lower=0> sigma;
  real<lower=0> sigma_p;
}
model {
  ba ~ normal(muba[site] + delta_ba[plot], sigma_ba);
  delta_ba ~ normal(0, sigma_ba_p);
  y ~ normal(mu[site] + delta[plot], sigma);
  delta ~ normal(0, sigma_p);
}

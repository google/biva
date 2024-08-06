/*
# Copyright 2024 Google LLC

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     https://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
*/

data {
    // The number of observations
    int<lower=0> N;
    // The number of covariates in the Y-models
    int<lower=0> P_ymodel;
    // The number of covariates in the S-models
    int<lower=0> P_smodel;
    // The treatment assigned vector (binary instrument variable)
    int<lower=0, upper=1> Z[N];
    // The treatment received vector (binary treatment variable)
    int<lower=0, upper=1> D[N];
    // The dependent variable vector (continuous outcome variable)
    real Y[N];
    // The covariate matrix in the Y-models (including intercept after data cleaning)
    matrix[N, P_ymodel] X_ymodel;
    // The covariate matrix in the S-models (including intercept after data cleaning)
    matrix[N, P_smodel] X_smodel;
    // The assumption of exclusion restriction (ER = 1 if assumed)
    int<lower=0, upper=1> ER;
    // The number of Y-models
    int<lower=3, upper=6> K_ymodel;
    // The number of S-models (equal to the number of noncompliance sides)
    int<lower=1, upper=2> K_smodel;
    // The prior means for coefficients in the Y-models (including intercept)
    matrix[K_ymodel, P_ymodel] beta_mean_ymodel;
    // The prior standard deviations for coefficients in the Y-models (including intercept)
    matrix[K_ymodel, P_ymodel] beta_sd_ymodel;
    // The prior means for coefficients in the S-models (including intercept)
    matrix[K_smodel, P_smodel] beta_mean_smodel;
    // The prior standard deviations for coefficients in the S-models (including intercept)
    matrix[K_smodel, P_smodel] beta_sd_smodel;
    // The prior shape for standard deviation of errors in the Y-models
    real sigma_shape_ymodel[K_ymodel];
    // The prior scale for standard deviation of errors in the Y-models
    real<lower=0> sigma_scale_ymodel[K_ymodel];
    // Flag for running estimation (0: no, 1: yes)
    int<lower=0, upper=1> run_estimation;
}

transformed data {
    // S[k] indexes Y-model [k] and maps from it to the strata s (s = 1 if complier, 2 if nt, 3 if at)
    int S[K_ymodel];
    // Mzd maps units with Z=z and D=d to a vector of eligible Y-models indicated by nonzero elements
    // Mzd[i] = 0 means no eligible Y-model at the ith element
    int<lower=0, upper=K_ymodel> M00[2];
    int<lower=0, upper=K_ymodel> M01[2];
    int<lower=0, upper=K_ymodel> M10[2];
    int<lower=0, upper=K_ymodel> M11[2];
    // lengthzd is the number of eligible Y-models with Z=z and D=d
    int length00;
    int length01;
    int length10;
    int length11;

    // compliers assigned to control
    S[1] = 1;
    // compliers assigned to treatment
    S[2] = 1;
    // never-takers assigned to control (and to treatment if assume ER)
    S[3] = 2;

    if (ER == 1 && K_smodel == 1) {
        // E.g., when we assume ER and there is only one-sided noncompliance
        // There are two eligible Y-models for those assigned control and observed control
        length00 = 2;
        // Eligible Y-models: the 1st Y-model for compliers assigned to control;
        //                    the 3rd Y-model for never-takers assigned to control
        M00[1] = 1;
        M00[2] = 3;

        length01 = 0;
        M01[1] = 0;
        M01[2] = 0;

        length10 = 1;
        M10[1] = 3;
        M10[2] = 0;

        length11 = 1;
        M11[1] = 2;
        M11[2] = 0;
    }
    else if (ER == 1 && K_smodel == 2) {
        // always-takers assigned to control and treatment
        S[4] = 3;

        length00 = 2;
        M00[1] = 1;
        M00[2] = 3;

        length01 = 1;
        M01[1] = 4;
        M01[2] = 0;

        length10 = 1;
        M10[1] = 3;
        M10[2] = 0;

        length11 = 2;
        M11[1] = 2;
        M11[2] = 4;
    }
    else if (ER == 0 && K_smodel == 1) {
        // never-takers assigned to treatment
        S[4] = 2;

        length00 = 2;
        M00[1] = 1;
        M00[2] = 3;

        length01 = 0;
        M01[1] = 0;
        M01[2] = 0;

        length10 = 1;
        M10[1] = 4;
        M10[2] = 0;

        length11 = 1;
        M11[1] = 2;
        M11[2] = 0;
    }
    else if (ER == 0 && K_smodel == 2) {
        // never-takers assigned to treatment
        S[4] = 2;
        // always-takers assigned to control
        S[5] = 3;
        // always-takers assigned to treatment
        S[6] = 3;

        length00 = 2;
        M00[1] = 1;
        M00[2] = 3;

        length01 = 1;
        M01[1] = 5;
        M01[2] = 0;

        length10 = 1;
        M10[1] = 4;
        M10[2] = 0;

        length11 = 2;
        M11[1] = 2;
        M11[2] = 6;
    }
}

parameters {
    // coefficients in the Y-models (including intercept)
    matrix[K_ymodel, P_ymodel] beta_ymodel;
    // coefficients in the S-models (including intercept)
    matrix[K_smodel, P_smodel] beta_smodel;
    // standard deviations of error in the Y-models
    real<lower=0> sigma[K_ymodel];
}

model {
    // prior
    for (k in 1:K_smodel) {
    	for (p in 1:P_smodel) {
    		beta_smodel[k, p] ~ normal(beta_mean_smodel[k,p], beta_sd_smodel[k,p]);
    	}
    }
    for (k in 1:K_ymodel) {
    	for (p in 1:P_ymodel) {
    		beta_ymodel[k, p] ~ normal(beta_mean_ymodel[k,p], beta_sd_ymodel[k,p]);
    	}
    }
    for (k in 1:K_ymodel) {
    	sigma[k] ~ inv_gamma(sigma_shape_ymodel[k],sigma_scale_ymodel[k]);
    }

    if (run_estimation==1) {
      // model
      for (n in 1:N) {
        int length;
          real log_prob[K_smodel+1];
          log_prob[1] = 0;
          for (k in 2:(K_smodel+1)) {
              log_prob[k] = X_smodel[n] * beta_smodel[k-1]';
          }

          if (Z[n] == 0 && D[n] == 0) {
              length = length00;
          }
          else if (Z[n] == 1 && D[n] == 0) {
              length = length10;
          }
          else if (Z[n] == 1 && D[n] == 1) {
              length = length11;
          }
          else if (Z[n] == 0 && D[n] == 1) {
              length = length01;
          }

          {
      real log_l[length];
      if (Z[n] == 0 && D[n] == 0) {
          for (l in 1:length) {
            log_l[l] = log_prob[S[M00[l]]] + normal_lpdf(Y[n] | X_ymodel[n] * beta_ymodel[M00[l]]', sigma[M00[l]]);
          }
      }
      else if (Z[n] == 1 && D[n] == 0) {
          for (l in 1:length) {
            log_l[l] = log_prob[S[M10[l]]] + normal_lpdf(Y[n] | X_ymodel[n] * beta_ymodel[M10[l]]', sigma[M10[l]]);
          }
      }
      else if (Z[n] == 1 && D[n] == 1) {
          for (l in 1:length) {
            log_l[l] = log_prob[S[M11[l]]] + normal_lpdf(Y[n] | X_ymodel[n] * beta_ymodel[M11[l]]', sigma[M11[l]]);
          }
      }
      else if (Z[n] == 0 && D[n] == 1) {
          for (l in 1:length) {
            log_l[l] = log_prob[S[M01[l]]] + normal_lpdf(Y[n] | X_ymodel[n] * beta_ymodel[M01[l]]', sigma[M01[l]]);
          }
      }
      target += log_sum_exp(log_l) - log_sum_exp(log_prob);
          }
      }
    }
}

generated quantities {
    // the probability of being in each stratum
    vector[K_smodel+1] strata_prob;
    // mean outcome
    vector[K_ymodel] mean_outcome;
    {
        // each individual's log probability of being in each stratum
        matrix[N, K_smodel+1] log_prob;
        // each individual's expected means for Y-models
        matrix[N, K_ymodel] expected_mean;
        // weighted sum of individual's expected means for Y-models by probability of being in the relevant stratum
        vector[K_ymodel] weighted_outcome_sum;

        log_prob[:, 1] = rep_vector(0, N);
        log_prob[:, 2:(K_smodel+1)] = X_smodel * beta_smodel';
        for (n in 1:N) {
            // standardization
            log_prob[n] -= log_sum_exp(log_prob[n]);
        }

        for (n in 1:N)
            for (k in 1:K_ymodel)
                expected_mean[n, k] = X_ymodel[n] * beta_ymodel[k]';

	      // aggregate individual level predictions for the average probability and outcomes
        for (k in 1:(K_smodel+1)) {
            strata_prob[k] = mean(exp(log_prob[:, k]));
        }
        for (k in 1:K_ymodel) {
            weighted_outcome_sum[k] = mean(expected_mean[:, k] .* exp(log_prob[:, S[k]]));
            mean_outcome[k] = weighted_outcome_sum[k] / strata_prob[S[k]];
        }
    }
}

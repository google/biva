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

library(testthat)
library(BIVA)

test_that("Hierarchical BIVA works and recovers heterogeneous effects", {
  set.seed(2026)
  n_per_group <- 100
  J <- 3
  n <- n_per_group * J
  groups <- rep(paste0("Group", 1:J), each = n_per_group)
  X <- rnorm(n)
  true.PS <- rbinom(n, 1, 0.5)
  Z <- rep(c(0, 1), n / 2)
  D <- rep(0, n)
  c.trt.ind <- (true.PS == 1) & (Z == 1)
  c.ctrl.ind <- (true.PS == 1) & (Z == 0)
  nt.ind <- (true.PS == 0)
  D[c.trt.ind] <- 1
  
  # Ground Truths: G1=0, G2=2, G3=4
  Y <- rep(0, n)
  Y[nt.ind] <- rnorm(sum(nt.ind), 5.0, 0.1)
  
  # Group 1
  idx1 <- (groups == "Group1")
  Y[idx1 & c.ctrl.ind] <- rnorm(sum(idx1 & c.ctrl.ind), 5.0, 0.1)
  Y[idx1 & c.trt.ind]  <- rnorm(sum(idx1 & c.trt.ind), 5.0, 0.1)
  
  # Group 2
  idx2 <- (groups == "Group2")
  Y[idx2 & c.ctrl.ind] <- rnorm(sum(idx2 & c.ctrl.ind), 5.0, 0.1)
  Y[idx2 & c.trt.ind]  <- rnorm(sum(idx2 & c.trt.ind), 7.0, 0.1)
  
  # Group 3
  idx3 <- (groups == "Group3")
  Y[idx3 & c.ctrl.ind] <- rnorm(sum(idx3 & c.ctrl.ind), 5.0, 0.1)
  Y[idx3 & c.trt.ind]  <- rnorm(sum(idx3 & c.trt.ind), 9.0, 0.1)
  
  df <- data.frame(Y = Y, Z = Z, D = D, X = X, G = groups)

  iv_obj <- biva$new(
    data = df, y = "Y", d = "D", z = "Z", group = "G",
    x_ymodel = "X", x_smodel = "X", side = 1, ER = 1,
    beta_mean_ymodel = matrix(5, 3, 2), 
    beta_sd_ymodel = matrix(5, 3, 2),
    beta_mean_smodel = matrix(0, 1, 2),
    beta_sd_smodel = matrix(2, 1, 2),
    sigma_shape_ymodel = rep(1, 3),
    sigma_scale_ymodel = rep(1, 3),
    chains = 1, iter = 500, seed = 2026
  )
  
  summary_df <- iv_obj$groupSummary()
  expect_equal(nrow(summary_df), J)
  expect_true(summary_df$point_estimate[1] < 1.0)
  expect_true(summary_df$point_estimate[2] > 1.0 && summary_df$point_estimate[2] < 3.0)
  expect_true(summary_df$point_estimate[3] > 3.0)
  
  # Test predict with grouping
  pred_df <- df[1:5, ]
  iv_obj$predict(new_data = pred_df, name = "hier_pred")
  # Check if predictions was stored
  expect_true("hier_pred" %in% iv_obj$predict_list)
})

test_that("groupSummary returns warning when not hierarchical", {
  set.seed(1997)
  n <- 50
  Z <- rbinom(n, 1, 0.5)
  D <- Z * rbinom(n, 1, 0.5) # One-sided noncompliance: D=1 only if Z=1
  df <- data.frame(Y = rnorm(n), Z = Z, D = D, X = rnorm(n))
  iv_obj <- biva$new(
    data = df, y = "Y", d = "D", z = "Z",
    x_ymodel = "X", x_smodel = "X", side = 1, ER = 1,
    beta_mean_ymodel = matrix(0, 3, 2), beta_sd_ymodel = matrix(1, 3, 2),
    beta_mean_smodel = matrix(0, 1, 2), beta_sd_smodel = matrix(1, 1, 2),
    sigma_shape_ymodel = rep(1, 3), sigma_scale_ymodel = rep(1, 3),
    chains = 1, iter = 200, fit = TRUE
  )
  expect_warning(res <- iv_obj$groupSummary(), "The model was not fit with groups")
  expect_null(res)
})

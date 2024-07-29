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
test_that("The Bayesian Instrumental Variable Analysis Works", {
  set.seed(1997)
  n <- 200
  X <- rnorm(n)
  true.PS <- rbinom(n, 1, 0.4)
  Z <- c(rep(0, n / 2), rep(1, n / 2))
  D <- rep(0, n)
  c.trt.ind <- (true.PS == 1) & (Z == 1)
  c.ctrl.ind <- (true.PS == 1) & (Z == 0)
  nt.ind <- (true.PS == 0)
  num.c.trt <- sum(c.trt.ind)
  num.c.ctrl <- sum(c.ctrl.ind)
  num.nt <- sum(nt.ind)
  D[c.trt.ind] <- rep(1, num.c.trt)
  Y <- rep(0, n)
  Y[c.ctrl.ind] <- rnorm(num.c.ctrl, 0.8, 0.1)
  Y[c.trt.ind] <- rnorm(num.c.trt, 0.4, 0.2)
  Y[nt.ind] <- rnorm(num.nt, 1.3, 0.3)
  df <- data.frame(Y = Y, Z = Z, D = D, X = X)

  iv_fake <- biva$new(
    data = df, y = "Y", d = "D", z = "Z",
    x_ymodel = c("X"),
    x_smodel = c("X"),
    ER = 1,
    side = 1,
    beta_mean_ymodel = matrix(0, 3, 2),
    beta_sd_ymodel = matrix(1, 3, 2),
    beta_mean_smodel = matrix(0, 1, 2),
    beta_sd_smodel = matrix(1, 1, 2),
    sigma_shape_ymodel = rep(1, 3),
    sigma_scale_ymodel = rep(1, 3),
    fit = TRUE
  )
  expect_s3_class(iv_fake, "biva")
})

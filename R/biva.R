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

#' @title Bayesian Instrument Variable Analysis
#' @docType class
#' @export
#' @field version BIVA package version used to fit model

biva <- R6::R6Class(
  classname = "biva",
  private = list(
    ..version = NULL,
    ..stan_data = NULL,
    ..y_type = NULL,
    ..ER = NULL,
    ..side = NULL,
    ..x_ymodel = NULL,
    ..x_smodel = NULL,
    ..name_ymodel = NULL,
    ..stanfit = NULL,
    ..mcmc_checks = NULL,
    ..weak_IV = NULL,
    ..beta_ymodel = NULL,
    ..beta_smodel = NULL,
    ..sigma = NULL,
    ..strata_prob = NULL,
    ..mean_outcome = NULL,
    ..CACE_draws = NULL,
    ..credible_interval = NULL,
    ..predict_list = NULL,
    ..predictions_s = NULL,
    ..predictions_y = NULL,
    ..prior_mean_outcome = NULL,
    ..prior_strata_prob = NULL,
    ..prior_CACE = NULL
  ),
  active = list(
    version = function() {
      return(private$..version)
    }
  ),
  public = list(
    #' @description
    #' Create a new Bayesian Instrumental Variable Analysis object.
    #'
    #' @param data A data.frame containing the data to be cleaned.
    #' @param y Name of the dependent variable (character; numeric).
    #' @param d Name of the treatment received variable (character; numeric 0 or 1).
    #' @param z Name of the treatment assigned variable (character; numeric 0 or 1).
    #' @param x_ymodel Names of the covariates to include in the Y-model (character vector, optional).
    #' @param x_smodel Names of the covariates to include in the S-model (character vector, optional).
    #' @param y_type The type of the dependent variable ("real" or "binary")
    #' @param ER The assumption of exclusion restriction (ER = 1 if assumed, 0 otherwise)
    #' @param side The number of noncompliance sides (1 or 2)
    #' @param beta_mean_ymodel The prior means for coefficients in the Y-models including intercept (numeric matrix; dimension = # outcome models x # covariates)
    #' @param beta_sd_ymodel The prior standard deviations for coefficients in the Y-models including intercept (numeric matrix; dimension = # outcome models x # covariates)
    #' @param beta_mean_smodel The prior means for coefficients in the S-models including intercept (numeric matrix; dimension = (# strata - 1) x # covariates)
    #' @param beta_sd_smodel The prior standard deviations for coefficients in the S-models including intercept (numeric matrix; dimension = (# strata - 1) x # covariates)
    #' @param sigma_shape_ymodel The prior shape for standard deviation of errors in the Y-models (numeric vector; length = # outcome models)
    #' @param sigma_scale_ymodel The prior scale for standard deviation of errors in the Y-models (numeric vector; length = # outcome models)
    #' @param seed Seed for Stan fitting
    #' @param fit Flag for fitting the data to the model or not (1 if fit, 0 otherwise)
    #' @param ... Additional arguments for Stan
    #' @return invisible
    initialize = function(data, y, d, z,
                          x_ymodel = NULL,
                          x_smodel = NULL,
                          y_type = "real",
                          ER = 1,
                          side = 2,
                          beta_mean_ymodel,
                          beta_sd_ymodel,
                          beta_mean_smodel,
                          beta_sd_smodel,
                          sigma_shape_ymodel,
                          sigma_scale_ymodel,
                          seed = 1997,
                          fit = TRUE,
                          ...) {
      private$..version <- packageVersion("biva")
      private$..y_type <- y_type
      private$..ER <- ER
      private$..side <- side
      private$..x_ymodel <- x_ymodel
      private$..x_smodel <- x_smodel
      private$..name_ymodel <- c()
      if (ER == 1 & side == 1) {
        private$..name_ymodel <- c(
          "Complier assigned ctrl",
          "Complier assigned trt",
          "Nevertaker"
        )
      } else if (ER == 1 & side == 2) {
        private$..name_ymodel <- c(
          "Compliers assigned ctrl",
          "Compliers assigned trt",
          "Nevertaker",
          "Alwaystaker"
        )
      } else if (ER == 0 & side == 1) {
        private$..name_ymodel <- c(
          "Compliers assigned ctrl",
          "Compliers assigned trt",
          "Nevertaker assigned ctrl",
          "Nevertaker assigned trt"
        )
      } else if (ER == 0 & side == 2) {
        private$..name_ymodel <- c(
          "Compliers assigned ctrl",
          "Compliers assigned trt",
          "Nevertaker assigned ctrl",
          "Nevertaker assigned trt",
          "Alwaystaker assigned ctrl",
          "Alwaystaker assigned trt"
        )
      }

      cleaned_data <- CleanData(
        data = data,
        y = y,
        d = d,
        z = z,
        x_ymodel = x_ymodel,
        x_smodel = x_smodel,
        ER = ER,
        side = side
      )
      if (y_type == "real") {
        stan_data <- list(
          N = cleaned_data$N,
          P_ymodel = cleaned_data$P_ymodel,
          P_smodel = cleaned_data$P_smodel,
          Z = cleaned_data$Z,
          D = cleaned_data$D,
          Y = cleaned_data$Y,
          X_ymodel = cleaned_data$X_ymodel,
          X_smodel = cleaned_data$X_smodel,
          ER = ER,
          K_ymodel = cleaned_data$K_ymodel,
          K_smodel = cleaned_data$K_smodel,
          beta_mean_ymodel = beta_mean_ymodel,
          beta_sd_ymodel = beta_sd_ymodel,
          beta_mean_smodel = beta_mean_smodel,
          beta_sd_smodel = beta_sd_smodel,
          sigma_shape_ymodel = sigma_shape_ymodel,
          sigma_scale_ymodel = sigma_scale_ymodel,
          run_estimation = 0
        )
      } else if (y_type == "binary") {
        stan_data <- list(
          N = cleaned_data$N,
          P_ymodel = cleaned_data$P_ymodel,
          P_smodel = cleaned_data$P_smodel,
          Z = cleaned_data$Z,
          D = cleaned_data$D,
          Y = cleaned_data$Y,
          X_ymodel = cleaned_data$X_ymodel,
          X_smodel = cleaned_data$X_smodel,
          ER = ER,
          K_ymodel = cleaned_data$K_ymodel,
          K_smodel = cleaned_data$K_smodel,
          beta_mean_ymodel = beta_mean_ymodel,
          beta_sd_ymodel = beta_sd_ymodel,
          beta_mean_smodel = beta_mean_smodel,
          beta_sd_smodel = beta_sd_smodel,
          run_estimation = 0
        )
      }

      private$..stan_data <- stan_data
      # Draw from the prior
      if (y_type == "real") {
        sim_out <- rstan::sampling(stanmodels$Gaussian,
          data = private$..stan_data
        )
      } else if (y_type == "binary") {
        sim_out <- rstan::sampling(stanmodels$logit,
          data = private$..stan_data
        )
      }
      stan_data$run_estimation <- 1
      private$..prior_mean_outcome <-
        rstan::extract(sim_out)$"mean_outcome"
      colnames(private$..prior_mean_outcome) <- private$..name_ymodel
      private$..prior_strata_prob <-
        rstan::extract(sim_out)$"strata_prob"
      private$..prior_CACE <-
        private$..prior_mean_outcome[, 2] - private$..prior_mean_outcome[, 1]

      if (fit) {
        message("Fitting model to the data")
        if (y_type == "real") {
          private$..stanfit <- rstan::sampling(stanmodels$Gaussian,
            data = stan_data, ...
          )
        } else if (y_type == "binary") {
          private$..stanfit <- rstan::sampling(stanmodels$logit,
            data = stan_data, ...
          )
        }

        private$..mcmc_checks <- im::mcmcChecks$new(
          fit = private$..stanfit,
          pars = "mean_outcome"
        )

        private$..beta_ymodel <-
          rstan::extract(private$..stanfit)$"beta_ymodel"

        private$..beta_smodel <-
          rstan::extract(private$..stanfit)$"beta_smodel"

        if (y_type == "real") {
          private$..sigma <-
            rstan::extract(private$..stanfit)$"sigma"
        }

        private$..strata_prob <-
          rstan::extract(private$..stanfit)$"strata_prob"

        private$..mean_outcome <-
          rstan::extract(private$..stanfit)$"mean_outcome"
        colnames(private$..mean_outcome) <- private$..name_ymodel

        private$..CACE_draws <-
          private$..mean_outcome[, 2] - private$..mean_outcome[, 1]
      }

      return(invisible())
    },

    #' @description
    #' Plot MCMC trace for the mean outcomes in each strata.
    #' @param ... Additional arguments for Stan
    #' @return A ggplot object.
    tracePlot = function(...) {
      return(
        bayesplot::mcmc_trace(private$..mean_outcome, ...) +
          ggplot2::labs(title = "MCMC Trace for Mean Outcomes by Strata")
      )
    },

    #' @description
    #' Perform a weak instrument test.
    #'
    #' @return A character string summarizing the test result and warnings
    #' if the instrument is tested weak.
    weakIVTest = function() {
      prob_estimate_c_digit <- colMeans(private$..strata_prob)[1]
      prob_estimate_c <- scales::percent(prob_estimate_c_digit)
      txt <- "The instrument is"
      if (prob_estimate_c_digit < 0.1) {
        private$..weak_IV <- TRUE
        statement <- glue::glue(
          "{txt} considered weak, resulting in an estimated ",
          "{prob_estimate_c} of compliers in the population. ",
          "Posterior point estimates of the CACE might be unreliable."
        )
      } else {
        private$..weak_IV <- FALSE
        statement <- glue::glue(
          "{txt} not considered weak, resulting in an estimated ",
          "{prob_estimate_c} of compliers in the population."
        )
      }
      return(statement)
    },

    #' @description
    #' Calculates the probability of a unit being in each of the assumed strata.
    #'
    #' @return  A character string summarizing the estimated probability
    #'
    #'
    strataProbEstimate = function() {
      prob_estimate_c <- scales::percent(colMeans(private$..strata_prob)[1])
      prob_estimate_nt <- scales::percent(colMeans(private$..strata_prob)[2])
      txt <- "Given the data, we estimate"
      if (private$..side == 2) {
        prob_estimate_at <- scales::percent(colMeans(private$..strata_prob)[3])
        statement <- glue::glue(
          "{txt} that there is ",
          "a {prob_estimate_c} probability that the unit is a complier, ",
          "a {prob_estimate_nt} probability that the unit is a never-taker, ",
          "and a {prob_estimate_at} probability that the unit is an always-taker."
        )
      } else {
        statement <- glue::glue(
          "{txt} that there is ",
          "a {prob_estimate_c} probability that the unit is a complier, ",
          "and a {prob_estimate_nt} probability that the unit is a never-taker."
        )
      }
      return(statement)
    },

    #' @description
    #' Calculates the posterior probability of
    #' the complier average causal effect (CACE)
    #' being greater than, less than,
    #' or within a range defined by thresholds.
    #'
    #' @param a Optional. Lower bound for the threshold.
    #' @param b Optional. Upper bound for the threshold.
    #' @param prior Logical. If TRUE, calculates probabilities based on
    #' the prior distribution.
    #'        If FALSE (default), uses the posterior distribution.
    #'
    #' @return  A character string summarizing the estimated probability
    #'
    #'
    calcProb = function(a = 0, b = NULL, prior = FALSE) {
      # Input validation (same as before)
      if (is.null(a) && is.null(b)) {
        stop("Either 'a' or 'b' must be provided.")
      }
      if (!is.null(a) && !is.null(b) && b <= a) {
        stop("'b' must be greater than 'a'.")
      }
      if (prior) {
        CACE_draws <- private$..prior_CACE
        txt <- "Our prior is"
      } else {
        CACE_draws <- private$..CACE_draws
        txt <- "Given the data, we estimate"
      }
      if (!is.null(a) && is.null(b)) {
        p <- scales::percent(mean(CACE_draws > a))
        statement <- glue::glue(
          "{txt} that the ",
          "probability that the effect is more than {a}",
          " is {p}."
        )
      } else if (is.null(a) && !is.null(b)) {
        p <- mean(CACE_draws < b)
        statement <- glue::glue(
          "{txt} that the probability that the",
          " effect is less than {b} is {p}."
        )
      } else { # both 'a' and 'b' are present
        p <- mean(CACE_draws > a & CACE_draws < b)
        statement <- glue::glue(
          "{txt} that the probability that the effect",
          " is between {a} and {b} is {p}."
        )
      }
      return(statement)
    },

    #' @description
    #' Calculates point estimate of the complier average causal effect (CACE).
    #'
    #' This R6 method calculates the point estimate of the effect size
    #' based on the posterior draws of the CACE parameter.
    #'
    #' @param median Logical value. If TRUE (default), the median of
    #'     the CACE draws is returned. If FALSE, the mean is returned.
    #'
    #' @return A numeric value representing the point estimate.
    #'
    #' @details This method uses the private$..mean_outcome internal variable
    #'     which contains MCMC draws of the mean outcome parameters in each
    #'     strata, based on which the CACE parameter is calculated.
    #'     Depending on the specified median argument, the method
    #'     calculates and returns either the median or the mean of the draws.
    pointEstimate = function(median = TRUE) {
      if (median) {
        return(median(private$..CACE_draws))
      } else {
        return(mean(private$..CACE_draws))
      }
    },

    #' @description
    #' Calculates credible interval for the complier average causal effect (CACE).
    #'
    #' This R6 method calculates and returns a formatted statement summarizing
    #' the credible interval of a specified width for the effect of the intervention.
    #'
    #' @param width Numeric value between 0 and 1 representing the desired
    #'     width of the credible interval (e.g., 0.95 for a 95% credible interval).
    #' @param round Integer value indicating the number of decimal places to round
    #'     the lower and upper bounds of the credible interval.
    #'
    #' @return A character string with the following information:
    #'   - The probability associated with the specified width
    #'   - The lower and upper bounds of the credible interval, rounded to the
    #'      specified number of decimal places
    #'
    #' @details This method uses the private$..mean_outcome internal variable
    #'     which contains MCMC draws of the mean outcome parameters in each
    #'     strata, based on which the CACE parameter is calculated.
    #'     It calculates the credible interval, stores it internally, and
    #'     returns a formatted statement summarizing the findings.
    #'
    credibleInterval = function(width = 0.75, round = 2) {
      private$..credible_interval <- im::credibleInterval(
        draws = private$..CACE_draws, width
      )
      statement <- glue::glue(
        "Given the data, we estimate that there is a ",
        "{scales::percent(width)} probability that the CACE is between ",
        "{round(private$..credible_interval$lower_bound, round)} and ",
        "{round(private$..credible_interval$upper_bound, round)}."
      )
      return(statement)
    },

    #' @description
    #' Plots the prior and posterior distributions
    #' for the complier average causal effect (CACE).
    #'
    #' For more details see [vizdraws::vizdraws()].
    #' @param ... other arguments passed to vizdraws.
    #' @return An interactive plot of the prior and posterior distributions.
    vizdraws = function(...) {
      p <- vizdraws::vizdraws(
        prior = private$..prior_CACE,
        posterior = private$..CACE_draws, ...
      )
      return(p)
    },

    #' @description
    #' Plots lollipop chart for the prior and posterior of
    #' the complier average causal effect (CACE) being
    #' greater or less than a threshold.
    #'
    #' For more details see [vizdraws::lollipops()].
    #' @param threshold cutoff used to calculate the probability. Defaults to
    #' zero
    #' @param ... other arguments passed to vizdraws.
    #' @return A lollipop chart with the prior and posterior probability of
    #' the CACE being above or below a threshold.
    lollipop = function(threshold = 0, ...) {
      data <- data.frame(
        Name = "CACE",
        Prior = mean(private$..prior_CACE > threshold),
        Posterior = mean(private$..CACE_draws > threshold)
      )
      p <- vizdraws::lollipops(data, ...)
      return(p)
    },

    #' @description
    #' Plots draws from the prior distribution of probability of
    #' a unit being in each of the assumed strata, and the
    #' complier average causal effect (CACE).
    plotPrior = function() {
      CACE <- ggplot2::ggplot(
        data = tibble::tibble(draws = private$..prior_CACE),
        ggplot2::aes(x = draws)
      ) +
        ggplot2::geom_histogram(bins = 30) +
        ggplot2::annotate("text",
          x = quantile(private$..prior_CACE, prob = 0.9),
          y = 200, label = paste0(
            "mean: ",
            round(mean(private$..prior_CACE), 1)
          ),
          hjust = 0
        ) +
        ggplot2::xlab("Complier average causal effect (CACE)") +
        ggplot2::ylab("N draws") +
        ggplot2::theme_minimal()

      c_prob <- ggplot2::ggplot(
        data = tibble::tibble(draws = private$..prior_strata_prob[, 1]),
        ggplot2::aes(x = draws)
      ) +
        ggplot2::geom_histogram(bins = 30) +
        ggplot2::scale_x_continuous(labels = scales::percent) +
        ggplot2::annotate("text",
          x = quantile(private$..prior_strata_prob[, 1], prob = 0.8),
          y = 200, label = paste0(
            "mean: ",
            round(mean(private$..prior_strata_prob[, 1] * 100), 1),
            "%"
          ),
          hjust = 0
        ) +
        ggplot2::xlab("proportion of compliers") +
        ggplot2::ylab("N draws") +
        ggplot2::theme_minimal()

      nt_prob <- ggplot2::ggplot(
        data = tibble::tibble(draws = private$..prior_strata_prob[, 2]),
        ggplot2::aes(x = draws)
      ) +
        ggplot2::geom_histogram(bins = 30) +
        ggplot2::scale_x_continuous(labels = scales::percent) +
        ggplot2::annotate("text",
          x = quantile(private$..prior_strata_prob[, 2], prob = 0.8),
          y = 200, label = paste0(
            "mean: ",
            round(mean(private$..prior_strata_prob[, 2] * 100), 1),
            "%"
          ),
          hjust = 0
        ) +
        ggplot2::xlab("proportion of never-takers") +
        ggplot2::ylab("N draws") +
        ggplot2::theme_minimal()

      if (private$..side == 2) {
        at_prob <- ggplot2::ggplot(
          data = tibble::tibble(draws = private$..prior_strata_prob[, 3]),
          ggplot2::aes(x = draws)
        ) +
          ggplot2::geom_histogram(bins = 30) +
          ggplot2::scale_x_continuous(labels = scales::percent) +
          ggplot2::annotate("text",
            x = quantile(private$..prior_strata_prob[, 3], prob = 0.8),
            y = 200, label = paste0(
              "mean: ",
              round(mean(private$..prior_strata_prob[, 3] * 100), 1),
              "%"
            ),
            hjust = 0
          ) +
          ggplot2::xlab("proportion of always-takers") +
          ggplot2::ylab("N draws") +
          ggplot2::theme_minimal()

        plots <- ggpubr::ggarrange(CACE,
          ggpubr::ggarrange(c_prob, nt_prob, at_prob,
            labels = c("c_prob", "nt_prob", "at_prob"), ncol = 3
          ),
          labels = "CACE", nrow = 2
        )
      } else {
        plots <- ggpubr::ggarrange(CACE,
          ggpubr::ggarrange(c_prob, nt_prob,
            labels = c("c_prob", "nt_prob"), ncol = 2
          ),
          labels = "CACE", nrow = 2
        )
      }

      plots <- ggpubr::annotate_figure(plots,
        top = ggpubr::text_grob("Draws from prior distributions",
          face = "bold", size = 14
        )
      )
      return(plots)
    },

    #' @description
    #' Predict new data

    #' @param new_data Data frame to be predicted
    #' @param name Group name of the prediction
    #' @param M Number of posterior draws to sample from
    predict = function(new_data, name = NULL, M = NULL, ...) {
      # Check new_data contains correct columns
      cols <- c(private$..x_smodel, private$..x_ymodel)
      if (!all(cols %in% colnames(new_data))) {
        missing_cols <- cols[which(!(cols %in% colnames(new_data)))]
        statement <- glue::glue(
          "{missing_cols} is missing."
        )
        stop(statement)
      }

      # Check number of posterior draws is not out of range
      if (is.null(M)) {
        M <- length(private$..CACE_draws)
      }
      if (M > length(private$..CACE_draws)) {
        M <- length(private$..CACE_draws)
        warning(
          "Number of posterior draws can not exceed ", M, ". ",
          "Setting number of posterior draws to ", M, ". "
        )
      }

      # Create prediction list if NULL
      if (is.null(private$..predictions)) {
        private$..predictions_s <- list()
        private$..predictions_y <- list()
        private$..predict_list <- list()
      }
      # Create group name if NULL
      if (is.null(name)) {
        name <- paste0("pred", length(private$..predictions) + 1)
        warning(
          "No name was supplied, assigning predictions to ", name, "."
        )
      }

      # Transform data
      X_spred <- new_data[, private$..x_smodel]
      X_ypred <- new_data[, private$..x_ymodel]

      # Sample posterior predictives for the S-model
      N_s <- nrow(X_spred)
      s_sim <- purrr::pmap(
        .l = list(
          beta_smodel = purrr::array_branch(private$..beta_smodel, 1)
        ),
        .f = function(beta_smodel, X, N) {
          exp_lin_pred <- exp(cbind(1, as.matrix(X)) %*% t(as.matrix(beta_smodel)))
          exp_lin_pred <- cbind(1, exp_lin_pred)
          return(exp_lin_pred / rowSums(exp_lin_pred))
        },
        X = X_spred, N = N_s
      )

      # Sample posterior predictives for the Y-model
      N_y <- nrow(X_ypred)
      if (private$..y_type == "real") {
        y_sim <- purrr::pmap(
          .l = list(
            beta_ymodel = purrr::array_branch(private$..beta_ymodel, 1)
          ),
          .f = function(beta_ymodel, X, N) {
            lin_pred <- cbind(1, as.matrix(X)) %*% t(as.matrix(beta_ymodel))
            return(lin_pred)
          },
          X = X_spred, N = N_s
        )
      } else if (private$..y_type == "binary") {
        y_sim <- purrr::pmap(
          .l = list(
            beta_ymodel = purrr::array_branch(private$..beta_ymodel, 1)
          ),
          .f = function(beta_ymodel, X, N) {
            lin_pred <- cbind(1, as.matrix(X)) %*% t(as.matrix(beta_ymodel))
            prob <- 1 / (1 + exp(-lin_pred))
            return(prob)
          },
          X = X_spred, N = N_s
        )
      }

      private$..predict_list[[length(private$..predict_list) + 1]] <- name
      private$..predictions_s[[name]] <- s_sim
      private$..predictions_y[[name]] <- y_sim
      return(invisible(self))
    }
  )
)

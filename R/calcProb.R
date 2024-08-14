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
#' Calculate Probability of Posterior Draws Falling Within a Range
#'
#' This function estimates the probability that a vector of posterior draws,
#' represented by the parameter `x`, falls within a specified range.
#' It provides flexibility to use either prior distributions or posterior draws,
#' and to specify one-sided or two-sided probability calculations.
#'
#' @param x A numeric vector containing either posterior draws (default)
#'          or prior samples of the treatment effect parameter.
#' @param a (Optional) The lower bound of the range.
#' @param b (Optional) The upper bound of the range.
#' @param prior A logical value indicating whether to use prior samples (`TRUE`)
#'              or posterior draws (`FALSE`, default) for calculation.
#' @param group_name A string describing the group for which the probability
#'                   is being calculated (default: "group average").
#'
#' @return A formatted string stating the calculated probability and the
#'         specified range. The probability is the proportion of samples
#'         (either prior or posterior) that fall within the defined range.
#'
#' @details This function checks the following cases:
#'  * If both `a` and `b` are `NULL`, it returns an empty string.
#'  * If `b` is less than or equal to `a`, it throws an error.
#'
#' The calculated probability and range are presented in a human-readable string
#' using the `glue` package for formatting.

calcProb <- function(
    x,
    a = NULL,
    b = NULL,
    prior = FALSE,
    group_name = "group average") {
  # Input validation (same as before)
  if (is.null(a) && is.null(b)) {
    statement <- ""
  }
  if (!is.null(a) && !is.null(b) && b <= a) {
    stop("'b' must be greater than 'a'.")
  }
  if (prior) {
    txt <- "Our prior is"
  } else {
    txt <- "Given the data, we estimate"
  }
  if (!is.null(a) && is.null(b)) {
    p <- scales::percent(mean(x > a))
    statement <- glue::glue(
      "{txt} that the ",
      "probability that the {group_name} is more than {a} ",
      "is {p}."
    )
  } else if (is.null(a) && !is.null(b)) {
    p <- scales::percent(mean(x < b))
    statement <- glue::glue(
      "{txt} that the ",
      "probability that the {group_name} is less than {b} is {p}."
    )
  } else { # both 'a' and 'b' are present
    p <- scales::percent(mean(x > a & x < b))
    statement <- glue::glue(
      "{txt} that the probability that the {group_name} ",
      "is between {a} and {b} is {p}."
    )
  }
  return(statement)
}

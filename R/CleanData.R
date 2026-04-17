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

#' Cleans and prepares data for analysis
#'
#' This function performs a series of data cleaning and preprocessing steps
#' to ensure the data is suitable for analysis.
#'
#' @param data A data.frame containing the data to be cleaned.
#' @param y Name of the dependent variable (character; numeric).
#' @param d Name of the treatment received variable (character; numeric 0 or 1).
#' @param z Name of the treatment assigned variable (character; numeric 0 or 1).
#' @param x_ymodel Names of the covariates to include in the Y-model (character vector, optional).
#' @param x_smodel Names of the covariates to include in the S-model (character vector, optional).
#' @param ER The assumption of exclusion restriction (numeric; 1 if assumed, 0 otherwise).
#' @param side The number of noncompliance sides (numeric; 1 or 2).
#' @param group Name of the grouping variable for hierarchical modeling (character, optional).
#'
#' @return A list containing the cleaned dataset and relevant metadata:
#'  * \code{N}: The number of observations after cleaning.
#'  * \code{P_ymodel} The number of covariates in the Y-model after cleaning.
#'  * \code{P_smodel} The number of covariates in the S-model after cleaning.
#'  * \code{Z} The treatment assigned vector (instrument variable).
#'  * \code{D} The treatment received vector (treatment variable).
#'  * \code{Y} The dependent variable vector (outcome variable).
#'  * \code{X_ymodel} The cleaned covariate matrix in the Y-model (including intercept).
#'  * \code{X_smodel} The cleaned covariate matrix in the S-model (including intercept).
#'  * \code{group_idx} The numeric index of the group for each unit (1 to J).
#'  * \code{J} The total number of unique groups.
#'  * \code{group_levels} The original levels of the grouping variable.
#'  * \code{use_hierarchical} Flag (0 or 1) indicating if hierarchical modeling is enabled.
#'  * \code{K_ymodel} The number of Y-models.
#'  * \code{K_smodel} The number of S-models.
#' @export

CleanData <- function(data,
                      y,
                      d,
                      z,
                      x_ymodel = NULL,
                      x_smodel = NULL,
                      ER = ER,
                      side = side) {
  stopifnot(is.data.frame(data))

  Y <- dplyr::pull(data, y)
  D <- dplyr::pull(data, d)
  Z <- dplyr::pull(data, z)
  X_ymodel <- dplyr::select(data, tidyselect::all_of(x_ymodel))
  X_smodel <- dplyr::select(data, tidyselect::all_of(x_smodel))

  if (!is.null(group)) {
    group_vec <- dplyr::pull(data, group)
    group_factor <- as.factor(group_vec)
    group_idx <- as.numeric(group_factor)
    J <- length(unique(group_idx))
    group_levels <- levels(group_factor)
    use_hierarchical <- 1
  } else {
    group_idx <- rep(1, nrow(data))
    J <- 1
    group_levels <- NULL
    use_hierarchical <- 0
  }

  if (!is.null(x_ymodel)) {
    xs_ymodel <- paste(x_ymodel, collapse = "+")
    # Build the model formula
    fml_ymodel <- paste0(y, " ~ 1 + ", xs_ymodel)
  } else {
    # Build the model formula
    fml_ymodel <- paste0(y, " ~ 1")
  }
  # Extract model matrix and keep intercept
  X_ymodel <- stats::model.matrix(stats::as.formula(fml_ymodel),
                                  data = data)

  if (!is.null(x_smodel)) {
    xs_smodel <- paste(x_smodel, collapse = "+")
    # Build the model formula
    fml_smodel <- paste0(y, " ~ 1 + ", xs_smodel)
  } else {
    # Build the model formula
    fml_smodel <- paste0(y, " ~ 1")
  }
  # Extract model matrix and keep intercept
  X_smodel <- stats::model.matrix(stats::as.formula(fml_smodel),
                                  data = data)

  # Compute the number of Y-models and S-models
  K_ymodel <- 2 + side * (2 - ER)
  K_smodel <- side

  df <- list(
    N = nrow(data),
    P_ymodel = ncol(X_ymodel),
    P_smodel = ncol(X_smodel),
    Z = Z,
    D = D,
    Y = Y,
    X_ymodel = X_ymodel,
    X_smodel = X_smodel,
    group_idx = group_idx,
    J = J,
    group_levels = group_levels,
    use_hierarchical = use_hierarchical,
    K_ymodel = K_ymodel,
    K_smodel = K_smodel
  )

  return(df)
}

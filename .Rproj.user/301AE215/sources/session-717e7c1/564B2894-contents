#' @title Main function fitting RECKMON (REluCtant Kernel-based Modeling On Non-linearity)
#' @description
#' Implements a three-stage modeling process to capture linear, interaction, and non-linear (Gaussian kernel) effects in the data.
#'
#' @details
#' The function follows a structured, sequential approach:
#'
#' 1.  **Main Effects**: A LASSO regression is first fit on the
#'     response variable `y` against the predictors `x` to model the primary
#'     linear effects. The residuals from this stage are passed to the next.
#'
#' 2.  **Interaction Effects**: The function screens for important interaction
#'     terms by fitting a polynomial kernel and then uses LASSO again to model
#'     the residuals from the first stage. This captures key second-order effects.
#'
#' 3.  **Non-linear Effects**: The remaining residuals are modeled using a
#'     Gaussian kernel regression. The Gaussian kernel's
#'     hyperparameter is tuned automatically using
#'     k-fold cross-validation. This internal CV process is parallelized for
#'     efficiency, assigning each fold-parameter combination to a separate core.
#'     The best parameter is then used to fit the final stepwise Gaussian model
#'     based on an information criterion (BIC).
#'
#' @param x A numeric matrix of predictors. Each row is an observation and
#'   each column is a feature.
#' @param y A numeric vector of the response variable.
#' @param poly_features The number of top-ranked interaction features to select
#'   in the second step. Defaults to $n/log(p)$
#' @param step_gaussian The number of Gaussian kernel features to add at each
#'   step during the final stepwise selection. Defaults to `1`.
#' @param guassian_para_list A numeric vector of candidate values for the
#'   Gaussian kernel hyperparameter (`gamma`). If not provided, a default value
#'   is calculated based on the median distance between data points.
#' @param cv The number of folds for the internal cross-validation used to
#'   tune the Gaussian hyperparameter. Defaults to `5`.
#' @param non_interactions An optional integer vector of feature indices to
#'   exclude from the interaction screening in Stage 2.
#' @param interactions An optional integer vector of feature indices to force
#'   into the interaction screening in Stage 2.
#' @param lambda A small numeric value for the regularization parameter used
#'   when inverting the kernel matrix. Defaults to `1e-3`.
#' @param mc.cores The number of CPU cores to use for parallelizing the internal
#'   cross-validation. Defaults to `1`.
#'
#' @return A list object containing the results of the
#'   modeling process:
#'   \item{BIC}{The final BIC value of the selected Gaussian model.}
#'   \item{poly_feature}{A vector of indices for the selected interaction features.}
#'   \item{gaussian_feature}{A vector of indices for the selected Gaussian kernel features.}
#'   \item{lasso_fit1}{The fitted `glmnet` object from the first stage.}
#'   \item{lasso_fit2}{The fitted `glmnet` object from the second stage.}
#'   \item{theta_hat}{The estimated variance components (`sigma_e^2`, `sigma_g^2`) from the final Gaussian model.}
#'   \item{guassian_para}{The best value selected for the Gaussian kernel hyperparameter.}
#'   \item{r2}{The final residuals after the first two stages, used as input for the third stage.}
#'
#' @export
main_RECKMON <- function(x, y, poly_features = min(ncol(x), floor(nrow(x)/log(ncol(x)))), step_gaussian = 1,
                               guassian_para_list = NULL, cv = 5, non_interactions = NULL, interactions = NULL, lambda = 1e-3, mc.cores = 1){

  # Ensure required packages are available and functions are called with ::
  # In a real package, dependencies are listed in the DESCRIPTION file.

  n <- length(y)

  # Initialize guassian_para_list
  if(is.null(guassian_para_list)){
    guassian_para_min = 1/(median(dist(x)))^2
    guassian_para_list = exp(seq(log(guassian_para_min), log(guassian_para_min*100), length.out = 10))
  }

  # =========================================================================
  # Stage 1: Initial LASSO on y to get first-level residuals (r1)
  # =========================================================================

  lasso_fit1_cv <- glmnet::cv.glmnet(x = x, y = y, family = "gaussian", intercept = TRUE)
  lasso_fit1 <- glmnet::glmnet(x = x, y = y, family = "gaussian", intercept = TRUE, lambda = lasso_fit1_cv$lambda.min)
  r1 <- y - predict(lasso_fit1, x, s = "lambda.min")
  attr(r1, "scaled:center") <- attr(r1, "scaled:scale") <- NULL

  # =========================================================================
  # Stage 2: Screen interaction terms & fit LASSO on r1 to get r2
  # =========================================================================
  K_mat_poly <- kernelMatrix(x, x, kernel = "poly", kparam = 2)
  alpha_hat_poly <- solve(K_mat_poly + nrow(K_mat_poly) * lambda * diag(rep(1, nrow(K_mat_poly)))) %*% r1
  derivative_poly <- pderiv(alpha_hat_poly, as.matrix(x), as.matrix(r1), kernel = "poly", kparam = 2)
  order_poly <- order(derivative_poly, decreasing = TRUE)

  interaction_features <- order_poly[1:poly_features]
  interactions_matrix <- model.matrix(~ .^2, data = as.data.frame(x[, interaction_features]))
  interactions_matrix <- interactions_matrix[, -1] # Remove intercept
  interactions_matrix[, 1:length(interaction_features)] <- x[, interaction_features]^2

  lasso_fit2_cv <- glmnet::cv.glmnet(x = cbind(x, interactions_matrix), y = r1, family = "gaussian", intercept = TRUE)
  lasso_fit2 <- glmnet::glmnet(x = cbind(x, interactions_matrix), y = r1, family = "gaussian", intercept = TRUE, lambda = lasso_fit2_cv$lambda.min)
  r2 <- r1 - predict(lasso_fit2_cv, newx = cbind(x, interactions_matrix), s = "lambda.min")

  # =========================================================================
  # Stage 3a: Parallel Internal CV for Gaussian Hyperparameter
  # =========================================================================
  if(length(guassian_para_list) > 1){
    cat("Starting internal CV for Gaussian parameter tuning...\n")

    chunks <- split(1:n, 1:cv)
    param_grid <- expand.grid(i = 1:length(guassian_para_list), j = 1:cv)

    # parallel::mclapply iterates over each task (a unique fold-parameter combo)
    mse_results <- parallel::mclapply(1:nrow(param_grid), function(task_idx) {
      i <- param_grid$i[task_idx]
      j <- param_grid$j[task_idx]
      current_para <- guassian_para_list[i]

      # Define train/test sets for the current fold
      train_idx <- setdiff(1:n, chunks[[j]])
      test_idx <- chunks[[j]]

      x_train <- x[train_idx, ]; x_test <- x[test_idx, ]
      r2_train <- r2[train_idx]; r2_test <- r2[test_idx]

      # --- Fit the stepwise model on the training fold ---
      # 1. Rank features on training data
      K_mat_g_train <- kernelMatrix(x_train, x_train, kernel = "gaussian", kparam = current_para)
      alpha_hat_g_train <- solve(K_mat_g_train + nrow(K_mat_g_train) * lambda * diag(rep(1, nrow(K_mat_g_train)))) %*% r2_train
      deriv_g_train <- pderiv(alpha_hat_g_train, as.matrix(x_train), as.matrix(r2_train), kernel = "gaussian", kparam = current_para)
      order_g_train <- order(deriv_g_train, decreasing = TRUE)

      # 2. Perform stepwise selection on the training fold
      IC_last_cv <- Inf
      output_p_gaussian_cv <- 0
      model_gaussian_output_cv <- NULL
      p_gaussian_cv <- 0

      # Limit max steps for efficiency during CV
      for(step_cv in 1:50){
        data_temp_cv <- data.frame(r2 = r2_train, id = 1:length(r2_train))

        if (p_gaussian_cv == 0) {
          model_glmmkin_cv <- GMMAT::glmmkin(r2 ~ 1, data = data_temp_cv, id = "id", family = gaussian(link = "identity"))
          K_list_cv <- list(diag(1, nrow = length(r2_train)))
        } else {
          features_to_use <- order_g_train[1:p_gaussian_cv]
          K_gaussian_temp_cv <- kernelMatrix(x_train[, features_to_use, drop = FALSE], x_train[, features_to_use, drop = FALSE], kernel = "gaussian", kparam = current_para)
          colnames(K_gaussian_temp_cv) = rownames(K_gaussian_temp_cv) = 1:length(r2_train)
          model_glmmkin_cv <- GMMAT::glmmkin(r2 ~ 1, data = data_temp_cv, kins = K_gaussian_temp_cv, id = "id", family = gaussian(link = "identity"), method.optim = "Brent")
          K_list_cv <- list(diag(1, nrow = length(r2_train)), K_gaussian_temp_cv)
        }

        Sigma_hat_cv <- Reduce("+", lapply(1:length(K_list_cv), function(k){K_list_cv[[k]] * model_glmmkin_cv$theta[k]}))
        log_likelihood_cv <- -0.5 * t(r2_train) %*% solve(Sigma_hat_cv) %*% r2_train - 0.5 * determinant(Sigma_hat_cv, log = TRUE)$modulus

        # Use BIC for CV model selection regardless of final IC choice
        IC_now_cv <- p_gaussian_cv * log(length(r2_train)) - 2 * log_likelihood_cv

        if (IC_now_cv > IC_last_cv) break

        model_gaussian_output_cv <- model_glmmkin_cv
        output_p_gaussian_cv <- p_gaussian_cv
        IC_last_cv <- IC_now_cv
        p_gaussian_cv <- step_gaussian + p_gaussian_cv
      }

      # --- Predict on the test fold ---
      if (length(model_gaussian_output_cv$theta) > 1 && output_p_gaussian_cv > 0) {
        feature_gaussian <- order_g_train[1:output_p_gaussian_cv]
        K_train_g <- kernelMatrix(x_train[, feature_gaussian, drop = FALSE], x_train[, feature_gaussian, drop = FALSE], kernel = "gaussian", kparam = current_para)
        sigma_Ya_train <- model_gaussian_output_cv$theta[1] * diag(1, nrow(x_train)) + model_gaussian_output_cv$theta[2] * K_train_g

        K_test_train_g <- kernelMatrix(x_test[, feature_gaussian, drop = FALSE], x_train[, feature_gaussian, drop = FALSE], kernel = "gaussian", kparam = current_para)
        sigma_Y_test_train <- model_gaussian_output_cv$theta[2] * K_test_train_g

        gaussian_effect_hat <- sigma_Y_test_train %*% solve(sigma_Ya_train) %*% r2_train
      } else {
        gaussian_effect_hat <- 0
      }

      return(mean((r2_test - gaussian_effect_hat)^2))
    }, mc.cores = mc.cores)

    # Aggregate results and find the best parameter
    mse_matrix <- matrix(unlist(mse_results), nrow = length(guassian_para_list), ncol = cv)
    avg_mse <- rowMeans(mse_matrix, na.rm = TRUE)
    best_para_idx <- which.min(avg_mse)
    guassian_para_best <- guassian_para_list[best_para_idx]

    cat(paste0("...CV complete. Best Gaussian parameter: ", round(guassian_para_best, 4), "\n"))
  }else{
    guassian_para_best = guassian_para_list
  }

  # =========================================================================
  # Stage 3b: FINAL Gaussian stepwise selection using best parameter
  # =========================================================================
  IC_last <- Inf
  output_p_gaussian <- p_gaussian <- 0

  # Rank features on the full dataset using the best parameter
  K_mat_gaussian <- kernelMatrix(x, x, kernel = "gaussian", kparam = guassian_para_best)
  alpha_hat_gaussian <- solve(K_mat_gaussian + n * lambda * diag(rep(1, n))) %*% r2
  derivative_gaussian <- pderiv(alpha_hat_gaussian, as.matrix(x), as.matrix(r2), kernel = "gaussian", kparam = guassian_para_best)
  order_gaussian <- order(derivative_gaussian, decreasing = TRUE)

  model_gaussian_output <- NULL

  for(step3 in 1:200){
    data_temp <- data.frame(r2 = r2, id = 1:n)

    if(p_gaussian == 0){
      model_glmmkin_step3 <- GMMAT::glmmkin(r2 ~ 1, data = data_temp, id = "id", family = gaussian(link = "identity"))
      K_list <- list(diag(1, nrow = n))
    } else {
      features_to_use <- order_gaussian[1:p_gaussian]
      K_gaussian_temp <- kernelMatrix(x[, features_to_use, drop = FALSE], x[, features_to_use, drop = FALSE], kernel = "gaussian", kparam = guassian_para_best)
      colnames(K_gaussian_temp) = rownames(K_gaussian_temp) = 1:n
      model_glmmkin_step3 <- GMMAT::glmmkin(r2 ~ 1, data = data_temp, kins = K_gaussian_temp, id = "id", family = gaussian(link = "identity"), method.optim = "Brent")
      K_list <- list(diag(1, nrow = n), K_gaussian_temp)
    }

    Sigma_hat <- Reduce("+", lapply(1:length(K_list), function(i){K_list[[i]] * model_glmmkin_step3$theta[i]}))
    log_likelihood <- -0.5 * t(r2) %*% solve(Sigma_hat) %*% r2 - 0.5 * determinant(Sigma_hat, logarithm = TRUE)$modulus

    IC_now <- p_gaussian * log(n) - 2 * log_likelihood
    if(IC_now > IC_last) break

    model_gaussian_output <- model_glmmkin_step3
    output_p_gaussian <- p_gaussian
    IC_last <- IC_now
    p_gaussian <- step_gaussian + p_gaussian
  }

  cat(paste0("Final model fitted. p_poly: ", length(interaction_features), "; p_gaussian: ", output_p_gaussian, "\n"))

  final_gaussian_features <- if(output_p_gaussian > 0) order_gaussian[1:output_p_gaussian] else NULL

  # Compile and return results
  results <- list(
    BIC = IC_last,
    poly_feature = interaction_features,
    gaussian_feature = final_gaussian_features,
    lasso_fit1 = lasso_fit1,
    lasso_fit2 = lasso_fit2,
    theta_hat = model_gaussian_output$theta,
    guassian_para = guassian_para_best,
    r2 = r2
  )

  # class(results) <- "msm_fit" # Assign a class for custom print methods, etc.
  return(results)
}

#' Predict from a Fitted RECKMON Model
#'
#' @description
#' This function generates predictions for new data using a model object
#' previously fitted by `main_RECKMON`.
#'
#' @param object A fitted model as returned by
#'   `main_RECKMON()`.
#' @param newx A numeric matrix containing the new data for which to generate
#'   predictions. Each row is an observation and columns must match the
#'   original training data.
#'
#' @return A numeric vector of predicted values.
#'
#' @export

predict_RECKMON <- function(newx, model_fit){
  if(is.vector(newx)) newx = matrix(newx, nrow = 1)
  x = model_fit$x
  y = model_fit$y
  n = nrow(x)
  newn = nrow(newx)
  x_merge = rbind(x, newx)
  train_seq = 1:n

  if(length(model_fit$theta_hat) == 1){
    gaussian_effect_hat = 0
  }else{
    feature_gaussian = model_fit$gaussian_feature
    OmicsKernelMatrix <- vector("list", length = 2)
    OmicsKernelMatrix[[1]] <- diag(1, nrow = nrow(x_merge))
    OmicsKernelMatrix[[2]] <- kernelMatrix(x_merge[,feature_gaussian], x_merge[,feature_gaussian], kernel = "gaussian", kparam = model_fit$guassian_para)
    sigma_Ya = Reduce("+", lapply(1:length(OmicsKernelMatrix), function(i){OmicsKernelMatrix[[i]] * model_fit$theta_hat[i]}))
    gaussian_effect_hat = sigma_Ya[-train_seq,train_seq] %*% solve(sigma_Ya[train_seq,train_seq]) %*% model_fit$r2
  }
  interactions_matrix_test <- model.matrix(~ .^2, data = as.data.frame(newx[, model_fit$poly_feature, drop = F ]))
  interactions_matrix_test <- interactions_matrix_test[,-1, drop = F]
  interactions_matrix_test[,1:length(model_fit$poly_feature) ] = newx[, model_fit$poly_feature, drop = F ]^2
  y_predict <- predict(model_fit$lasso_fit1, newx) +
    predict(model_fit$lasso_fit2, newx = cbind(newx, interactions_matrix_test)) +
    gaussian_effect_hat

  return(y_predict)
}





###
cv_asgl <- function(x, y, index, family = c("gaussian", "binomial"), offset = NULL, 
                    alpha = 0.95, lambda = NULL, lambda_min = 0.1, nlambda = 20, 
                    maxit = 1000, thresh = 0.001, gamma = 0.8, step = 1, standardize = FALSE, 
                    grp_weights = NULL, ind_weights = NULL, nfolds = 5) {
  if (!is.matrix(x)) {
    stop("the argument 'x' must be a matrix.", call. = FALSE)
  }
  dimx <- dim(x)
  if (is.null(dimx) || dimx[2] < 2) {
    stop("the argument 'x' must be a matrix with 2 or more columns.", call. = FALSE)
  }
  nobs <- dimx[1]
  nvar <- dimx[2]
  if (!is.vector(y)) {
    stop("the argument 'y' must be a vector.", call. = FALSE)
  }
  leny <- length(y)
  if (leny != nobs) {
    stop(paste("the length of 'y' (", leny, ") is not equal to the number of ", 
               "rows of 'x' (", nobs, ").", sep = ""), call. = FALSE)
  }
  if (!is.vector(index)) {
    stop("the argument 'index' must be a vector.", call. = FALSE)
  }
  leni <- length(index)
  if (leni != nvar) {
    stop(paste("the length of 'index' (", leni, ") is not equal to the number ", 
               "of columns of 'x' (", nvar, ").", sep = ""), call. = FALSE)
  }
  if (is.null(offset)) {
    offset <- rep.int(0, leny)
  }
  if (!is.vector(offset)) {
    stop("the argument 'offset' must be a vector.", call. = FALSE)
  }
  leno <- length(offset)
  if (leno != nobs) {
    stop(paste("the length of 'offset' (", leno, ") is not equal to the ", 
               "number of rows of 'x' (", nobs, ").", sep = ""), call. = FALSE)
  }
  lenui <- length(unique(index))
  if (is.null(grp_weights)) {
    grp_weights <- rep.int(1, lenui)
  }
  if (!is.vector(grp_weights)) {
    stop("the argument 'grp_weights' must be a vector.", call. = FALSE)
  }
  if (is.null(ind_weights)) {
    ind_weights <- rep.int(1, nvar)
  }
  if (!is.vector(ind_weights)) {
    stop("the argument 'ind_weights' must be a vector.", call. = FALSE)
  }
  lengw <- length(grp_weights)
  if (lengw != lenui) {
    stop(paste("the length of 'grp_weights' (", lengw, ") is not equal to the ", 
               "number of unique elements of 'index' (", lenui, 
               ").", sep = ""), call. = FALSE)
  }
  leniw <- length(ind_weights)
  if (leniw != nvar) {
    stop(paste("the length of 'ind_weights' (", leniw, ") is not equal to the ", 
               "number of columns of 'x' (", nvar, ").", sep = ""), call. = FALSE)
  }
  family <- match.arg(family)
  if (alpha < 0 || alpha > 1) {
    stop("the argument 'alpha' must be between 0 and 1.")
  }
  if (is.null(lambda)) {
    lambda <- get_lambda_sequence(x, y, index, family, lambda_min, 
                                  nlambda, alpha, grp_weights, ind_weights)
  } else {
    nlambda <- length(lambda)
  }
  foldid <- sample(rep(1:nfolds, length.out = nobs))
  cv_errors <- matrix(NA, nrow = nlambda, ncol = nfolds)
  for (i in 1:nfolds) {
    train <- (foldid != i)
    test <- (foldid == i)
    fit <- asgl(x[train, ], y[train], index, family = family, offset = offset[train], 
                alpha = alpha, lambda = lambda, lambda_min = lambda_min, nlambda = nlambda, 
                maxit = maxit, thresh = thresh, gamma = gamma, step = step, 
                standardize = standardize, grp_weights = grp_weights, ind_weights = ind_weights)   
    for (j in 1:nlambda) {
      beta <- fit$beta[, j]
      intercept <- fit$intercept[j]
      if (standardize && !is.null(fit$X.transform)) {
        x_test <- scale(x[test, ], center = fit$X.transform$X.means, scale = fit$X.transform$X.scale)
      } else {
        x_test <- x[test, ]
      }
      y_pred <- x_test %*% beta + intercept + offset[test]
      if (family == "gaussian") {
        cv_errors[j, i] <- mean((y[test] - y_pred)^2)
      } else if (family == "binomial") {
        prob <- 1 / (1 + exp(-y_pred))
        eps <- 1e-5  
        prob <- pmin(pmax(prob, eps), 1 - eps)
        cv_errors[j, i] <- -mean(y[test] * log(prob) + (1 - y[test]) * log(1 - prob))
      }
    }
  }
  mean_cv_errors <- rowMeans(cv_errors)
  best_lambda_index <- which.min(mean_cv_errors)
  best_lambda <- lambda[best_lambda_index]
  best_fit <- asgl(x, y, index, family = family, offset = offset, 
                   alpha = alpha, lambda = lambda, lambda_min = lambda_min, 
                   nlambda = nlambda, maxit = maxit, thresh = thresh, 
                   gamma = gamma, step = step, standardize = standardize, 
                   grp_weights = grp_weights, ind_weights = ind_weights)
  return(list(fit = best_fit, lambda = best_lambda,
              lambda_index = best_lambda_index, cv_errors = cv_errors,
              mean_cv_errors = mean_cv_errors))
}
plot_cv_errors <- function(cv_result) {
  lambda <- cv_result$fit$lambda
  cv_errors <- cv_result$cv_errors
  mean_cv_errors <- rowMeans(cv_errors)
  plot(log(lambda), mean_cv_errors, type = "b", xlab = "log(lambda)", ylab = "Mean CV Error",
       main = "Cross-Validation Error vs. log(lambda)")
  abline(v = log(cv_result$lambda), col = "red", lty = 2)
  legend("topright", legend = paste("Best lambda =", round(cv_result$lambda, 4)), 
         col = "red", lty = 2)
}
plot_coefficients <- function(cv_result) {
  best_lambda_index <- which.min(rowMeans(cv_result$cv_errors))
  beta <- cv_result$fit$beta[, best_lambda_index]
  feature_names <- colnames(cv_result$fit$x)
  if (is.null(feature_names)) {
    feature_names <- paste("X", 1:length(beta), sep = "")
  }  
  plot(beta, type = "h", lwd = 2, xlab = "Features", ylab = "Coefficients",
       main = paste("Coefficients at Best Lambda (", round(cv_result$lambda, 4), ")", sep = ""))
  axis(1, at = 1:length(beta), labels = feature_names, las = 2, cex.axis = 0.7)
}


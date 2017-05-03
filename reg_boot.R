# ################################################################################# #
# SCRIPT_NAME: Bag of little bootstraps (BLB)
#		The BLB is a procedure to produce robust, computationally efficient 
#		bootstrap estimates for a statistic of interest.
# INPUTS
#	- data: x_1, x_2,..., x_n
#	- theta_hat: estimator of interest
# 	- m: subset size
#	- s: number of sampled subsets
#	- B: number of Monte Carlo iterations
#	- q: estimator of quality assessment
# OUTPUTS
#	- estimate of q, i.e., q(Q_n(P))
# ################################################################################# #
## BEGIN SCRIPT ------------------------------------------------------------------- #
# --------------------------------------------------------------------------------- #
# regular bootstrap
reg_boot <- function(X, y, B){
      set.seed(5)
      d <- dim(X)[2]
      beta_ls_boots <- matrix(rep(0, B*d), ncol = d)
      ls_fit <- lm(y ~ X - 1)
      beta_ls <- ls_fit$coefficients
      res_ls <- ls_fit$residuals
      std_err <- summary(ls_fit)$coefficients[,2]
      for(i in 1:B){
            idx <- sample(n, n, replace = TRUE)
            y_boot <- X[idx,] %*% beta_ls + res_ls[idx]
            beta_ls_boots[i,] <- lsfit(X[idx,], y_boot, 
                                       intercept = FALSE)$coefficients
      }
      beta_ls_mean <- colMeans(beta_ls_boots)
      q_boot <- apply(beta_ls_boots, 2, sd)
      q_hi <- apply(beta_ls_boots, 2, quantile, probs = 0.975)
      q_lo <- apply(beta_ls_boots, 2, quantile, probs = 0.025)
      return(list(q_boot = q_boot, q_hi = q_hi, q_lo = q_lo, std_err = std_err,
                  beta_ls_boots = beta_ls_boots, beta_ls = beta_ls))
}
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
log_reg_boot <- function(X, y, B){
      set.seed(5)
      d <- dim(X)[2]
      n <- dim(X)[1]
      beta_glm_boots <- matrix(rep(0, B*d), ncol = d)
      glm_fit <- glm(y ~ X - 1, family = "binomial")
      beta_glm <- glm_fit$coefficients
      pred_glm <- glm_fit$fitted.values
      for(i in 1:B){
            idx <- sample(n, n, replace = TRUE)
            y_boot <- rbinom(n, 1, pred_glm)
            beta_glm_boots[i,] <- glm(y_boot[idx] ~ X[idx,] - 1, 
                                      family = "binomial")$coeff
      }
      beta_glm_mean <- colMeans(beta_glm_boots)
      q_boot <- apply(beta_glm_boots, 2, sd)
      q_hi <- apply(beta_glm_boots, 2, quantile, probs = 0.975)
      q_lo <- apply(beta_glm_boots, 2, quantile, probs = 0.025)
      return(list(q_boot = q_boot, q_hi = q_hi, q_lo = q_lo, 
                  beta_glm_mean = beta_glm_mean, beta_glm = beta_glm,
                  beta_glm_boots = beta_glm_boots, pred_glm = pred_glm))
}
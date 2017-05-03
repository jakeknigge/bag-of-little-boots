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
# BLB function
log_reg_blb <- function(X, y, B, ss, m, s){
      set.seed(3)
      n <- dim(X)[1]
      d <- dim(X)[2]
      # multinomial weighting vector
      w <- rep(1/m, m)
      beta_boot <- matrix(rep(0, B*d), ncol = d)
      beta_blb_boots <- matrix(rep(0, s*d), ncol = d)
      preds <- rep(0, m)
      q <- matrix(rep(0, s*d), ncol = d)
      q_h <- matrix(rep(0, s*d), ncol = d)
      q_l <- matrix(rep(0, s*d), ncol = d)
      
# --------------------------------------------------------------------------------- #
      # BLB
      for(j in 1:s){
            # subsample the data
            idx <- sample(n, m, replace = FALSE)
            # approximate quality assessment of estimator
            for(k in 1:B){
                  if(k == 1){
                        glm_fit <- glm(y[idx] ~ X[idx,] - 1, family = "binomial")
                        beta_glm <- glm_fit$coefficients
                        pred_glm <- glm_fit$fitted.values
                  }
                  y_boot <- rbinom(m, 1, pred_glm)
                  counts <- as.vector(rmultinom(1, n, w))
                  beta_boot[k,] <-glm(y_boot ~ X[idx,] - 1, 
                                    family = "binomial", 
                                    weights = counts)$coeff
            }
            #preds <- (1/j) * pred_glm + ((j-1)/j) * preds
            beta_blb_boots[j, ] <- colMeans(beta_boot)
            q[j, ] <- apply(beta_boot, 2, sd)
            q_l[j, ] <- apply(beta_boot, 2, quantile, probs = 0.025)
            q_h[j, ] <- apply(beta_boot, 2, quantile, probs = 0.975)
      }
      q_means <- colMeans(q)
      q_hi <- colMeans(q_h)
      q_lo <- colMeans(q_l)
      beta_blb_mean <- colMeans(beta_blb_boots)
      return(list(q_means = q_means, q_lo = q_lo, q_hi = q_hi,
                  beta_blb_boots = beta_blb_boots, beta_glm = beta_glm,
                  beta_blb_mean = beta_blb_mean,
                  preds = preds))
}
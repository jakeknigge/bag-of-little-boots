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
reg_blb <- function(X, y, B, ss, m, s){
      set.seed(3)
      d <- dim(X)[2]
      # multinomial weighting vector
      w <- rep(1/m, m)
      beta_boot <- matrix(rep(0, B*d), ncol = d)
      beta_blb_boots <- matrix(rep(0, s*d), ncol = d)
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
                        ls_fit <- lsfit(X[idx,], y[idx], intercept = FALSE)
                        beta_ls <- ls_fit$coefficients
                        res_ls <- ls_fit$residuals
                  } 
                  counts <- as.vector(rmultinom(1, n, w))
                  wX <- sqrt(counts) * X[idx,]
                  wy <- wX %*% beta_ls + res_ls * sqrt(counts)
                  beta_boot[k,] <- lsfit(wX, wy, intercept = FALSE)$coefficients
            }
            beta_blb_boots[j, ] <- colMeans(beta_boot)
            q[j, ] <- apply(beta_boot, 2, sd)
            q_l[j, ] <- apply(beta_boot, 2, quantile, probs = 0.025)
            q_h[j, ] <- apply(beta_boot, 2, quantile, probs = 0.975)
      }
      q_means <- colMeans(q)
      q_hi <- colMeans(q_h)
      q_lo <- colMeans(q_l)
      return(list(q_means = q_means, q_lo = q_lo, q_hi = q_hi,
                  beta_blb_boots = beta_blb_boots))
}
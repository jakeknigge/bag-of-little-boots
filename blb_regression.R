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
# generate data
rm(list = ls())
set.seed(2)
n <- 20000
d <- 200
ndr <- n / d
X <- matrix(rt(n*d, 3), nrow = n)
#X <- X %*% diag(1/colSums(X)) # normalize columns of data matrix
X <- t(t(X) * (1/colSums(X)))
beta <- rep(1, d) # runif(d, 0, 5) # rnorm(d)
w <- sqrt(0.1^2)*rnorm(n, 0, 1) # 0.05 for 20000 by 100
y <- X %*% beta + w
snr <- sum((X %*% beta)^2) / sum(w^2)

#ls_fit <- lsfit(X, y, intercept = FALSE)
#beta_ls <- ls_fit$coefficients
#res_ls <- ls_fit$residuals

lm_fit <- lm(y ~ X - 1)
std_err <- summary(lm_fit)$coefficients[,2]
rm(lm_fit)

#hist(y)
#plot(1:d, beta, type = "h")
#points(1:d, beta, pch = 16)
#points(1:d, beta_ls, pch = 16, col = 'red3')
#hist(res_ls)

# --------------------------------------------------------------------------------- #
# BLB parameters
# subset size
ss <- 0.7
m <- ceiling(n^ss)
#m <- max(ceiling(n^ss), ceiling(d /(1 - exp(-1))))
# number of subsets
s <- ceiling(n / m)
# bootstrap iterations
B <- 50
# multinomial weighting vector
w <- rep(1/m, m)
beta_boot <- matrix(rep(0, B*d), ncol = d)
beta_boots <- matrix(rep(0, s*d), ncol = d)
q <- matrix(rep(0, s*d), ncol = d)

# --------------------------------------------------------------------------------- #
# BLB
tic <- proc.time()
set.seed(3)
for(j in 1:s){
      # subsample the data
      idx <- sample(n, m, replace = FALSE)
      # approximate quality assessment of estimator
      for(k in 1:B){
            if(k == 1){
                  ls_fit <- lsfit(X[idx,], y[idx], intercept = FALSE)
                  beta_ls <- ls_fit$coefficients
                  res_ls <- ls_fit$residuals
                  #counts <- as.vector(rmultinom(1, n, w))
                  #wX <- sqrt(counts) * X[idx,]
                  #wy <- wX %*% beta_ls + res_ls * sqrt(counts)
                  #beta_boot[k,] <- lsfit(wX, wy, intercept = FALSE)$coefficients
            } 
            counts <- as.vector(rmultinom(1, n, w))
            wX <- sqrt(counts) * X[idx,]
            wy <- wX %*% beta_ls + res_ls * sqrt(counts)
            beta_boot[k,] <- lsfit(wX, wy, intercept = FALSE)$coefficients
      }
      beta_boots[j, ] <- colMeans(beta_boot)
      q[j, ] <- apply(beta_boot, 2, sd)
}
q_means <- colMeans(q)
toc_blb <- (proc.time() - tic)[3]

res_reg_blb <- reg_blb(X, y, B, ss, m, s)

# --------------------------------------------------------------------------------- #
# regular bootstrap
beta_ls_boots <- matrix(rep(0, B*d), ncol = d)
tic <- proc.time()
ls_fit <- lsfit(X, y, intercept = FALSE)
beta_ls <- ls_fit$coefficients
res_ls <- ls_fit$residuals
for(i in 1:B){
      idx <- sample(n, n, replace = TRUE)
      y_boot <- X[idx,] %*% beta_ls + res_ls[idx]
      beta_ls_boots[i,] <- lsfit(X[idx,], y_boot, intercept = FALSE)$coefficients
}
beta_ls_mean <- colMeans(beta_ls_boots)
q_boot <- apply(beta_ls_boots, 2, sd)
toc_boot <- (proc.time() - tic)[3]

res_reg_boot <- reg_boot(X, y, B)

# --------------------------------------------------------------------------------- #
# parameter vector comparisons
plot(1:d, beta, pch = 1, cex = 1.5)
points(1:d, colMeans(beta_boots) + q_means, pch = 18, col = 'orangered1')
points(1:d, colMeans(beta_boots) - q_means, pch = 18, col = 'orangered1')
segments(1:d, colMeans(beta_boots) - q_means, 1:d, colMeans(beta_boots) + q_means,
            col = 'orangered1')
#segments(1:d, beta_ls - std_err, 1:d, beta_ls + std_err,
#         col = 'blue4', lty = 2)
#points(1:d, beta_ls + std_err, pch = 23, col = 'blue4')
#points(1:d, beta_ls - std_err, pch = 23, col = 'blue4')
points(1:d, colMeans(beta_boots), pch = 19, col = 'darkorange')
points(1:d, beta_ls, pch = 23, col = 'royalblue1', lwd = 1.15, cex = 1.15)

# --------------------------------------------------------------------------------- #
# parameter vector norm differences
sqrt(sum(beta - beta_ls)^2) # full least squares fit
sqrt(sum(beta - colMeans(beta_boots))^2) # blb fit
sqrt(sum(beta - beta_ls_boots)^2) # bootstrap averaged fit

# --------------------------------------------------------------------------------- #
# standard error comparisons
plot(1:d, std_err, pch = 1, lwd = 1.15, cex = 1.15)
points(1:d, q_means, pch = 19, col = 'darkorange')
points(1:d, q_boot, pch = 23, col = 'royalblue3')
# ################################################################################# #
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
# clear out the garbage
rm(list = ls())
source("reg_boot.R")
source("reg_blb.R")

# generate data
set.seed(2)
n <- 20000
d <- 100
ndr <- n / d
X <- matrix(rt(n*d, 3), nrow = n)
X <- t(t(X) * (1/colSums(X)))
beta <- rnorm(d) # rep(1, d) # runif(d, 0, 5)
w <- sqrt(0.01)*rnorm(n, 0, 1)
y <- X %*% beta + w
snr <- sum((X %*% beta)^2) / sum(w^2)

# --------------------------------------------------------------------------------- #
# BLB parameters
# subset size
ss <- 0.7
m <- ceiling(n^ss)
#m <- max(ceiling(n^ss), ceiling(d /(1 - exp(-1))))
# number of subsets
s <- min(ceiling(n / m), 15)
# bootstrap iterations
B <- 100

# --------------------------------------------------------------------------------- #
tic <- proc.time()
res_reg_blb <- reg_blb(X, y, B, ss, m, s)
toc_blb <- (proc.time() - tic)[3]
q_means <- res_reg_blb$q_means
q_hi_blb <- res_reg_blb$q_hi
q_lo_blb <- res_reg_blb$q_lo
beta_blb_boots <- res_reg_blb$beta_blb_boots
beta_blb_means <- colMeans(beta_blb_boots)

tic <- proc.time()
res_reg_boot <- reg_boot(X, y, B)
toc_boot <- (proc.time() - tic)[3]
q_boot <- res_reg_boot$q_boot
q_hi_boot <- res_reg_boot$q_hi
q_lo_boot <- res_reg_boot$q_lo
std_err <- res_reg_boot$std_err
beta_ls_boots <- res_reg_boot$beta_ls_boots
beta_ls <- res_reg_boot$beta_ls

tic <- proc.time()
ls_fit <- lm(y ~ X - 1)
toc_ls_fit <- (proc.time() - tic)[3]

# --------------------------------------------------------------------------------- #
par(mfrow = c(1,1))
# parameter vector comparisons
plot(1:d, beta, pch = 0, cex = 1)
points(1:d, beta_blb_means + q_means, pch = 18, col = 'orangered1')
points(1:d, beta_blb_means - q_means, pch = 18, col = 'orangered1')
segments(1:d, beta_blb_means - q_means, 1:d, 
         beta_blb_means + q_means, col = 'orangered1')
#segments(1:d, beta_ls - std_err, 1:d, beta_ls + std_err,
#         col = 'blue4', lty = 2)
#points(1:d, beta_ls + std_err, pch = 23, col = 'blue4')
#points(1:d, beta_ls - std_err, pch = 23, col = 'blue4')
points(1:d, beta_blb_means, pch = 19, col = 'darkorange')
points(1:d, beta_ls, pch = 23, col = 'royalblue1', lwd = 1.15, cex = 1.15)

# --------------------------------------------------------------------------------- #
# parameter vector norm differences
sqrt(sum(beta - beta_ls)^2) # full least squares fit
sqrt(sum(beta - beta_blb_means)^2) # blb fit
sqrt(sum(beta - beta_ls_boots)^2) # bootstrap averaged fit

#plot(1:d, q_hi_blb, pch = 15, col = 'darkorange3',
#     ylim=c(min(q_lo_blb, q_lo_boot), max(q_hi_blb, q_hi_boot)))
plot(1:d, beta, ylim=c(min(q_lo_blb), max(q_hi_blb)), pch = 0)
points(1:d, q_hi_blb, pch = 15, col = 'orangered3')
points(1:d, q_lo_blb, pch = 15, col = 'orangered3')
segments(1:d, q_lo_blb, 1:d, q_hi_blb, col = 'orangered3')
points(1:d, beta_blb_means, pch = 16, col = 'darkorange')
points(1:d, beta_blb_means + q_means, pch = 18, col = 'orangered1')
points(1:d, beta_blb_means - q_means, pch = 18, col = 'orangered1')
points(1:d, q_hi_boot, pch = 15, col = 'royalblue3')
points(1:d, q_lo_boot, pch = 15, col = 'royalblue3')
points(1:d, beta_ls, pch = 16, col = 'royalblue1')
points(1:d, beta_ls + std_err, pch = 18, col = 'royalblue1')
points(1:d, beta_ls - std_err, pch = 18, col = 'royalblue1')

plot(1:d, std_err, pch = 0)
points(1:d, q_means, pch = 16, col = 'orangered')
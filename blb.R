# ################################################################################# #
# SCRIPT_NAME: Bag of little bootstraps (BLB)
#		The BLB is a procedure to produce robust, computationally efficient 
#		bootstrap estimates for a statistic of interest.
# INPUTS
#	- data: x_1, x_2,..., x_n
#	- theta_hat: estimator of interest
# 	- b: subset size
#	- s: number of sampled subsets
#	- r: number of Monte Carlo iterations
#	- q: estimator of quality assessment
# OUTPUTS
#	- estimate of q, i.e., q(Q_n(P))
# ################################################################################# #
## BEGIN SCRIPT ------------------------------------------------------------------- #
library(boot)
## BAG OF LITTLE BOOTSTRAPS FUNCTION
# --------------------------------------------------------------------------------- #
blb <- function(data, theta_hat, r, b, s, q){
# data dimensions
q <- rep(0, s)
Q <- rep(0, r)
theta_boots <- rep(0, s)
n <- NROW(data)
d <- NCOL(data)
w <- rep(1/b, b)
	# loop for little boots
	for(i in 1:s){
		# subsample the data
		idx <- sample(n, b)
		# tboot <- boot(data = data[idx], statistic = theta_hat, R = r)
		# Q <- tboot$t / r
		# theta_boots[i] <- tboot$t0
		for(j in 1:r){
			w_counts <- as.vector(rmultinom(1, n, w))
			Q[j] <- median(w_counts/n) * sum(data[idx])
		} # j for loop
		q[i] <- sd(Q)
	} # i for loop
q_boots <- q
q_boot <- mean(q)
return(list(q_boot = q_boot, q_boots = q_boots, theta_boots = theta_boots))
} # function

## END OF SCRIPT ------------------------------------------------------------------ #
# ################################################################################# #

# MEDIAN EXAMPLE
n_samples <- 100000
x <- rnorm(n_samples, 5, 2)
b <- ceiling(n_samples^0.7)
s <- ceiling(n_samples / b)
r <- 100
median_boot <- function(x, i){median(x[i])}
blb_out <- blb(data = x, theta_hat = median_boot, r = r, 
			b = b, s = s, q = NA)
boot_out <- boot(x, median_boot, R = r)

blb_out$q_boot
boot_out

# LINEAR REGRESSION - BOOTSTRAP WEIGHTING PROOF OF CONCEPT

set.seed(2)
n <- 10
d <- 5
X <- matrix(rnorm(n*d), nrow = n)
beta <- runif(d, 5, 20)
y <- X %*% beta + rnorm(n, 0, 0.1)
fit <- lm(y ~ X - 1)
beta_ls <- fit$coeff
residuals <- fit$resid

idx <- c(1, 1, 2, 3, 4, 4, 6, 8, 8, 9)
u_idx <- unique(idx)
counts <- c(2, 1, 1, 2, 1, 2, 1)
#counts <- c(2, 1, 1, 2, 0, 1, 0, 2, 1, 0)

BX <- X[idx,]
y_boot <- BX %*% beta_ls + residuals[idx]
# y_boot <- y[idx]
beta_boot <- lm(y_boot ~ BX - 1)$coeff

WX <- sqrt(counts) * X[u_idx,]
#WX <- counts * X
#WX <- WX[u_idx,]
y_alt <- WX %*% beta_ls + residuals[u_idx] * sqrt(counts)
#y_alt <- WX %*% beta_ls + residuals[u_idx] * counts
#y_alt <- sqrt(counts) * y[u_idx]
beta_alt <- lm(y_alt ~ WX - 1)$coeff

plot(beta, pch = 0, col = "black", main = "estimate comparisons")
points(beta_ls, pch = 16, col = "royalblue")
points(beta_boot, pch = 1, col = "green4")
points(beta_alt, pch = 4, col = "red")
beta
beta_ls
beta_boot
beta_alt

rm(list = ls())
set.seed(2)
n <- 20000
d <- 20
ndr <- n / d
X <- matrix(rnorm(n*d), nrow = n)
beta <- rnorm(d)
locs <- as.vector(X %*% beta)
probs <- 1/(1 + exp(-locs))
u <- runif(n)
y <- as.numeric(u < probs)

# --------------------------------------------------------------------------------- #
# BLB parameters
# subset size
ss <- 0.7
m <- ceiling(n^ss)
# number of subsets
s <- min(ceiling(n / m), 15)
# bootstrap iterations
B <- 100

source("log_reg_boot.R")
tic <- proc.time()
out_boot <- log_reg_boot(X, y, B)
toc_boot <- (proc.time() - tic)[3]

source("log_reg_blb.R")
tic <- proc.time()
out_blb <- log_reg_blb(X, y, B, ss, m, s)
toc_blb <- (proc.time() - tic)[3]

sqrt(sum(beta - out_blb$beta_blb_mean)^2)
sqrt(sum(beta - out_boot$beta_glm)^2)

par(mfrow = c(1,1))
plot(beta, pch = 0, 
     ylim=c(min(out_blb$beta_blb_mean, out_boot$beta_glm), 
            max(out_blb$beta_blb_mean, out_boot$beta_glm)))
points(out_blb$beta_blb_mean, pch = 16, col = 'darkorange')
points(out_boot$beta_glm, pch = 18, col = 'royalblue')

plot(beta, ylim=c(min(out_blb$q_lo), max(out_blb$q_hi)), type = "n")
points(out_blb$beta_blb_mean, pch = 16, col = 'darkorange')
segments(1:d, out_blb$q_lo, 1:d, out_blb$q_hi, col = 'darkorange3')
points(out_blb$q_lo, pch = 15, col = 'darkorange3')
points(out_blb$q_hi, pch = 15, col = 'darkorange3')
#points(out_blb$beta_blb_mean + out_blb$q_means, pch = 18, col = 'darkorange2')
#points(out_blb$beta_blb_mean - out_blb$q_means, pch = 18, col = 'darkorange2')
points(beta, pch = 0)
set.seed(2)
n <- 2000
d <- 2
X <- matrix(rnorm(n*d), nrow = n)
beta <- rnorm(d)
locs <- as.vector(X %*% beta)
probs <- 1/(1 + exp(-locs))
u <- runif(n)
y <- as.numeric(u < probs)
fit <- glm(y ~ X-1, family = "binomial")
beta_glm <- fit$coeff
pred <- fit$fitted.values

idx <- sort(sample(n, n, replace = TRUE))
u_idx <- unique(idx)
counts <- rle(idx)$lengths

y_boot <- rbinom(n, 1, pred)
beta_boot <- glm(y_boot[idx] ~ X[idx,]-1, family = "binomial")$coeff
beta_alt <- glm(y_boot[u_idx] ~ X[u_idx,] - 1, 
                family = "binomial", weights = counts)$coeff

plot(beta, pch = 0, col = "black", main = "estimate comparisons")
points(beta_glm, pch = 16, col = "royalblue")
points(beta_boot, pch = 1, col = "green4")
points(beta_alt, pch = 4, col = "red")
beta
beta_glm
beta_boot
beta_alt

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

source("log_reg_boot.R")
out_boot <- log_reg_boot(X, y, B)

source("log_reg_blb.R")
out_blb <- log_reg_blb(X, y, B, ss, m, s)

plot(beta, pch = 0)
points(out_blb$beta_blb_mean, pch = 16, col = 'darkorange')
points(out_boot$beta_glm, pch = 18, col = 'royalblue')
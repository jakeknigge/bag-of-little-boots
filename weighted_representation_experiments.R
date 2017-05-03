# Weighted representation expirements

# Median
set.seed(2)
n <- 35
b <- 11
counts <- rmultinom(1, n, rep(1/b, b))
weights <- counts / sum(counts)
x <- rnorm(b)
ord <- order(x)
plot(sort(x), weights[ord], type = "h")
points(sort(x), weights[ord], pch = 20)
median(weights)
median(x)
which(sort(x) == median(x))

# Regression
set.seed(22)
n <- 100
d <- 5
X <- matrix(c(rep(1, n), rt(n*d, 3)), nrow = n)
beta <- rnorm(d+1)
y <- X %*% beta + rnorm(100,0,0.5)
plot(y)
fit <- lm(y ~ X-1)
boot_idx <- sample(n, n, replace = TRUE)


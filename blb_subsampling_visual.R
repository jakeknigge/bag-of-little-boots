# ################################################################################# #
# SCRIPT_NAME:	BLB SUBSAMPLING VISUAL
#	provides visuals that illustrate what the bootstrap and bag of little 
#	bootstrap are doing at each sampling iteration
# ################################################################################# #
## BEGIN SCRIPT ------------------------------------------------------------------- #
# library(weights)
rm(list = ls())
set.seed(2)

# parameters
n <- 1000 		  # data size
m <- ceiling(n^0.7) # subsample size
s <- ceiling(n/m)	  # number of subsamples

# generate data and bin it
x <- c(rnorm(n/2, 2, 1), rnorm(n/2, -2, 1))
h <- hist(x, plot = FALSE)

# sample / subsample indices
idx_blb <- sample(n, m)
idx_boot <- sample(n, n, replace = TRUE)
counts <- rmultinom(1, n, prob = rep(1/m, m))

# sanity checks
length(counts) 	# same as m
sum(counts) 	# same as n

# --------------------------------------------------------------------------------- #
# visualizations
#dev.new(height = 12, width = 12)
par(mfcol=c(2,2))
#par(mfrow=c(4,1))
hist(x, breaks = h$breaks)
rug(x, lwd = 1)
rug(x[idx_blb], col = 'darkorange2', lwd = 1)
hist(x[idx_boot], breaks = h$breaks, 
     main = "Histogram of bootstrap sample")
rug(x, lwd = 1)
rug(x[idx_boot], col = 'green3', lwd = 1)
hist(x[idx_blb], breaks = h$breaks, 
     main = "Histogram of 'little' bootstrap sample")
rug(x[idx_blb], col = 'darkorange2', lwd = 1)
wtd.hist(x[idx_blb], weight = counts, breaks = h$breaks,
         main = "Histogram of weighted 'little' bootstrap sample")
rug(x[idx_blb], col = 'darkorange2', lwd = 1)

# --------------------------------------------------------------------------------- #
# little bootstrap samples
par(mfrow=c(2,1))
p <- 20
plot(1:p, seq(1, n, length.out = p), type = "n", xlab = "sample number",
     ylab = "selected indices", main = "Little bootstraps - indices")
for(i in 1:p){
      points(rep(i, m), sample(n, m), pch = 16, 
             cex = 0.65, col = 'darkorange2')
}
plot(1:p, seq(1, n, length.out = p), type = "n", xlab = "sample number",
     ylab = "x values", main = "Little bootstraps - locations",
     ylim = c(min(x), max(x)))
for(i in 1:p){
      points(rep(i, m), x[sample(n, m)], pch = 16, 
             cex = 0.65, col = 'darkorange2')
}
# --------------------------------------------------------------------- END SCRIPT ##
# ################################################################################# #
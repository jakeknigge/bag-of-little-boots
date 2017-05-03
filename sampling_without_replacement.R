# ################################################################################# #
# SCRIPT_NAME:	SAMPLING WITHOUT REPLACEMENT
#			Performs weighted sampling without replacement using R's built
#			in 'sample' function; returns plots of the conditional
#			distributions for each sampled index.
# ################################################################################# #
## BEGIN SCRIPT ------------------------------------------------------------------- #
# generate data
set.seed(3)
d = 10
x = 1:d
size = d

#probs = sample(1:4, d, replace = TRUE)
#probs = probs / sum(probs)

# use same weights as Matlab script
w = c(2,1,3,2,2,2,1,3,2,2)
w = w / sum(w)

n_draws = 10000
draws = matrix(rep(0, d*n_draws), ncol = d)

# sample without replacement from weighted distribution
for(i in 1:n_draws){
	draws[i,] <- sample(x, d, replace = FALSE, prob = w)
}

# compute empirical distributions
d_dist = matrix(rep(0, d*d), ncol = d)
d_counts = matrix(rep(0, d), ncol = d)
for(i in 1:d){
	d_dist[,i] <- colSums(draws == i) / n_draws
	d_counts[i] <- sum(draws == i)
}

# plot empirical distributions
dev.new(height = 16, width = 12)
par(mfrow=c(5,2))
for(i in 1:d){
	barplot(d_dist[i,], main = paste("sample index", i))
}

# plot actual weight distribution
dev.new(height = 4, width = 6)
par(mfrow=c(1,1))
barplot(w, main = "actual weights")

dev.new(height = 4, width = 6)
barplot(d_counts/sum(d_counts), main = "counts")

## END OF SCRIPT ------------------------------------------------------------------ #
# ################################################################################# #

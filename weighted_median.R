r <- c(rep(-7,10),-5,rep(-2,2),rep(1,3),rep(2,3),rep(5,3),4,rep(7,2),10)
r <- sort(c(rpois(1000, 25), rnorm(1000, -25, 10)))
w <- rle(r)$lengths
w_cs <- cumsum(w)
l <- sum(w)
u <- rle(r)$values
if(l %% 2 == 1){
	l_med <- (l + 1)/2
	r_med <- r[l_med]
	} else {
	l_med <- l/2
	r_med <- (r[l_med] + r[l_med+1])/2
}
median(r)
r_med

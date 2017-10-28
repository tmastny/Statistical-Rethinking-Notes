library(rethinking)
data("Howell1")
d <- Howell1
d2 <- d[d$age >= 18,]

# R 4.13
sample_mu <- rnorm(1e4, 178, 20)
sample_sigma <- runif(1e4, 0, 50)
prior_h <- rnorm(1e4, sample_mu, sample_sigma)
dens(prior_h)
simplehist(prior_h)

# R Code 4.25
flist <- alist(
    height ~ dnorm(mu, sigma),
    mu ~ dnorm(178, 20),
    sigma ~ dunif(0,50)
)
m4.1 <- map(flist, data=d2)
precis(m4.1)

# R Code 4.32
post <- extract.samples(m4.1, n=1e4)

# R Code 4.38
m4.3 <- map(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*weight,
    a ~ dnorm(178, 100),
    b ~ dnorm(0,10),
    sigma ~ dunif(0, 50)
  ),
  data=d2)

# R Code 4.42
d2$weight.c <- d2$weight - mean(d2$weight)

m4.4 <- map(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*weight.c,
    a ~ dnorm(178, 100),
    b ~ dnorm(0,10),
    sigma ~ dunif(0, 50)
  ),
  data=d2)

precis(m4.4, corr=TRUE)
# now mean of a is the expected height for the mean of weight

# R Code 4.45

plot(height~weight, data=d2, col='blue')
abline(a=coef(m4.3)['a'], b=coef(m4.3)['b'])

# R Code 4.46
post <- extract.samples(m4.3)

# R Code 4.48

N <- 352
dN <- d2[1:N,]
mN <- map(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*weight,
    a ~ dnorm(178, 100),
    b ~ dnorm(0,10),
    sigma ~ dunif(0, 50)
  ),
  data=dN)

post <- extract.samples(mN, n=20)

plot(dN$weight, dN$height, xlim=range(d2$weight), ylim=range(d2$height), col=rangi2)
mtext(concat("N = ", N))

for (i in 1:20)
  abline(a=post$a[i], b=post$b[i], col=col.alpha('black',0.3))

# R Code 4.54
# samples of various heights, given sequence of weights
weight.seq <- seq(from=25, to=70, by=1)
mu <- link(m4.3, data=data.frame(weight=weight.seq))
mu.mean <- apply(mu, 2, mean)
mu.HPDI <- apply(mu, 2, HPDI, prob=0.89)

plot(height~weight, data=d2, col=col.alpha(rangi2,0.5))
lines(weight.seq, mu.mean)
shade(mu.HPDI, weight.seq)

# R Code 4.59
# Simulated heights, accounting for uncertainity in
# likelihood and posterior
sim.height <- sim(m4.3, data=list(weight=weight.seq))
height.PI <- apply(sim.height, 2 , PI, prob=0.89)

plot(height~weight, data=d2, col=col.alpha(rangi2,0.5))
lines(weight.seq, mu.mean)
shade(mu.HPDI, weight.seq)
shade(height.PI, weight.seq)

# 4.5 Polynomial Regression
plot(height~weight, data=d)
ggplot(d, aes(x=weight, y=height, col=age>=18)) + geom_point()

# Standardizing (see the difference from centralizing)
d$weight.s <- (d$weight - mean(d$weight))/sd(d$weight)
ggplot(d, aes(x=weight.s, y=height, col=age>=18)) + geom_point()

# R Code 4.66
d$weight.s2 <- d$weight.s^2
m4.5 <- map(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b1*weight.s + b2*weight.s2,
    a ~ dnorm(178, 100),
    b1 ~ dnorm(0,10),
    b2 ~ dnorm(0,10),
    sigma ~ dunif(0, 50)
  ),
  data=d)
precis(m4.5)

# R Code 4.68
weight.seq <- seq(from=-2.2, to=2, length.out=30)
pred_dat <- list(weight.s=weight.seq, weight.s2=weight.seq^2)
mu <- link(m4.5, data=pred_dat)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI, prob=0.89)
sim.height <- sim(m4.5, data=pred_dat)
height.PI <- apply(sim.height, 2, PI, prob=0.89)

plot(height~weight.s, data=d, col=col.alpha(rangi2,0.5))
lines(weight.seq, mu.mean)
shade(mu.PI, weight.seq)
shade(height.PI, weight.seq)

# Plot with real values of predictor
plot(height~weight.s, data=d, col=col.alpha(rangi2,0.5), xaxt='n')
at <- c(-2,1,0,1,2)
labels <- at*sd(d$weight) + mean(d$weight)
axis(side=1, at=at, labels=round(labels,1))



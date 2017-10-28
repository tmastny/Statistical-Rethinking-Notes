# 4H1
library(rethinking)
data("Howell1")

d <- Howell1
d2 <- d[d$age >= 18,]

m4H1 <- map(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*weight,
    a ~ dnorm(178, 100),
    b ~ dnorm(0,10),
    sigma ~ dunif(0, 50)
  ),
  data=d2)

precis(m4H1)

weight.sim <- c(46.95, 43.72, 64.78, 32.59, 54.63)
height.sim <- sim(m4H1, data=list(weight=weight.sim))
height.PI <- apply(height.sim, 2 , PI, prob=0.89)
colMeans(height.sim)


# 4H2
d3 = d[d$age<18,]

m4H2 <- map(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*weight,
    a ~ dnorm(100, 100),
    b ~ dnorm(0,10),
    sigma ~ dunif(0, 50)
  ),
  data=d3)

weight.seq <- seq(from=0, to=50, by=1)
mu <- link(m4H2, data=data.frame(weight=weight.seq))
mu.mean <- apply(mu, 2, mean)
mu.HPDI <- apply(mu, 2, HPDI, prob=0.89)

sim.height <- sim(m4H2, data=list(weight=weight.seq))
height.HPDI <- apply(sim.height, 2 , HPDI, prob=0.89)

plot(height~weight, data=d3, col=col.alpha(rangi2,0.5))
lines(weight.seq, mu.mean)
shade(mu.HPDI, weight.seq)
shade(height.PI, weight.seq)

# for ten unit change in weight, how much change in height?
precis(m4H2)
10*coef(m4H2)['b']

# 4H3
library(rethinking)
data("Howell1")
d <- Howell1

m4H3 <- map(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + b*log2(weight),
    a ~ dnorm(178, 100),
    b ~ dnorm(0,140),
    sigma ~ dunif(0, 50)
  ),
  data=d)
precis(m4H3)

plot(height~weight, data=d, col=col.alpha(rangi2,0.5))

weight.seq <- seq(from=1, to=70, by=1)
mu <- link(m4H3, data=data.frame(weight=weight.seq))
mu.mean <- colMeans(mu)
mu.HPDI <- apply(mu, 2, HPDI, prob=0.89)

sim.height <- sim(m4H3, data=list(weight=weight.seq))
height.HPDI <- apply(sim.height, 2 , HPDI, prob=0.89)


lines(weight.seq, mu.mean)
shade(mu.HPDI, weight.seq)
shade(height.HPDI, weight.seq)

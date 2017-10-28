library(tidyverse)
library(rethinking)
data(foxes)
d <- foxes

# 5H1
# fit weight ~ area and weight ~ groupsize
plot(weight~area, data=d)
plot(weight~groupsize, data=d)

# weight~area

plot(weight~area.s, data=d)

d$area.s = (d$area - mean(d$area))/sd(d$area)
m.w.a <- map(
  alist(
    weight ~ dnorm(mu, sigma),
    mu <- a + bA * area.s,
    a ~ dnorm(0,10),
    bA ~ dnorm(0,10),
    sigma ~ dunif(0,10)
  ),
  data = d
)
precis(m.w.a)

a.seq <- seq(from=-3, to=3, length.out=100)
mu <- link(m.w.a, data=data.frame(area.s=a.seq))
mu.PI <-apply(mu, 2, PI)

# no apparent effect
plot(weight~area.s, data=d, col=rangi2)
abline(m.w.a)
shade(mu.PI, a.seq)

# with lm:
plot(weight~area.s, data=d)
abline(lm(weight~area.s, data=d))

# weight~groupsize
m.w.gs <- map(
  alist(
    weight ~ dnorm(mu, sigma),
    mu <- a + bGS * groupsize,
    a ~ dnorm(0,10),
    bGS ~ dnorm(0,10),
    sigma ~ dunif(0,10)
  ),
  data = d
)
precis(m.w.gs)

gs.seq <- seq(from=2, to=8, length.out=100)
mu <- link(m.w.gs, data=data.frame(groupsize=gs.seq))
mu.PI <-apply(mu, 2, PI)

# no apparent effect
plot(weight~groupsize, data=d, col=rangi2)
abline(m.w.gs)
shade(mu.PI, gs.seq)






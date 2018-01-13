library(rethinking)
data("WaffleDivorce")
d = WaffleDivorce

dlist = list(
  div_obs = d$Divorce,
  div_sd = d$Divorce.SE,
  R = d$Marriage,
  A = d$MedianAgeMarriage
)

m14.1 = map2stan(
  alist(
    div_est ~ dnorm(mu, sigma),
    mu <- a + bA * A + bR * R,
    div_obs ~ dnorm(div_est, div_sd),
    a ~ dnorm(0,100),
    bA ~ dnorm(0,10),
    bR ~ dnorm(0,10),
    sigma ~ dcauchy(0, 2.5)
  ),
  data = dlist,
  start = list(div_est = dlist$div_obs),
  WAIC = FALSE,
  iter = 5000,
  warmup = 1000,
  chains = 4,
  cores = 4,
  control = list(adapt_delta = 0.95)
)

m_old = map2stan(
  alist(
    Divorce ~ dnorm(mu, sigma),
    mu <- a + bA * MedianAgeMarriage + bR * Marriage,
    a ~ dnorm(0,100),
    bA ~ dnorm(0,10),
    bR ~ dnorm(0,10),
    sigma ~ dcauchy(0, 2.5)
  ),
  data = d,
  WAIC = FALSE,
  iter = 5000,
  warmup = 1000,
  chains = 4,
  cores = 4,
  control = list(adapt_delta = 0.95)
)
precis(m_old)

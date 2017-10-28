# 4M1
trials <- 1e4
mu.prior.samples <- rnorm(n=trials, mean=0, sd=10)
sigma.prior.samples <- runif(n=trials, min=0, max=10)
y.simulated.samples <- rnorm(n=trials, mean=mu.prior.samples, sd=sigma.prior.samples)
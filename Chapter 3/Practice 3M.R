n <- 15
w <- 8
p_grid <- seq(from=0, to=1, length.out = 1000)
prior <- ifelse(p_grid < 0.5, 0, 1)
likelihood <- dbinom(w, size=n, prob=p_grid)
posterior <- likelihood*prior
posterior <- posterior/sum(posterior)
samples <- sample(p_grid, prob=posterior, size = 1e4, replace=TRUE)
plot(posterior ~ p_grid)

# #we want a posterior predictive check for the model and data
prob_of_8_15 <- rbinom(1e4, size=n, prob=samples)
simplehist(prob_of_8_15)

prob_of_6_9 <- rbinom(1e4, size=9, prob=samples)
simplehist(prob_of_6_9)
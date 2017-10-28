library(rethinking)
data(homeworkch3)

b <- sum(birth1) + sum(birth2) #[1] 111. Total number of boys
n <- length(birth1) + length(birth2)

# ratio of boys increases in second birth
sum(birth1==1)/length(birth1) #[1] 0.51
sum(birth2==1)/length(birth2) #[1] 0.6

#[1] 0.555
(sum(birth1) + sum(birth2))/(length(birth1) + length(birth2))

# 3H1
p_grid <- seq(from=0, to=1, length.out=1000)
prior <- rep(1, length(p_grid))
# likelihood: for each possibly probably value in p_grid
# what is the probability that our data could be generated
likelihood <- dbinom(b, size=n, prob=p_grid)
posterior <- likelihood*prior
posterior <- posterior/sum(posterior)

p_grid[which.max(posterior)]

# 3H2
samples <- sample(p_grid, prob=posterior, size=1e4, replace=TRUE)
HPDI(samples, prob=c(0.5, 0.89, 0.97))

# 3H3
sim_births <- rbinom(1e4, size=n, prob=samples)
simplehist(sim_births)
abline(v=111, col='red')
mean(sim_births==111) #[1] 0.038
# Value in simulation with highest frequency
counts <- table(sim_births)
counts[which.max(counts)] # returns 112

# 3H4
num_boys_birth1 = sum(birth1)
birth_count <- 100
sim_births_low <- rbinom(1e4, size=birth_count, prob=samples)
simplehist(sim_births_low)
abline(v=num_boys_birth1, col='red')
# simulation is overestimating true proporation of boys in birth1
# likely that the each birth is not independent from last
# see proporation of boys in birth1 vs. birth2 above

# 3H5
# we want to compare simulated counts of boys, compared to
# the count of boys after the first born was a girl
# we need count of boys born after a girl
girls_birth1 <- length(birth1) - sum(birth1) # sum(1 - birth1)
boys_after_girls = sum(birth2[birth1==0])
sim_girl_births <- rbinom(1e4, size=girls_birth1, prob=samples)
simplehist(sim_girl_births)
abline(v=boys_after_girls, col='red') 
# the actual number of boys is inconsistent with the simulation




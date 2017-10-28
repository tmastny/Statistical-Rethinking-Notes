# 5M1
# Invent an example of a suprious correlation
# Think the correlation of marriage rate on divorce and
# median age and divorce rate. When considered independently, they
# both have an effect. With multiple regression, only 
# median age matters


# what is spurious correlation? When there two predictors
# show a positive correlation to an outcome, but when
# one predictor's effect is removed (by considering the 
# residuals after the regression) there is no longer a correlation
library(rethinking)

x1 <- seq(from=-5, to=5, length.out = 100)
y <- 100 + 2*x1 + rnorm(length(x1), 0, 15)
plot(y~x1)

d <- data.frame(y,x1)
m5M1 = lm(y~x1, data=d)

mu = coef(summary(m5M1))[1,1] + coef(summary(m5M1))[2,1]*d$x1
lines(x1,mu)
lines(x1, 100 + 2*x1, col='blue')
x1.resid = d$x1 - mu


# generate observations
trials <- 1e3
x.real.correlation <- rnorm(n = trials, mean = 1, sd = 1)
x.spurious.correlation <- rnorm(n = trials, mean = x.real.correlation)
y <- rnorm(n = trials, mean = x.real.correlation)

# plot correlations
df <- data.frame(y = y, x.real.correlation = x.real.correlation, x.spurious.correlation = x.spurious.correlation)
pairs(df)

# buid model, inspect results
model <- lm(y ~ x.real.correlation + x.spurious.correlation)
precis(model)


# 5M2
# Masked relation. Outcome variable should be correlated
# in opposite directions with each predictor. Likewise, the
# two predictors should be correlated

N <- 100
rho <- 0.7
x_pos <- rnorm(N)
x_neg <- rnorm(N, rho*x_pos, sqrt(1-rho^2))
y <- rnorm(N, x_pos - x_neg)
d <- data.frame(
  y,
  x_pos,
  x_neg
)
pairs(d)

# 5M3
# Does divorce predict marriage rates
library(rethinking)
data("WaffleDivorce")

d <- WaffleDivorce

m.D <- map(
  alist(
    Marriage ~ dnorm(mu, sigma),
    mu <- a + bD * Divorce,
    a ~ dnorm(10, 50),
    bD ~ dnorm(0, 5),
    sigma ~ dunif(0,10)
  ) , data=d
)


m.A <- map(
  alist(
    Marriage ~ dnorm(mu, sigma),
    mu <- a + bA * MedianAgeMarriage,
    a ~ dnorm(10, 50),
    bA ~ dnorm(0, 5),
    sigma ~ dunif(0,10)
  ) , data=d
)


m.A.D <- map(
  alist(
    Marriage ~ dnorm(mu, sigma),
    mu <- a + bA * MedianAgeMarriage + bD * Divorce,
    a ~ dnorm(10, 50),
    bA ~ dnorm(0, 5),
    bD ~ dnorm(0, 5),
    sigma ~ dunif(0,10)
  ) , data=d
)
d.pair <- data.frame(d$Marriage, d$MedianAgeMarriage, d$Divorce)
pairs(d.pair)
precis(m.D)
precis(m.A)
precis(m.A.D)

# We learn less about the marriage rate from divorce rates
# than we do from median age.

#Counterfactual plots
A.avg <- mean(d$MedianAgeMarriage)
D.seq <- seq(from=5, to=14, length.out = 30)
pred.data <- data.frame(
  Divorce = D.seq,
  MedianAgeMarriage=A.avg
)

#counterfactual mean divorce (mu)
mu <- link(m.A.D, data=pred.data)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

D.sim <- sim(m.A.D, data=pred.data, n=1e4)
D.PI <- apply(D.sim, 2, PI)

plot(Marriage~Divorce, data=d, type='n')
mtext("MedianAgeMarriage is average")
lines(D.seq, mu.mean)
shade(mu.PI, D.seq)
shade(D.PI, D.seq)

# an increase in divorce rate may result in an increase
# in marriage rate if divorce people get remarried.
# Recall that Marriage is marriage rate per 1000 adults
# and Divorce is divorce rate per 1000 adults.
# If people who get divorced remarry, then they will remarry 
# at a later age.
#
# To test this hypothesis, we want to see if states 
# with higher divorce rates also have a higher
# than expected median marriage age. 

plot(d$MedianAgeMarriage~d$Divorce)

m.a.on.d <- map(
  alist(
    MedianAgeMarriage ~ dnorm(mu, sigma),
    mu <- a + bD * Divorce,
    a ~ dnorm(10, 10),
    bD ~ dnorm(0, 10),
    sigma ~ dunif(0,10)
  ) , data=d
)
mu <- coef(m.a.on.d)['a'] + coef(m.a.on.d)['bD']*d$Divorce
m.a.resid <- d$MedianAgeMarriage - mu

plot(MedianAgeMarriage~Divorce, data=d, col=rangi2)
abline(m.a.on.d)
for (i in 1:length(m.a.resid)){
  y <- d$MedianAgeMarriage[i]
  x <- d$Divorce[i]
  lines(c(x,x), c(mu[i],y), lwd=0.5, col=col.alpha("black",0.7))
}
plot(d$Marriage~m.a.resid)
abline(lm(d$Marriage~m.a.resid))

# To the right: People who get married older than expected,
# given a state's divorce rate also tend to get married
# at a slower rate.

m.d.on.a <- map(
  alist(
    Divorce ~ dnorm(mu, sigma),
    mu <- a + bD * MedianAgeMarriage,
    a ~ dnorm(10, 50),
    bD ~ dnorm(0, 5),
    sigma ~ dunif(0,10)
  ) , data=d
)
mu <- coef(m.d.on.a)['a'] + coef(m.d.on.a)['bD']*d$MedianAgeMarriage
m.d.resid <- d$Divorce - mu

plot(Divorce~MedianAgeMarriage, data=d, col=rangi2)
abline(m.d.on.a)
for (i in 1:length(m.d.resid)){
  y <- d$Divorce[i]
  x <- d$MedianAgeMarriage[i]
  lines(c(x,x), c(mu[i],y), lwd=0.5, col=col.alpha("black",0.7))
}
plot(d$Marriage~m.d.resid)
abline(lm(d$Marriage~m.d.resid))

# On the Right: People who divorce more than expected
# based on their age do not tend to get married more frequently.

# Note: Without standardizing measures, extra care needs to go
# into choosing priors, because they won't be centered at 
# zero with a small SD

# 5M4
library(rethinking)
data("WaffleDivorce")
d <- WaffleDivorce

d$pct_LDS <- c(0.75, 4.53, 6.18, 1, 2.01, 2.82, 0.43, 0.55, 0.38,
               0.75, 0.82, 5.18, 26.35, 0.44, 0.66, 0.87, 1.25, 0.77, 0.64, 0.81,
               0.72, 0.39, 0.44, 0.58, 0.72, 1.14, 4.78, 1.29, 0.61, 0.37, 3.34,
               0.41, 0.82, 1.48, 0.52, 1.2, 3.85, 0.4, 0.37, 0.83, 1.27, 0.75,
               1.21, 67.97, 0.74, 1.13, 3.99, 0.92, 0.44, 11.5 )

# standardize variables
d$Marriage.standardized <- (d$Marriage - mean(d$Marriage)) / sd(d$Marriage)
d$MedianAgeMarriage.standardized <- (d$MedianAgeMarriage - mean(d$MedianAgeMarriage)) / sd(d$MedianAgeMarriage)
d$pct_LDS.standardized <- (d$pct_LDS - mean(d$pct_LDS)) / sd(d$pct_LDS)

# build a model, inspect results
model <- map(
  alist(
    Divorce ~ dnorm(mean = mu, sd = sigma),
    mu <- alpha + beta.marriage.rate * Marriage.standardized + beta.median.age.at.marriage * MedianAgeMarriage.standardized + beta.lds * pct_LDS.standardized,
    alpha ~ dnorm(mean = 0, sd = 100),
    c(beta.marriage.rate, beta.median.age.at.marriage, beta.lds) ~ dnorm(mean = 0, sd = 10),
    sigma ~ dunif(min = 0, 10)
  ),
  data = d
)
precis(model)
model.1 = lm(Divorce~Marriage.standardized +
               MedianAgeMarriage.standardized + 
               pct_LDS.standardized, data =d)

summary(model.1)
plot(d$Divorce~d$pct_LDS.standardized)
abline(lm(d$Divorce~d$pct_LDS.standardized))

# The bivariate fit of divorce regressing on pct_LDS does
# not seem to make a relevant prediction.
# But recall that multivariate regression is seeing if 
# the predictor is able to account for variable in the model
# outside the other predictor variables.
# Recall the large residuals from the model. According to the
# multivariate model, those large residuals shrink after 
# taking into account LDS states. 

m.a.on.lds <- map(
  alist(
    pct_LDS.standardized ~ dnorm(mean = mu, sd = sigma),
    mu <- alpha + beta.median.age.at.marriage * MedianAgeMarriage.standardized,
    alpha ~ dnorm(mean = 0, sd = 100),
    beta.median.age.at.marriage ~ dnorm(mean = 0, sd = 10),
    sigma ~ dunif(min = 0, 10)
  ),
  data = d
)
mu <- coef(m.a.on.lds)['alpha'] + coef(m.a.on.lds)['beta.median.age.at.marriage']*d$MedianAgeMarriage.standardized
lds.resid <- d$pct_LDS.standardized - mu

plot(pct_LDS.standardized~MedianAgeMarriage.standardized, data=d, col=rangi2)
abline(m.a.on.lds)
for (i in 1:length(lds.resid)){
  y <- d$pct_LDS.standardized[i]
  x <- d$MedianAgeMarriage.standardized[i]
  lines(c(x,x), c(mu[i],y), lwd=0.5, col=col.alpha("black",0.7))
}
plot(d$Divorce~lds.resid)
abline(lm(d$Divorce~lds.resid))
# Last plot shows that if the number of LDS members is 
# greater than expected given a median marriage age,
# they are less likely to divorce, which is what the model shows








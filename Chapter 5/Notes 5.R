library(rethinking)
data("WaffleDivorce")

d <- WaffleDivorce

d$MedianAgeMarriage.s <- (d$MedianAgeMarriage - mean(d$MedianAgeMarriage))/sd(d$MedianAgeMarriage)
d$Marriage.s <- (d$Marriage - mean(d$Marriage))/sd(d$Marriage)

m5.1 <- map(
  alist(
    Divorce ~ dnorm(mu, sigma),
    mu <- a + bA * MedianAgeMarriage.s,
    a ~ dnorm(10, 10),
    bA ~ dnorm(0, 10),
    sigma ~ dunif(0,10)
  ) , data=d
)

precis(m5.1)

MAM.seq <- seq(from=-3, to=3.5, length.out=30)
mu <- link(m5.1, data=data.frame(MedianAgeMarriage.s=MAM.seq))
mu.PI <-apply(mu, 2, PI)

plot(Divorce~MedianAgeMarriage.s, data=d, col=rangi2)
abline(m5.1)
shade(mu.PI, MAM.seq)


m5.3 <- map(
  alist(
    Divorce ~ dnorm(mu, sigma),
    mu <- a + bR * Marriage.s + bA * MedianAgeMarriage.s,
    a <- dnorm(10, 10),
    bR ~ dnorm(0 ,1),
    bA ~ dnorm(0, 1),
    sigma ~ dunif(0, 10)
  ),
  data=d
)

precis(m5.3)


m5.4 <- map(
  alist(
    Marriage.s ~ dnorm(mu, sigma),
    mu <- a + b * MedianAgeMarriage.s,
    a ~ dnorm(0 ,10),
    b ~ dnorm(0 , 1),
    sigma ~ dunif(0,10)
  ),
  data = d
)
mu <- coef(m5.4)['a'] + coef(m5.4)['b']*d$MedianAgeMarriage.s
m.resid <- d$Marriage.s - mu

plot(Marriage.s~MedianAgeMarriage.s, data=d, col=rangi2)
abline(m5.4)
for (i in 1:length(m.resid)){
  x <- d$MedianAgeMarriage.s[i]
  y <- d$Marriage.s[i]
  lines(c(x,x), c(mu[i],y), lwd=0.5, col=col.alpha("black",0.7))
}
plot(d$Divorce~m.resid)
abline(lm(d$Divorce~m.resid))

# residuals controlling for marriage rate
m5.4b <- map(
  alist(
    MedianAgeMarriage.s ~ dnorm(mu, sigma),
    mu <- a + b * Marriage.s,
    a ~ dnorm(0 ,10),
    b ~ dnorm(0 , 1),
    sigma ~ dunif(0,10)
  ),
  data = d
)
mu <- coef(m5.4b)['a'] + coef(m5.4b)['b']*d$Marriage.s
m.resid <- d$MedianAgeMarriage.s - mu

plot(MedianAgeMarriage.s~Marriage.s, data=d, col=rangi2)
abline(m5.4b)
for (i in 1:length(m.resid)){
  y <- d$MedianAgeMarriage.s[i]
  x <- d$Marriage.s[i]
  lines(c(x,x), c(mu[i],y), lwd=0.5, col=col.alpha("black",0.7))
}
plot(d$Divorce~m.resid)
abline(lm(d$Divorce~m.resid))


#Counterfactual plots
A.avg <- mean(d$MedianAgeMarriage.s)
R.seq <- seq(from=-3, to=3, length.out = 30)
pred.data <- data.frame(
  Marriage.s = R.seq,
  MedianAgeMarriage.s=A.avg
)

#counterfactual mean divorce (mu)
mu <- link(m5.3, data=pred.data)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

R.sim <- sim(m5.3, data=pred.data, n=1e4)
R.PI <- apply(R.sim, 2, PI)

plot(Divorce~Marriage.s, data=d, type='n')
mtext("MedianAgeMarriage.s=0")
lines(R.seq, mu.mean)
shade(mu.PI, R.seq)
shade(R.PI, R.seq)

# Posterior Prediction Plots
mu <- link(m5.3)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

divorce.sim <- sim(m5.3, n=1e4)
divorce.PI <- apply(divorce.sim, 2, PI)

plot(mu.mean~d$Divorce, col=rangi2, ylim=range(mu.PI), 
     xlab="Observed Divorce", ylab="Predicted Divorce")
abline(a = 0, b=1, lty=2)
for (i in 1:nrow(d)) {
  lines( rep(d$Divorce[i], 2), c(mu.PI[1,i], mu.PI[2,i]), 
         col=rangi2)
}
identify(x=d$Divorce, y=mu.mean, labels=d$Loc, cex=0.8)

# compute residuals
divorce.resid <- d$Divorce - mu.mean
o <- order(divorce.resid)
dotchart(divorce.resid[o], labels=d$Loc[o], xlim=c(-6,5), cex=0.5)
abline(v=0, col=col.alpha("black",0.2))
for (i in 1:nrow(d)) {
  j <- o[i]
  lines(d$Divorce[j] - c(mu.PI[1,j], divorce.PI[2,j]), rep(i,2))
  points(d$Divorce[j] - c(mu.PI[1,j], divorce.PI[2,j]), rep(i,2),
         pch=3, cex=0.6, col='gray')
}

d$WaffleHouses.per.cap = d$WaffleHouses/d$Population
plot(divorce.resid~d$WaffleHouses.per.cap)
abline(lm(divorce.resid~d$WaffleHouses.per.cap))

# spurious
N <- 100
x_real = rnorm(N)
x_spur <- rnorm(N, x_real)
y = rnorm(N, x_real)
d <- data.frame(y, x_real, x_spur)

pairs(d)
m1 = lm(y~x_real+x_spur, data=d)
m2 = lm(y~x_real, data=d)

# Chapter 5.2 Masked Relation
library(rethinking)
data(milk)
d <- milk
str(d)

# Want to find the relation between
# calories per gram of milk,
# mass of female
# percent of brain mass that is neocortex

dcc <- d[complete.cases(d),]
# Model where neocortex percentage predicts calories in milk
m5.5 <- map(
  alist(
    kcal.per.g ~ dnorm(mu, sigma),
    mu <- a + bn*neocortex.perc,
    a ~ dnorm(0, 100),
    bn ~ dnorm(0,1),
    sigma ~ dunif(0,1)
  ),
  data = dcc
)
precis(m5.5, digits=3)

np.seq <- 0:100
pred.data <- data.frame(neocortex.perc=np.seq)
mu <- link(m5.5, data=pred.data, n=1e4)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

plot(kcal.per.g~neocortex.perc, data=dcc, col=rangi2)
lines(np.seq, mu.mean)
lines(np.seq, mu.PI[1,], lty=2)
lines(np.seq, mu.PI[2,], lty=2)

# Predict calories based on log(mass)
# log mass corresponds to magnitude instead of measure
dcc$log.mass <- log(dcc$mass)

m5.6 <- map(
  alist(
    kcal.per.g ~ dnorm(mu, sigma),
    mu <- a + bn*log.mass,
    a ~ dnorm(0, 100),
    bn ~ dnorm(0,1),
    sigma ~ dunif(0,1)
  ),
  data = dcc
)
precis(m5.6, digits=3)

np.seq <- -5:5
pred.data <- data.frame(log.mass=np.seq)
mu <- link(m5.6, data=pred.data, n=1e4)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

plot(kcal.per.g~log.mass, data=dcc, col=rangi2)
lines(np.seq, mu.mean)
lines(np.seq, mu.PI[1,], lty=2)
lines(np.seq, mu.PI[2,], lty=2)

# Consider the multivariate model with both 
# neocortex percentage and female body mass
m5.7 <- map(
  alist(
    kcal.per.g ~ dnorm(mu, sigma),
    mu <- a + bn*neocortex.perc + bm*log.mass,
    a ~ dnorm(0, 100),
    bn ~ dnorm(0,1),
    bm ~ dnorm(0,1),
    sigma ~ dunif(0,1)
  ),
  data = dcc
)
precis(m5.7, digits=3)

# controlling for mass, how does this model predict
# calories from neocortex percentage
mean.log.mass <- mean(log(dcc$mass))
np.seq <- 0:100
pred.data <- data.frame(
  neocortex.perc=np.seq,
  log.mass=mean.log.mass
)
mu <- link(m5.7, data=pred.data, n=1e4)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

plot(kcal.per.g~neocortex.perc, data=dcc, type="n")
lines(np.seq, mu.mean)
lines(np.seq, mu.PI[1,], lty=2)
lines(np.seq, mu.PI[2,], lty=2)


# Controlling for neocortex percentage, how well does the
# model predict calories from female body mass
mean.np <- mean(dcc$neocortex.perc)
lgm.seq <- -5:5
pred.data <- data.frame(
  log.mass=lgm.seq,
  neocortex.perc=mean.np)
mu <- link(m5.7, data=pred.data, n=1e4)
mu.mean <- apply(mu, 2, mean)
mu.PI <- apply(mu, 2, PI)

plot(kcal.per.g~log.mass, data=dcc, type='n')
lines(lgm.seq, mu.mean)
lines(lgm.seq, mu.PI[1,], lty=2)
lines(lgm.seq, mu.PI[2,], lty=2)

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

# Chapter 5.3
# Multicollinearity, when predictors that are correlated
# produce a worse regression than when considered individually

N <- 100
height <- rnorm(N, 10, 2)
leg_prop <- runif(N, 0.4, 0.5)
leg_left <- leg_prop*height + rnorm(N, 0, 0.02)
leg_right <- leg_prop*height + rnorm(N, 0, 0.02)
d <- data.frame(
  height,
  leg_left,
  leg_right
)
pairs(d)
# as expected, since our model legs have a consistent proportion
# to total height, height is well correlated with both
# left and right leg length, which are also correlated with
# each other. For a given height, the right and left leg are the
# same, with a small error of 0.02

m5.8 <- map(
  alist(
    height ~ dnorm(mu, sigma),
    mu <- a + bl*leg_left + br*leg_right,
    a ~dnorm(10, 100),
    bl ~ dnorm(2, 10),
    br ~ dnorm(2, 10),
    sigma ~ dunif(0,10)
  ),
  data=d
)
precis(m5.8)

post <- extract.samples(m5.8)
plot(bl~br, post, col=col.alpha(rangi2,0.1), pch=16)

sum_blbr <- post$bl + post$br
dens(sum_blbr, col=rangi2, lwd=2, xlab='sum of bl and br')

# Multicollinear Milk
library(rethinking)
data(milk)
d <- milk

# calories regressed on percent fat
m5.10 <- map(
  alist(
    kcal.per.g ~ dnorm(mu, sigma),
    mu <- a + bf*perc.fat,
    a ~ dnorm(0.6, 10),
    bf ~ dnorm(0,1),
    sigma ~ dunif(0,10)
  ),
  data=d
)

# calories regressed on percent lactose
m5.11 <- map(
  alist(
    kcal.per.g ~ dnorm(mu, sigma),
    mu <- a + bf*perc.lactose,
    a ~ dnorm(0.6, 10),
    bf ~ dnorm(0,1),
    sigma ~ dunif(0,10)
  ),
  data=d
)
precis(m5.10, digits=3)
precis(m5.11, digits=3)

m5.12 <- map(
  alist(
    kcal.per.g ~ dnorm(mu, sigma),
    mu <- a + bf*perc.lactose + bl*perc.lactose,
    a ~ dnorm(0.6, 10),
    bf ~ dnorm(0,1),
    bl ~ dnorm(0,1),
    sigma ~ dunif(0,10)
  ),
  data=d
)
precis(m5.12)

pairs(~kcal.per.g + perc.fat + perc.lactose, data=d, col=rangi2)

# Post treatment bias
# number of plants
N <- 100
# initial heights, simulated
h0 <- rnorm(N, 10, 2)

# assign treatment and simulate fungus and growth
treatment <- rep(0:1, each=N/2)
fungus <- rbinom(N, size = 1, prob=0.5 - treatment*0.4)
h1 <- h0 + rnorm(N, 5 - 3*fungus)

d <- data.frame(
  h0=h0,
  h1=h1,
  treatment=treatment,
  fungus=fungus
)

m5.13 <- map(
  alist(
    h1 ~ dnorm(mu,sigma),
    mu <- a + bh*h0 + bt*treatment + bf*fungus,
    a ~ dnorm(0,100),
    c(bh,bt,bf) ~ dnorm(0,10),
    sigma ~ dunif(0,10)
  ),
  data=d
)
precis(m5.13)

# problem: multivariate models asks a very specific question:
# once we know whether or not a plant had fungus, 
# does the regression improve if we also know treatment?
# the answer is not really

m5.14 <- map(
  alist(
    h1 ~ dnorm(mu,sigma),
    mu <- a + bh*h0 + bt*treatment,
    a ~ dnorm(0,100),
    c(bh,bt) ~ dnorm(0,10),
    sigma ~ dunif(0,10)
  ),
  data=d
)
precis(m5.14)

# post-treatment bias: literally including a variable that
# changes as a result of another variable. Here the soil treatment
# implicitly predicted the amount of fungus, so a regression with
# both does not provide more info than one without


# Chapter 5.4 Categorical Variables
data("Howell1")
d <- Howell1
str(d)

m5.15 <- map(
  alist(
    height ~ dnorm(mu,sigma),
    mu <- a + bm*male,
    a ~ dnorm(178,100),
    bm ~ dnorm(0,10),
    sigma ~ dunif(0,50)
  ),
  data=d
)
precis(m5.15)

post <- extract.samples(m5.15)
mu.male <- post$a + post$bm
PI(mu.male)

# Multiple Categories
data(milk)
d <- milk
unique(d$clade)

d$clade.NWM <- ifelse(d$clade=="New World Monkey", 1, 0)
d$clade.OWM <- ifelse(d$clade=="Old World Monkey", 1, 0)
d$clade.S <- ifelse(d$clade=="Strepsirrhine", 1, 0)


m5.16 <- map(
  alist(
    kcal.per.g ~ dnorm(mu,sigma),
    mu <- a + b.NWM*clade.NWM + b.OWM*clade.OWM + b.S*clade.S,
    a ~ dnorm(0.6,10),
    b.NWM ~ dnorm(0,1),
    b.OWM ~ dnorm(0,1),
    b.S ~ dnorm(0,1),
    sigma ~ dunif(0,10)
  ),
  data=d
)
precis(m5.16)

post <- extract.samples(m5.16)
mu.ape <- post$a
mu.NWM <- post$a + post$b.NWM
mu.OWM <- post$a + post$b.OWM
mu.S <- post$a + post$b.S

precis(data.frame(mu.ape,mu.NWM,mu.OWM,mu.S))







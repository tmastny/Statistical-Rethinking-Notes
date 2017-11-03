library(brms)
library(rethinking)

data("chimpanzees")
d <- chimpanzees
d$recipient <- NULL

mod <- brm( pulled_left ~ 1 + (1 | actor) +
              prosoc_left*condition - condition,
            data = d, family = bernoulli(), iter = 5000,
            prior = c(set_prior("normal(0,10)", class = 'Intercept'),
                      set_prior("normal(0,10)", class = 'b'),
                      set_prior("cauchy(0,1)", class = 'sd',
                                group = 'actor')))


d.pred <- data.frame(prosoc_left = c(0,1,0,1),
           condition = c(0,0,1,1),
           actor = c(8,8,8,8))

# doesn't work
fitted(mod, newdata=d.pred, allow_new_levels = TRUE, 
       sample_new_levels = "gaussian")

  # Error in cor[, k] : subscript out of bounds
  # In addition: Warning messages:
  # 1: In min(cor) : no non-missing arguments to min; returning Inf
  # 2: In max(cor) : no non-missing arguments to max; returning -Inf

# however, this work
fitted(mod, newdata=d.pred, allow_new_levels = TRUE, 
       sample_new_levels = "uncertainty")

  #    Estimate Est.Error   2.5%ile  97.5%ile
  # 1 0.5116447 0.2771832 0.1981633 0.9965754
  # 2 0.6587207 0.2183250 0.3462795 0.9985525
  # 3 0.5116447 0.2771832 0.1981633 0.9965754
  # 4 0.6342009 0.2297492 0.3207487 0.9982886

# and this
fitted(mod, newdata=d.pred, allow_new_levels = TRUE, 
       sample_new_levels = "old_levels")

  #    Estimate  Est.Error   2.5%ile  97.5%ile
  # 1 0.3315611 0.05830818 0.2246771 0.4499873
  # 2 0.5275584 0.07385115 0.3822419 0.6702285
  # 3 0.3315611 0.05830818 0.2246771 0.4499873
  # 4 0.4943078 0.07544043 0.3456309 0.6399698


# These work, but I don't think they are what I want
# I want new level intercept calculated from normal(0, sd(Intercept))
# which I thought sample_new_levels = 'gaussian' did
fitted(mod, newdata=d.pred, re_formula = NA, allow_new_levels = TRUE, 
       sample_new_levels = 'gaussian')

  #    Estimate Est.Error   2.5%ile  97.5%ile
  # 1 0.5885639 0.1818542 0.1991871 0.9108210
  # 2 0.7450321 0.1534072 0.3610278 0.9601266
  # 3 0.5885639 0.1818542 0.1991871 0.9108210
  # 4 0.7216129 0.1601085 0.3377181 0.9536617

d.pred2 <- data.frame(prosoc_left = c(0,1,0,1),
                     condition = c(0,0,1,1),
                     actor = c(1,1,1,1))

# gives same result as above
fitted(mod, newdata=d.pred2, re_formula = NA)

  #    Estimate Est.Error   2.5%ile  97.5%ile
  # 1 0.5885639 0.1818542 0.1991871 0.9108210
  # 2 0.7450321 0.1534072 0.3610278 0.9601266
  # 3 0.5885639 0.1818542 0.1991871 0.9108210
  # 4 0.7216129 0.1601085 0.3377181 0.9536617


# For example, if I calculate the population level effect by hand
# I get the same results as fitted, which tells me they do not
# include the sd(Intercept)
sample <- fixef(mod, summary=FALSE)

linear_link <- function(prosoc_left, condition, sample) {
  logodds <- sample[,'Intercept'] + sample[,'prosoc_left'] * prosoc_left +
    sample[,'prosoc_left:condition'] * prosoc_left * condition
  return(logistic(logodds))
}
prosoc_left <- c(0,1,0,1)
condition <- c(0,0,1,1)
pred.table <- sapply(1:4, function(i) linear_link(prosoc_left[i], condition[i], sample))
pred.mean <- apply(pred.table, 2, mean)

  # 0.5885639 0.7450321 0.5885639 0.7216129

pred.PI <- apply(pred.table, 2, PI, prob=0.95)

  #          [,1]      [,2]      [,3]      [,4]
  # 3%  0.1991871 0.3610278 0.1991871 0.3377181
  # 98% 0.9108210 0.9601266 0.9108210 0.9536617



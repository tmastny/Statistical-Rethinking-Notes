# Chapter 12 Notes
Tim  
10/22/2017  



## 12.1 Multilevel Tadpoles


```r
library(rethinking)
library(brms)
library(tidyverse)
library(tidybayes)

rstan_options (auto_write=TRUE)
options (mc.cores=parallel::detectCores ()) # Run on multiple cores

data(reedfrogs)
d <- reedfrogs
d$tank <- 1:NROW(d)
d %>% as.tibble()
```

```
## # A tibble: 48 x 6
##    density   pred   size  surv propsurv  tank
##      <int> <fctr> <fctr> <int>    <dbl> <int>
##  1      10     no    big     9      0.9     1
##  2      10     no    big    10      1.0     2
##  3      10     no    big     7      0.7     3
##  4      10     no    big    10      1.0     4
##  5      10     no  small     9      0.9     5
##  6      10     no  small     9      0.9     6
##  7      10     no  small    10      1.0     7
##  8      10     no  small     9      0.9     8
##  9      10   pred    big     4      0.4     9
## 10      10   pred    big     9      0.9    10
## # ... with 38 more rows
```

To fit the multilevel model in `brms` as described in we need to explicitly remove the population parameter with `-1` as shown below:


```r
mod.intercept <- brm(surv | trials(density) ~ -1 + (1 | tank),
                     family = binomial(), data=d,
                     prior = c(set_prior("normal(0,1)", class = 'sd',
                                         group = 'tank', coef='Intercept'),
                               set_prior("cauchy(0,1)", class = 'sd',
                                         group = 'tank')))
```

```r
summary(mod.intercept)
```

```
##  Family: binomial(logit) 
## Formula: surv | trials(density) ~ -1 + (1 | tank) 
##    Data: d (Number of observations: 48) 
## Samples: 4 chains, each with iter = 2000; warmup = 1000; thin = 1; 
##          total post-warmup samples = 4000
##     ICs: LOO = NA; WAIC = NA; R2 = NA
##  
## Group-Level Effects: 
## ~tank (Number of levels: 48) 
##               Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
## sd(Intercept)     2.03      0.25     1.61     2.56       1045 1.00
## 
## Samples were drawn using sampling(NUTS). For each parameter, Eff.Sample 
## is a crude measure of effective sample size, and Rhat is the potential 
## scale reduction factor on split chains (at convergence, Rhat = 1).
```

How do we interpret the summary of this hierarchal model? Well, in terms of the rethinking model on page 359, `sd(Intercept)` is the adaptive standard deviation from the Normal distrubtion that each tank intercept is draw from. In other words, each intercept is draw from a normal distribution with `sd(Intercept)` standard deviation with an adaptive mean. 

Let's reconstruct the plot on page 361:


```r
p <- ggplot(d, aes(x=tank)) +
  geom_point(aes(y=propsurv)) +
  geom_hline(yintercept = mean(d$propsurv), linetype=2) +
  geom_vline(xintercept = c(16.5, 32.5))
```

Now let's take a look at the predicted survival rates:


```r
d.mean <- d %>%
  add_fitted_samples(mod.intercept) %>%
  mean_qi() %>%
  mutate(propsurv_pred = estimate/density)
```

`brms` doesn't give the adapative population prior $\alpha$ directly. We need to calculate it by sampling from the posterior parameters of the group-level intercepts. The function `tidybayes::spread_samples` is nice tool to easily sample parameters from the posterior into tidy data frames. `r_tank` is the name of the group-level intercepts. There is a standard convention to `brms` name, but you can always find the names with `parnames` as shown below. `tidybayes::spread_samples` is flexible enough to allow syntax matching for the parameter of interest. 


```r
parnames(mod.intercept)
```

```
##  [1] "sd_tank__Intercept"   "r_tank[1,Intercept]"  "r_tank[2,Intercept]" 
##  [4] "r_tank[3,Intercept]"  "r_tank[4,Intercept]"  "r_tank[5,Intercept]" 
##  [7] "r_tank[6,Intercept]"  "r_tank[7,Intercept]"  "r_tank[8,Intercept]" 
## [10] "r_tank[9,Intercept]"  "r_tank[10,Intercept]" "r_tank[11,Intercept]"
## [13] "r_tank[12,Intercept]" "r_tank[13,Intercept]" "r_tank[14,Intercept]"
## [16] "r_tank[15,Intercept]" "r_tank[16,Intercept]" "r_tank[17,Intercept]"
## [19] "r_tank[18,Intercept]" "r_tank[19,Intercept]" "r_tank[20,Intercept]"
## [22] "r_tank[21,Intercept]" "r_tank[22,Intercept]" "r_tank[23,Intercept]"
## [25] "r_tank[24,Intercept]" "r_tank[25,Intercept]" "r_tank[26,Intercept]"
## [28] "r_tank[27,Intercept]" "r_tank[28,Intercept]" "r_tank[29,Intercept]"
## [31] "r_tank[30,Intercept]" "r_tank[31,Intercept]" "r_tank[32,Intercept]"
## [34] "r_tank[33,Intercept]" "r_tank[34,Intercept]" "r_tank[35,Intercept]"
## [37] "r_tank[36,Intercept]" "r_tank[37,Intercept]" "r_tank[38,Intercept]"
## [40] "r_tank[39,Intercept]" "r_tank[40,Intercept]" "r_tank[41,Intercept]"
## [43] "r_tank[42,Intercept]" "r_tank[43,Intercept]" "r_tank[44,Intercept]"
## [46] "r_tank[45,Intercept]" "r_tank[46,Intercept]" "r_tank[47,Intercept]"
## [49] "r_tank[48,Intercept]" "lp__"
```

```r
# group parameter samples using tidybayes
pop.intercept <- mod.intercept %>% spread_samples(r_tank[tank,])
pop.proportion <- logistic(mean(pop.intercept$r_tank))

p +
  geom_point(aes(y=propsurv_pred), data=d.mean, shape=1) +
  geom_hline(yintercept = pop.proportion) # predicted population mean (intercept)
```

![](Chapter_12_Notes_files/figure-html/unnamed-chunk-6-1.png)<!-- -->

Here we can see the classic shrinkage on each tank using a multilevel model. Moreover the tanks are sorted from smallest to larger: therefore as we move to the right we have less shrinkage. Lastly, also note that the new estimated population mean is different than the observed mean.

As an aside, let's try to fit the model with a predictor. McElreath withholds varying slopes until the Chapter 13, so let's try a population predictor. I would expect that predators would decrease the the probability of survival.


```r
# indicator variable for predator
d$pred <- ifelse(d$pred == 'pred', 1, 0)
mod.pred <- brm(surv | trials(density) ~ pred + (1 | tank), data=d,
                family = binomial(),
                prior = c(set_prior("normal(0,1)", class = 'sd',
                                     group = 'tank', coef='Intercept'),
                           set_prior("cauchy(0,1)", class = 'sd',
                                     group = 'tank')))
```

```r
summary(mod.pred)
```

```
##  Family: binomial(logit) 
## Formula: surv | trials(density) ~ pred + (1 | tank) 
##    Data: d (Number of observations: 48) 
## Samples: 4 chains, each with iter = 2000; warmup = 1000; thin = 1; 
##          total post-warmup samples = 4000
##     ICs: LOO = NA; WAIC = NA; R2 = NA
##  
## Group-Level Effects: 
## ~tank (Number of levels: 48) 
##               Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
## sd(Intercept)     0.83      0.15     0.58     1.15       1825 1.00
## 
## Population-Level Effects: 
##           Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
## Intercept     2.71      0.25     2.25     3.23       3475 1.00
## pred         -2.68      0.32    -3.30    -2.07       2894 1.00
## 
## Samples were drawn using sampling(NUTS). For each parameter, Eff.Sample 
## is a crude measure of effective sample size, and Rhat is the potential 
## scale reduction factor on split chains (at convergence, Rhat = 1).
```

As expected, the presence of a predator has a strong effect on the survival of tadpoles.

It would also makes sense that survival depends on the number of tadpoles relative to the size of a tank, whether a predator is present or not. If there are many tadpoles in a small tank, the survival rate should decrease when a predator is present.


```r
mod.interaction <- brm(surv | trials(density) ~ pred*density*size
                       + (1 | tank),
                       data=d, family = binomial(),
                       prior = c(set_prior("normal(0,1)", class = 'sd',
                                     group = 'tank', coef='Intercept'),
                                 set_prior("cauchy(0,1)", class = 'sd',
                                     group = 'tank')))
```

```r
summary(mod.interaction)
```

```
##  Family: binomial(logit) 
## Formula: surv | trials(density) ~ pred * density * size + (1 | tank) 
##    Data: d (Number of observations: 48) 
## Samples: 4 chains, each with iter = 2000; warmup = 1000; thin = 1; 
##          total post-warmup samples = 4000
##     ICs: LOO = NA; WAIC = NA; R2 = NA
##  
## Group-Level Effects: 
## ~tank (Number of levels: 48) 
##               Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
## sd(Intercept)     0.71      0.15     0.45     1.02       1307 1.00
## 
## Population-Level Effects: 
##                        Estimate Est.Error l-95% CI u-95% CI Eff.Sample
## Intercept                  2.41      0.97     0.63     4.40       1100
## pred                      -1.31      1.19    -3.63     0.93        912
## density                    0.01      0.04    -0.06     0.08       1060
## sizesmall                  0.03      1.39    -2.66     2.79        852
## pred:density              -0.08      0.04    -0.16     0.01        885
## pred:sizesmall             0.14      1.73    -3.22     3.44        835
## density:sizesmall         -0.01      0.05    -0.11     0.09        844
## pred:density:sizesmall     0.04      0.06    -0.08     0.16        828
##                        Rhat
## Intercept              1.00
## pred                   1.00
## density                1.00
## sizesmall              1.00
## pred:density           1.00
## pred:sizesmall         1.00
## density:sizesmall      1.00
## pred:density:sizesmall 1.00
## 
## Samples were drawn using sampling(NUTS). For each parameter, Eff.Sample 
## is a crude measure of effective sample size, and Rhat is the potential 
## scale reduction factor on split chains (at convergence, Rhat = 1).
```

```r
LOO(mod.intercept, mod.pred, mod.interaction)
```

```
##                                  LOOIC    SE
## mod.intercept                   232.31  6.48
## mod.pred                        213.46 11.32
## mod.interaction                 214.79 11.27
## mod.intercept - mod.pred         18.85  8.58
## mod.intercept - mod.interaction  17.52  7.99
## mod.pred - mod.interaction       -1.33  5.40
```

Based on the LOO information criteria, the predator predictor model seems like the best fit. 

We can also do some posterior predictor checks on the response density, as advocated by Gelman et al in Chapter 6 of Bayesian Data Analysis. `pp_check` is a method in `brms` that calls the package `bayesplot`, a tie in Stan for visualizing the posterior:


```r
pp_check(mod.pred)
```

![](Chapter_12_Notes_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

```r
pp_check(mod.intercept)
```

![](Chapter_12_Notes_files/figure-html/unnamed-chunk-12-2.png)<!-- -->

A standard warning with these checks, as noted in the `pp_check` `brms` documentation, a graphical fit may look good for both models. Indeed, here both seem to fit okay, with `mod.pred` being a bit better as expected. Information criteria like LOO help us select the model in light of clear graphical errors.

## Multilevel chimps

Next we return to the chimp data and consider multiple cluster types.


```r
data("chimpanzees")
d <- chimpanzees
d$recipient <- NULL
d %>% as.tibble()
```

```
## # A tibble: 504 x 7
##    actor condition block trial prosoc_left chose_prosoc pulled_left
##    <int>     <int> <int> <int>       <int>        <int>       <int>
##  1     1         0     1     2           0            1           0
##  2     1         0     1     4           0            0           1
##  3     1         0     1     6           1            0           0
##  4     1         0     1     8           0            1           0
##  5     1         0     1    10           1            1           1
##  6     1         0     1    12           1            1           1
##  7     1         0     2    14           1            0           0
##  8     1         0     2    16           1            0           0
##  9     1         0     2    18           0            1           0
## 10     1         0     2    20           0            1           0
## # ... with 494 more rows
```

#### One Cluster

First, we'll fit one cluster:


```r
mod <- brm( pulled_left ~ 1 + (1 | actor) +
                      prosoc_left*condition - condition,
                    data = d, family = bernoulli(),
                    prior = c(set_prior("normal(0,10)", class = 'Intercept'),
                              set_prior("normal(0,10)", class = 'b'),
                              set_prior("cauchy(0,1)", class = 'sd',
                                        group = 'actor')))
```

```r
summary(mod)
```

```
##  Family: bernoulli(logit) 
## Formula: pulled_left ~ 1 + (1 | actor) + prosoc_left * condition - condition 
##    Data: d (Number of observations: 504) 
## Samples: 4 chains, each with iter = 2000; warmup = 1000; thin = 1; 
##          total post-warmup samples = 4000
##     ICs: LOO = NA; WAIC = NA; R2 = NA
##  
## Group-Level Effects: 
## ~actor (Number of levels: 7) 
##               Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
## sd(Intercept)     2.27      0.95     1.10     4.53        835 1.01
## 
## Population-Level Effects: 
##                       Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
## Intercept                 0.43      0.96    -1.53     2.39        891 1.00
## prosoc_left               0.82      0.26     0.32     1.32       2028 1.00
## prosoc_left:condition    -0.13      0.30    -0.71     0.48       2466 1.00
## 
## Samples were drawn using sampling(NUTS). For each parameter, Eff.Sample 
## is a crude measure of effective sample size, and Rhat is the potential 
## scale reduction factor on split chains (at convergence, Rhat = 1).
```

And to get the get the total intercepts for each actor as per R Code 12.22, we use `brms:coef`, which is the sum of the population and group level effects per level.


```r
coef(mod)$actor[,,'Intercept']
```

```
##     Estimate Est.Error    2.5%ile   97.5%ile
## 1 -0.7151930 0.2710388 -1.2610062 -0.2096185
## 2  4.6289749 1.6508541  2.5556927  8.7364989
## 3 -1.0254793 0.2786557 -1.5779986 -0.4861515
## 4 -1.0184635 0.2739890 -1.5746687 -0.4992832
## 5 -0.7134285 0.2622040 -1.2502631 -0.2140437
## 6  0.2272997 0.2732106 -0.2982452  0.7733302
## 7  1.7638911 0.3983002  1.0368972  2.6044048
```

Alternatively, we can use `tidybayes`. One reason to prefer `tidybayes` is that it has consistent `tidyverse` style syntax and always outputs tidy tibbles, grouped by `spread_sample` parameters for quick summaries.


```r
mod %>% 
  spread_samples(r_actor[actor,], b_Intercept) %>%
  mean_qi(r_actor + b_Intercept) # no group_by necessary, already included
```

```
## # A tibble: 7 x 5
## # Groups:   actor [7]
##   actor `r_actor + b_Intercept`   conf.low  conf.high .prob
##   <int>                   <dbl>      <dbl>      <dbl> <dbl>
## 1     1              -0.7151930 -1.2610062 -0.2096185  0.95
## 2     2               4.6289749  2.5556927  8.7364989  0.95
## 3     3              -1.0254793 -1.5779986 -0.4861515  0.95
## 4     4              -1.0184635 -1.5746687 -0.4992832  0.95
## 5     5              -0.7134285 -1.2502631 -0.2140437  0.95
## 6     6               0.2272997 -0.2982452  0.7733302  0.95
## 7     7               1.7638911  1.0368972  2.6044048  0.95
```

#### Two Clusters

The study was organized into different blocks, where each monkey pulled their levels once per day as opposed to one monkey doing all their pulls at once. This technique called cross-classification is a useful design feature to eliminate temporal effects on the treatment.

Thus we can also provide unique intercepts for each blocks. Ideally, we want to see that there is little to no variation within each blot: that's the entire design goal of the blocks. If there is added variation in different blocks, we can measure that variation and see if the treatment appears after controlling for the block variation.


```r
mod.cluster <- brm(pulled_left ~ 1 + (1 | actor) + (1 | block) + 
                     prosoc_left + prosoc_left:condition,
                   data=d, family=bernoulli(),
                   prior = c(set_prior("normal(0,10)", class = 'Intercept'),
                             set_prior("normal(0,10)", class = 'b'),
                             set_prior("cauchy(0,1)", class = 'sd',
                                       group = 'actor')))
```

```r
summary(mod.cluster)
```

```
##  Family: bernoulli(logit) 
## Formula: pulled_left ~ 1 + (1 | actor) + (1 | block) + prosoc_left + prosoc_left:condition 
##    Data: d (Number of observations: 504) 
## Samples: 4 chains, each with iter = 2000; warmup = 1000; thin = 1; 
##          total post-warmup samples = 4000
##     ICs: LOO = NA; WAIC = NA; R2 = NA
##  
## Group-Level Effects: 
## ~actor (Number of levels: 7) 
##               Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
## sd(Intercept)     2.26      0.90     1.12     4.71        769 1.01
## 
## ~block (Number of levels: 6) 
##               Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
## sd(Intercept)     0.24      0.23     0.01     0.78       1482 1.00
## 
## Population-Level Effects: 
##                       Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
## Intercept                 0.45      1.00    -1.40     2.46        544 1.00
## prosoc_left               0.83      0.27     0.31     1.36       2452 1.00
## prosoc_left:condition    -0.14      0.30    -0.75     0.44       3015 1.00
## 
## Samples were drawn using sampling(NUTS). For each parameter, Eff.Sample 
## is a crude measure of effective sample size, and Rhat is the potential 
## scale reduction factor on split chains (at convergence, Rhat = 1).
```
These results match the output of R Code 12.24. 

For the charter in Figure 12.4:


```r
# it would be nice if spread_samples() with no args just spread every
# parameter available. 
parnames(mod.cluster)
```

```
##  [1] "b_Intercept"             "b_prosoc_left"          
##  [3] "b_prosoc_left:condition" "sd_actor__Intercept"    
##  [5] "sd_block__Intercept"     "r_actor[1,Intercept]"   
##  [7] "r_actor[2,Intercept]"    "r_actor[3,Intercept]"   
##  [9] "r_actor[4,Intercept]"    "r_actor[5,Intercept]"   
## [11] "r_actor[6,Intercept]"    "r_actor[7,Intercept]"   
## [13] "r_block[1,Intercept]"    "r_block[2,Intercept]"   
## [15] "r_block[3,Intercept]"    "r_block[4,Intercept]"   
## [17] "r_block[5,Intercept]"    "r_block[6,Intercept]"   
## [19] "lp__"
```

```r
mod.cluster %>%
  gather_samples(r_actor[actor,], r_block[block,],
                 b_Intercept, b_prosoc_left, `b_prosoc_left:condition`,
                 sd_block__Intercept, sd_actor__Intercept) %>%
  mean_qi() %>%
  replace_na(list(actor = "", block = "")) %>%
  unite(variable, term, actor, block) %>%
  ggplot(aes(y = variable, x = estimate)) +
  geom_point() +
  geom_segment(aes(x=conf.low, xend=conf.high, yend=variable))
```

![](Chapter_12_Notes_files/figure-html/unnamed-chunk-20-1.png)<!-- -->

And the LOO comparison:


```r
LOO(mod, mod.cluster)
```

```
##                    LOOIC    SE
## mod               531.64 19.52
## mod.cluster       533.13 19.73
## mod - mod.cluster  -1.49  1.81
```





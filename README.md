# baconr: An rstan implementation of the Bayesian Age-Depth model *Bacon* (Blaauw and Christen, 2011).

-------------------------------

**baconr** implements the Bayesian Age~Depth model, **Bacon**, in the **Stan** probabilistic programming language. It is at a very early stage of development and at the moment implements only the core non-Gaussian AR1 model described in Blaauw and Christen (2011). Functions from the R package **Bchron** (Parnell 2016) can be used to calibrate ^14^C ages to calendar ages.


There are currently just three exported functions: 

* `make_stan_dat` prepares data and parameter values into the correct format to be passed to the Stan model 
* `stan_bacon` calls a pre-compiled Stan implementation of Bacon and estimates the age model
* `plot_stan_bacon` plots a sample from the posterior distribution of the estimated age model

**baconr** also includes the example data, core MSB2K, included in the existing C++ implementation of [Bacon v2.2](http://www.chrono.qub.ac.uk/blaauw/bacon.html)

The motivation for creating this package is to make use of the Bacon age modelling routine more flexible, and to make further development 
of the age model itself easier, by coding it in a widely used higher-level probabilistic programming language.


*  Parnell, Andrew. 2016. Bchron: Radiocarbon Dating, Age-Depth Modelling, Relative Sea Level Rate Estimation, and Non-Parametric Phase Modelling. R package version 4.2.6. https://CRAN.R-project.org/package=Bchron

*  Blaauw, Maarten, and J. Andrés Christen. 2011. Flexible Paleoclimate Age-Depth Models Using an Autoregressive Gamma Process. Bayesian Analysis 6 (3): 457-74. doi:10.1214/ba/1339616472.

*  Stan Development Team. 2016. Stan Modeling Language Users Guide and Reference Manual, Version 2.14.0.   http://mc-stan.org




## Installation

**baconr** can be installed directly from Github


```r
if (!require("devtools")) {
  install.packages("devtools")
}

devtools::install_github("andrewdolman/baconr")

opts_chunk$set(cache = TRUE, autodep = TRUE)
```


## Development

We aim to develop this package following the [Guidelines for R packages providing interfaces to Stan](https://cran.r-project.org/web/packages/rstantools/vignettes/developer-guidelines.html)


## Example



```r
library(baconr)
library(knitr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(Bchron)

opts_chunk$set(echo=TRUE, message = FALSE, warning = FALSE, cache = TRUE,
               fig.width = 6, fig.pos = "H", dpi = 300, autodep = TRUE)
```

### Data and parameters

First convert ^14^C ages to calendar ages. Functions from the R package **Bchron** can be used for this.


```r
cal.ages <- Bchron::BchronCalibrate(ages=MSB2K$age,
                            ageSds=MSB2K$error,
                            calCurves=rep("intcal13", nrow(MSB2K)),
                            ids=paste0("Date-", 1:nrow(MSB2K)))

age.samples <- Bchron::SampleAges(cal.ages)

MSB2K$age.cal <- apply(age.samples, 2, median)
MSB2K$age.cal.sd <- apply(age.samples, 2, sd)
```


the output from `make_stan_dat` is shown here to illustrate the required data and parameter format, but a separate call to `make_stan_dat` is not normally necessary as it is used internally by `stan_bacon`


```r
# Get number of sections K, so that they will be ~ 5cm
K_for_5cm <- round(diff(range(MSB2K$depth)) / 5)

stan_dat <- make_stan_dat(depth = MSB2K$depth, 
  obs_age = MSB2K$age.cal, 
  obs_err = MSB2K$age.cal.sd,
  K = K_for_5cm, nu = 6,
  acc_mean = 20, acc_alpha = 1.5,
  mem_mean = 0.1, mem_strength = 4)

stan_dat
```

```
## $depth
##  [1]  1.5  4.5  8.5 12.5 14.5 14.5 14.5 17.5 20.5 21.5 21.5 22.5 28.5 31.5
## [15] 32.5 33.5 34.5 37.5 38.5 41.5 43.5 46.5 47.5 48.5 49.5 50.5 52.5 53.5
## [29] 54.5 55.5 58.5 59.5 64.5 70.5 71.5 73.5 75.5 77.5 79.5 99.5
## 
## $obs_age
##  [1] 4665.0 4635.0 4538.0 4713.0 4590.0 4639.0 4621.0 4711.0 4737.0 4851.0
## [11] 4960.0 5146.0 5104.0 5354.0 5404.0 5480.0 5389.0 5526.0 5487.0 5575.0
## [21] 5651.0 5622.0 5730.0 5811.0 5734.0 5846.0 5784.0 6019.0 5918.0 5861.0
## [31] 6013.5 6077.0 6080.0 6155.0 6300.0 6367.0 6316.0 6337.0 6431.0 6707.0
## 
## $obs_err
##  [1] 100.81027 104.52430 110.08822  84.04567 112.45944 104.45101 105.65239
##  [8]  79.41174  85.20993 104.20485 114.31259 108.72639 116.54580 126.51788
## [15]  94.55640  85.34468 117.13311  83.24009  93.85686  75.56124  76.94347
## [22]  66.29263  85.49959  73.66715  83.26338  98.04250  73.31915  89.25769
## [29]  75.20711  86.46860  90.11399  63.40706  75.86692  82.40875  75.74652
## [36]  50.19286  55.69356  45.54537  85.53695  50.63419
## 
## $K
## [1] 20
## 
## $nu
## [1] 6
## 
## $acc_mean
## [1] 20
## 
## $acc_alpha
## [1] 1.5
## 
## $mem_mean
## [1] 0.1
## 
## $mem_strength
## [1] 4
## 
## $N
## [1] 40
## 
## $acc_beta
## [1] 0.075
## 
## $mem_alpha
## [1] 0.4
## 
## $mem_beta
## [1] 3.6
## 
## $c
##  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
## 
## $delta_c
## [1] 5.22
## 
## $c_depth_bottom
##  [1]   5.22  10.44  15.66  20.88  26.10  31.32  36.54  41.76  46.98  52.20
## [11]  57.42  62.64  67.86  73.08  78.30  83.52  88.74  93.96  99.18 104.40
## 
## $c_depth_top
##  [1]  0.00  5.22 10.44 15.66 20.88 26.10 31.32 36.54 41.76 46.98 52.20
## [12] 57.42 62.64 67.86 73.08 78.30 83.52 88.74 93.96 99.18
## 
## $which_c
##  [1]  1  1  2  3  3  3  3  4  4  5  5  5  6  7  7  7  7  8  8  8  9  9 10
## [24] 10 10 10 11 11 11 11 12 12 13 14 14 15 15 15 16 20
```




### Fit the bacon model with `stan_bacon`


```r
fit <- stan_bacon(
  depth = MSB2K$depth, 
  obs_age = MSB2K$age.cal, 
  obs_err = MSB2K$age.cal.sd,
  K = K_for_5cm, nu = 6,
  acc_mean = 20, acc_alpha = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  iter = 2000, chains = 4)
```

```
## 
## SAMPLING FOR MODEL 'bacon' NOW (CHAIN 1).
## 
## Gradient evaluation took 0 seconds
## 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
## Adjust your expectations accordingly!
## 
## 
## Iteration:    1 / 2000 [  0%]  (Warmup)
## Iteration:  200 / 2000 [ 10%]  (Warmup)
## Iteration:  400 / 2000 [ 20%]  (Warmup)
## Iteration:  600 / 2000 [ 30%]  (Warmup)
## Iteration:  800 / 2000 [ 40%]  (Warmup)
## Iteration: 1000 / 2000 [ 50%]  (Warmup)
## Iteration: 1001 / 2000 [ 50%]  (Sampling)
## Iteration: 1200 / 2000 [ 60%]  (Sampling)
## Iteration: 1400 / 2000 [ 70%]  (Sampling)
## Iteration: 1600 / 2000 [ 80%]  (Sampling)
## Iteration: 1800 / 2000 [ 90%]  (Sampling)
## Iteration: 2000 / 2000 [100%]  (Sampling)
## 
##  Elapsed Time: 1.966 seconds (Warm-up)
##                1.353 seconds (Sampling)
##                3.319 seconds (Total)
## 
## 
## SAMPLING FOR MODEL 'bacon' NOW (CHAIN 2).
## 
## Gradient evaluation took 0 seconds
## 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
## Adjust your expectations accordingly!
## 
## 
## Iteration:    1 / 2000 [  0%]  (Warmup)
## Iteration:  200 / 2000 [ 10%]  (Warmup)
## Iteration:  400 / 2000 [ 20%]  (Warmup)
## Iteration:  600 / 2000 [ 30%]  (Warmup)
## Iteration:  800 / 2000 [ 40%]  (Warmup)
## Iteration: 1000 / 2000 [ 50%]  (Warmup)
## Iteration: 1001 / 2000 [ 50%]  (Sampling)
## Iteration: 1200 / 2000 [ 60%]  (Sampling)
## Iteration: 1400 / 2000 [ 70%]  (Sampling)
## Iteration: 1600 / 2000 [ 80%]  (Sampling)
## Iteration: 1800 / 2000 [ 90%]  (Sampling)
## Iteration: 2000 / 2000 [100%]  (Sampling)
## 
##  Elapsed Time: 1.992 seconds (Warm-up)
##                1.368 seconds (Sampling)
##                3.36 seconds (Total)
## 
## 
## SAMPLING FOR MODEL 'bacon' NOW (CHAIN 3).
## 
## Gradient evaluation took 0 seconds
## 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
## Adjust your expectations accordingly!
## 
## 
## Iteration:    1 / 2000 [  0%]  (Warmup)
## Iteration:  200 / 2000 [ 10%]  (Warmup)
## Iteration:  400 / 2000 [ 20%]  (Warmup)
## Iteration:  600 / 2000 [ 30%]  (Warmup)
## Iteration:  800 / 2000 [ 40%]  (Warmup)
## Iteration: 1000 / 2000 [ 50%]  (Warmup)
## Iteration: 1001 / 2000 [ 50%]  (Sampling)
## Iteration: 1200 / 2000 [ 60%]  (Sampling)
## Iteration: 1400 / 2000 [ 70%]  (Sampling)
## Iteration: 1600 / 2000 [ 80%]  (Sampling)
## Iteration: 1800 / 2000 [ 90%]  (Sampling)
## Iteration: 2000 / 2000 [100%]  (Sampling)
## 
##  Elapsed Time: 1.875 seconds (Warm-up)
##                1.253 seconds (Sampling)
##                3.128 seconds (Total)
## 
## 
## SAMPLING FOR MODEL 'bacon' NOW (CHAIN 4).
## 
## Gradient evaluation took 0 seconds
## 1000 transitions using 10 leapfrog steps per transition would take 0 seconds.
## Adjust your expectations accordingly!
## 
## 
## Iteration:    1 / 2000 [  0%]  (Warmup)
## Iteration:  200 / 2000 [ 10%]  (Warmup)
## Iteration:  400 / 2000 [ 20%]  (Warmup)
## Iteration:  600 / 2000 [ 30%]  (Warmup)
## Iteration:  800 / 2000 [ 40%]  (Warmup)
## Iteration: 1000 / 2000 [ 50%]  (Warmup)
## Iteration: 1001 / 2000 [ 50%]  (Sampling)
## Iteration: 1200 / 2000 [ 60%]  (Sampling)
## Iteration: 1400 / 2000 [ 70%]  (Sampling)
## Iteration: 1600 / 2000 [ 80%]  (Sampling)
## Iteration: 1800 / 2000 [ 90%]  (Sampling)
## Iteration: 2000 / 2000 [100%]  (Sampling)
## 
##  Elapsed Time: 1.962 seconds (Warm-up)
##                1.285 seconds (Sampling)
##                3.247 seconds (Total)
```



```r
options(width = 85)
print(fit$fit, par = c("R", "w", "c_ages"))
```

```
## Inference for Stan model: bacon.
## 4 chains, each with iter=2000; warmup=1000; thin=1; 
## post-warmup draws per chain=1000, total post-warmup draws=4000.
## 
##               mean se_mean    sd    2.5%     25%     50%     75%   97.5% n_eff Rhat
## R             0.72    0.00  0.17    0.28    0.62    0.76    0.86    0.93  2796    1
## w             0.28    0.00  0.21    0.00    0.09    0.24    0.45    0.69  2447    1
## c_ages[1]  4567.59    0.79 49.99 4463.54 4535.53 4570.24 4602.12 4659.68  4000    1
## c_ages[2]  4606.96    0.68 43.32 4515.18 4579.19 4607.88 4636.32 4689.70  4000    1
## c_ages[3]  4651.76    0.64 40.22 4571.15 4625.35 4651.46 4678.88 4732.42  4000    1
## c_ages[4]  4712.63    0.65 40.97 4634.73 4683.69 4711.58 4739.78 4792.30  4000    1
## c_ages[5]  4857.44    0.87 54.78 4754.81 4819.81 4856.49 4895.04 4963.61  4000    1
## c_ages[6]  5097.94    1.54 87.40 4913.94 5044.28 5100.84 5154.81 5269.12  3242    1
## c_ages[7]  5329.73    0.82 51.83 5228.53 5294.78 5329.17 5364.54 5432.69  4000    1
## c_ages[8]  5478.37    0.67 42.25 5392.39 5450.83 5478.97 5505.40 5562.13  4000    1
## c_ages[9]  5594.06    0.60 37.74 5515.20 5570.12 5595.35 5619.85 5665.06  4000    1
## c_ages[10] 5709.27    0.60 38.16 5630.72 5684.06 5710.70 5735.19 5780.46  4000    1
## c_ages[11] 5848.49    0.62 39.12 5771.03 5822.33 5849.25 5874.10 5924.19  4000    1
## c_ages[12] 5980.32    0.63 39.60 5901.90 5954.62 5980.48 6006.89 6056.17  4000    1
## c_ages[13] 6085.54    0.71 44.60 5997.77 6055.65 6085.18 6114.97 6172.67  4000    1
## c_ages[14] 6184.60    0.75 47.20 6084.85 6154.66 6187.19 6217.79 6270.47  4000    1
## c_ages[15] 6296.17    0.49 31.29 6233.76 6276.39 6296.59 6316.98 6355.62  4000    1
## c_ages[16] 6372.68    0.51 32.16 6312.76 6350.59 6372.07 6393.71 6437.64  4000    1
## c_ages[17] 6457.16    0.81 51.50 6363.12 6421.48 6454.96 6490.09 6566.69  4000    1
## c_ages[18] 6539.56    0.95 60.06 6427.78 6499.49 6538.99 6578.53 6659.82  4000    1
## c_ages[19] 6623.03    0.96 60.53 6500.78 6583.89 6624.27 6661.41 6746.27  4000    1
## c_ages[20] 6707.74    0.86 54.66 6606.31 6673.15 6705.39 6740.18 6824.48  4000    1
## c_ages[21] 6807.96    1.32 83.73 6669.56 6750.57 6799.48 6855.56 6995.30  4000    1
## 
## Samples were drawn using NUTS(diag_e) at Fri Oct 13 21:59:07 2017.
## For each parameter, n_eff is a crude measure of effective sample size,
## and Rhat is the potential scale reduction factor on split chains (at 
## convergence, Rhat=1).
```


### Plot the estimated age model


```r
set.seed(20170406)
plot_stan_bacon(fit, 1000)
```

![](readme_files/figure-html/bacon_defaults-1.png)<!-- -->


### Fit again with stronger prior on accumulation rates



```r
fit2 <- stan_bacon(
  depth = MSB2K$depth, 
  obs_age = MSB2K$age.cal, 
  obs_err = MSB2K$age.cal.sd,
  K = K_for_5cm, nu = 6,
  acc_mean = 20, acc_alpha = 25,
  mem_mean = 0.7, mem_strength = 4,
  iter = 2000, chains = 4)
```



```r
plot_stan_bacon(fit2, 1000)
```

![](readme_files/figure-html/stronger_prior-1.png)<!-- -->


### Fit again with strong prior for lower memory



```r
fit3 <- stan_bacon(
  depth = MSB2K$depth, 
  obs_age = MSB2K$age.cal, 
  obs_err = MSB2K$age.cal.sd,
  K = K_for_5cm, nu = 6,
  acc_mean = 20, acc_alpha = 1.5,
  mem_mean = 0.1, mem_strength = 4,
  iter = 2000, chains = 4)
```



```r
plot_stan_bacon(fit3, 1000)
```

![](readme_files/figure-html/weaker_prior-1.png)<!-- -->


### Add tephras at 25 and 90 cm




```r
teph <- data.frame(depth = c(25, 87.5), age.cal = c(5000, 6500), age.cal.sd = c(3, 5))
MSB2K.2 <- bind_rows(MSB2K, teph) %>% 
  arrange(age.cal)
  

fit3 <- stan_bacon(
  depth = MSB2K.2$depth, 
  obs_age = MSB2K.2$age.cal, 
  obs_err = MSB2K.2$age.cal.sd,
  K = K_for_5cm, nu = 6,
  acc_mean = 20, acc_alpha = 1.5,
  mem_mean = 0.7, mem_strength = 4,
  iter = 2000, chains = 4)
```



```r
plot_stan_bacon(fit3, 1000) +
  geom_pointrange(data = teph,
                  aes(x = depth, y = age.cal, ymax = age.cal + age.cal.sd, ymin = age.cal - age.cal.sd),
                  group = NA, colour = "Blue", alpha = 0.5) 
```

![](readme_files/figure-html/tephras-1.png)<!-- -->


### Distribution of sediment accumulation rates


```r
age.fit <- fit
age.mod <- rstan::extract(age.fit$fit)

bp.dat <- age.mod$x[1:4000,] %>% 
  as_tibble() %>% 
  tibble::rownames_to_column("Rep") %>% 
  gather(Depth, value, -Rep) %>% 
  mutate(Depth = as.numeric(gsub("V", "", Depth)),
         Depth = Depth * age.fit$data$delta_c,
         Rep = as.numeric(Rep))

# bp.dat.bw <- bp.dat %>% 
#   group_by(Depth) %>% 
#   summarise(mean.val = mean(value),
#             upr = quantile(value, 0.9),
#             lwr = quantile(value, 0.1))
# 
# bp.dat.bw %>% 
#   ggplot(aes(x = Depth)) + 
#   geom_pointrange(aes(y = mean.val, ymax = upr, ymin = lwr))
```


```r
bp.dat %>% 
  ggplot(aes(x = Depth, y = 1/value, group = Depth)) + 
  geom_violin() +
  #geom_boxplot() +
  #geom_point(alpha = 0.15)
  scale_x_continuous("Depth [cm]")+
  scale_y_continuous("Sediment accumulation rate [cm/yr]",
                     trans = "log10", breaks = c(0.01, 0.05, 0.1, 0.5)) +
  annotation_logticks(sides = "l") +
  expand_limits(y = 0.01) + 
  theme_bw()
```

![](readme_files/figure-html/acc_rates-1.png)<!-- -->


```r
bp.dat %>% 
  ggplot(aes(x = Depth, y = value, group = Depth)) + 
  geom_violin() +
  #geom_boxplot() +
  #geom_point(alpha = 0.15)
  scale_x_continuous("Depth [cm]")+
  scale_y_continuous("Sedimentation interval [yr/cm]",
                     trans = "log10", breaks = c(1, 5, 10, 20, 50)) +
  annotation_logticks(sides = "l") +
#  expand_limits(y = 0.01) + 
  theme_bw()
```

![](readme_files/figure-html/acc_intervals-1.png)<!-- -->







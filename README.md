# Bias in occupancy estimate for a static model

Here we provide some R code to calculate bias in occupancy estimate as a function of the detection probability given various levels of occupancy probability, various number of sites and surveys.

Load package `unmarked` to carry out occupancy analyses
```r
library(unmarked)
```

Define function to carry out simulations
```r
occu_par <- function(
  nb_sites = 50, # number of sites
  nb_surveys = 5, # number of surveys
  occpr = 0.3, # occupancy prob 
  detpr = 0.5, # detection prob
  n_sim = 500){ # number of simulations

# preallocate memory for storing occupancy estimates
res <- rep(NA, n_sim)

# simulate data from a static occupancy model n_sim times
for (j in 1:n_sim){
  
  # define state process
  z <- rbinom(nb_sites, 1, occpr) # occupancy state
  
  # pre-allocate memory for matrix of detection/non-detections
  y <- matrix(NA, nrow = nb_sites, ncol = nb_surveys) # detection histories
  
  # define observation process
  for(i in 1:nb_sites){
      prob <- z[i] * detpr
      y[i,1:nb_surveys] <- rbinom(nb_surveys, 1, prob)
  }
  
  # format data
  dat <- unmarkedFrameOccu(y)
  
  # fit static occupancy model w/ constant parameters
  fm <- occu(~ 1 ~ 1, dat)
  
  # get estimate of occupancy prob
  res[j] <- backTransform(fm, type = 'state')@estimate
  
}

# return bias
bias <- round(mean((res - occpr)/occpr) * 100,1)
print(bias)
}
```

Grid on detection
```r
det <- seq(0.1, 0.9, by = 0.1)
```

Consider several scenarios, and carry out simulations (this is ugly):
```r
# nsites = 50, nsurveys = 5, occ = 0.2
sim1 <- rep(NA, length(det))
index <- 1
for (s in det){
  sim1[index] <- occu_par(nb_sites = 50, 
         nb_surveys = 5, 
         occpr = 0.2, 
         detpr = s,
         n_sim = 100)
  index <- index + 1
}
```

```
## [1] 175.6
## [1] 29.6
## [1] 5.2
## [1] 5.2
## [1] -0.1
## [1] -0.9
## [1] -2.2
## [1] 2.6
## [1] -0.3
```

```r
# nsites = 50, nsurveys = 5, occ = 0.5
det <- seq(0.1, 0.9, by = 0.1)
sim2 <- rep(NA, length(det))
index <- 1
for (s in det){
  sim2[index] <- occu_par(nb_sites = 50, 
                          nb_surveys = 5, 
                          occpr = 0.5, 
                          detpr = s,
                          n_sim = 100)
  index <- index + 1
}
```

```
## [1] 22
## [1] 7.5
## [1] 2.2
## [1] 0.6
## [1] 1.9
## [1] -1
## [1] 1.1
## [1] 0.8
## [1] -0.3
```

```r
# nsites = 100, nsurveys = 5, occ = 0.2
det <- seq(0.1, 0.9, by = 0.1)
sim3 <- rep(NA, length(det))
index <- 1
for (s in det){
  sim3[index] <- occu_par(nb_sites = 100, 
                          nb_surveys = 5, 
                          occpr = 0.2, 
                          detpr = s,
                          n_sim = 100)
  index <- index + 1
}
```

```
## [1] 79
## [1] 15.9
## [1] 5.9
## [1] 0.2
## [1] 2.9
## [1] 1.2
## [1] -0.3
## [1] -0.3
## [1] 0.6
```

```r
# nsites = 100, nsurveys = 5, occ = 0.5
det <- seq(0.1, 0.9, by = 0.1)
sim4 <- rep(NA, length(det))
index <- 1
for (s in det){
  sim4[index] <- occu_par(nb_sites = 100, 
                          nb_surveys = 5, 
                          occpr = 0.5, 
                          detpr = s,
                          n_sim = 100)
  index <- index + 1
}
```

```
## [1] 9.6
## [1] 3.6
## [1] 1.5
## [1] -1.8
## [1] 0.8
## [1] -0.6
## [1] 0.3
## [1] -0.6
## [1] 0
```

```r
# nsites = 100, nsurveys = 10, occ = 0.2
det <- seq(0.1, 0.9, by = 0.1)
sim5 <- rep(NA, length(det))
index <- 1
for (s in det){
  sim5[index] <- occu_par(nb_sites = 100, 
                          nb_surveys = 10, 
                          occpr = 0.2, 
                          detpr = s,
                          n_sim = 100)
  index <- index + 1
}
```

```
## [1] 6.3
## [1] 2.1
## [1] 1.2
## [1] 1.2
## [1] -0.1
## [1] -2
## [1] -0.7
## [1] 3.2
## [1] -0.2
```

```r
# nsites = 100, nsurveys = 10, occ = 0.5
det <- seq(0.1, 0.9, by = 0.1)
sim6 <- rep(NA, length(det))
index <- 1
for (s in det){
  sim6[index] <- occu_par(nb_sites = 100, 
                          nb_surveys = 10, 
                          occpr = 0.5, 
                          detpr = s,
                          n_sim = 100)
  index <- index + 1
}
```

```
## [1] 4.2
## [1] -0.3
## [1] 0.9
## [1] -1.1
## [1] -0.4
## [1] 0.8
## [1] 0.2
## [1] 0.7
## [1] -2
```

```r
# nsites = 50, nsurveys = 10, occ = 0.2
det <- seq(0.1, 0.9, by = 0.1)
sim7 <- rep(NA, length(det))
index <- 1
for (s in det){
  sim7[index] <- occu_par(nb_sites = 50, 
                          nb_surveys = 10, 
                          occpr = 0.2, 
                          detpr = s,
                          n_sim = 100)
  index <- index + 1
}
```

```
## [1] 39.3
## [1] 5.6
## [1] 6.3
## [1] -0.6
## [1] 1.5
## [1] 4.6
## [1] -3.1
## [1] -0.4
## [1] -2
```

```r
# nsites = 50, nsurveys = 10, occ = 0.5
det <- seq(0.1, 0.9, by = 0.1)
sim8 <- rep(NA, length(det))
index <- 1
for (s in det){
  sim8[index] <- occu_par(nb_sites = 50, 
                          nb_surveys = 10, 
                          occpr = 0.5, 
                          detpr = s,
                          n_sim = 100)
  index <- index + 1
}
```

```
## [1] 9.4
## [1] 1.1
## [1] -1.8
## [1] 1.1
## [1] 0.1
## [1] -1
## [1] 1.5
## [1] -0.1
## [1] 1.5
```

Visualize bias
```r
library(tidyverse)
res <- tibble(bias = c(sim1,sim2,sim3,sim4,sim5,sim6,sim7,sim8),
       nsites = as_factor(c(rep('50 sites', 18), rep('100 sites', 36), rep('50 sites', 18))),
       nsurveys = as_factor(c(rep('5 surveys', 36), rep('10 surveys', 36))),
       procc = as_factor(rep(c(rep(0.2, 9) , rep(0.5, 9)), 4)),
       prdet = rep(as_factor(det), 8))

res %>%
  ggplot() +
  aes(x = prdet, y = bias, fill = procc) + 
  geom_col(position = "dodge", width=.6) + 
  labs(x = 'detection probability',
       y = '% bias in occupancy estimate') + 
    scale_fill_manual(name = NULL,
                      values = c('blue', 'red'),
                      labels = c('Pr occupancy is 0.2','Pr occupancy is 0.5')) +
  facet_wrap(~ nsites + nsurveys, labeller = label_wrap_gen(multi_line = FALSE)) + 
  theme_bw(base_size = 16)
```

![](biasoccupancy.png)<!-- -->



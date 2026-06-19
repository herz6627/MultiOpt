## Overview

`MultiOpt` is an R package for running simulated annealing to select individuals best suited to the users priorities..

The package was developed to support analyses presented in:

> Author A, Author B, and Author C. (Year). *Title of manuscript*. Journal Name. DOI: xxx


## Installation

### Development version

```r
install.packages("remotes")

remotes::install_github("herz6627/MultiOpt")


```


## Getting Started

Load the package:

```r
library(MultiOpt)
```

Basic example:

```r

set.seed(12345)
n = 100
x = rnorm(n = n, mean = 120, sd = 2)
y = x * 3 + rnorm(n = n, mean = 0, sd = 20)
dat = data.frame(x = x, y = y)

trait_list = list(
x = as.matrix(dat$x),
y = as.matrix(dat$y)
)

measure_list = list(
 x = weighted_mean_of_vector,
 y = weighted_mean_of_vector
)

args_list = list(
 x = NULL,
 y = list(direction = -1)
)

trait_list_scaled <- scale_traits(trait_list)

sa_args = list(
trait_list = trait_list_scaled,
measure_list = measure_list,
measure_args_list = args_list,
n_t = 10,
weights_max = 1,
max_steps = 50000,
nda = T,
nd_samples = 500
)

test = rand_multiopt(n_runs = 5, multiopt_args = sa_args, parallel = F)
str(test)
```

---

## Overview

`MultiOpt` is an R package for running simulated annealing to select individuals 
best suited to the users priorities. This package has been designed as a tool
for ex situ conservation collections (seed banks, living collections, introductions, etc.)
but is applicable across many different fields. 

The package was developed to support analyses presented in:

> Author A, Author B, and Author C. (Year). *Title of manuscript*. Journal Name. DOI: xxx


Website: 
[Here](https://herz6627.github.io/MultiOpt/index.html)


## Installation

This R package is very similar to and utilizes the core functionality of the 
R package OptGenMix (jasongbragg/OptGenMix). MultiOpt expands OptGenMix to allow
for more than two traits to be considered at once, which can be helpful when
considering multiple forms of variation for your ex situ collection. Some of the 
publications OptGenMix has been used in include:

> Bragg, Jason G., et al. "Optimizing the genetic composition of a translocation population: Incorporating constraints and conflicting objectives." Molecular Ecology Resources 20.1 (2020): 54-65.

> Dimon, Richard J., et al. "Integrating Genomic and Climate Data to Design Representative Seed Production Areas: A Pragmatic Workflow for Climate‐Adjusted Provenancing." Ecology and Evolution 16.1 (2026): e72658.

> Bragg, Jason G., et al. "Conserving the genetic diversity of condemned populations: Optimizing collections and translocation." Evolutionary Applications 14.5 (2021): 1225-1238.


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



Generate some trait data to use

```r
set.seed(12345)
n = 100
x = rnorm(n = n, mean = 120, sd = 2)
y = x * 3 + rnorm(n = n, mean = 0, sd = 20)
dat = data.frame(x = x, y = y)
```


Format trait data to a list of matrices

```r
trait_list = list(
x = as.matrix(dat$x),
y = as.matrix(dat$y)
)

```


Here is where we specify how we want to measure success. This needs to be a list of functions as listed in ?multiopt_sa().
Full details are can be found in the Measures vignette. 

```r
measure_list = list(
 x = weighted_mean_of_vector,
 y = weighted_mean_of_vector
)
```

And we can specify and additional arguments to the selected functions in `measure_list`.
Some of the functions require these additional arguments, which can be found in the respective 
function descriptions, but a key thing to note is that this is the location to specify which 
direction you want to optimize (maximize or minimize). 
By default, MultiOpt maximizes a trait (direction = 1). 
Use direction = -1 to get the minimum instead as we do here for trait `y`.
```r
args_list = list(
 x = NULL,
 y = list(direction = -1)
)
```
Note: across these three lists, all trait names need to match for list elements.



It is important for simulated annealing that all traits being optimized are on the same scale. 
If trait values are not scaled, larger values will be prioritized by the algorithm.
trait_list_scaled() will scale all provided traits to 0-1.

```r
trait_list_scaled <- scale_traits(trait_list)
```


Simulated annealing is a stochastic algorithm. As such, it is good to run the algorithm a few times to get an idea of how stable the results are, which we can do with `rand_multiopt()` which requires a list of inputs to specify our simulated annealing parameters. `multiopt_sa()` will run a single round of simulated annealing and has slightly different formatting. 

```r
sa_args = list(
trait_list = trait_list_scaled, # our list of traits
measure_list = measure_list, # our list of measures
measure_args_list = args_list, # our list of arguments for each measure
n_t = 10, # how many individuals we want in our final population
weights_max = 1, # the maximum number of times an individual can be selected
max_steps = 5000 # how many time steps to run each simulated annealing
)
```


And now we can run simulated annealing
```r
test = rand_multiopt(
n_runs = 5, # how many times to replicate the algorithm
multiopt_args = sa_args# all of our arguments
)
str(test)
```

---

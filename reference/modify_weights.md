# Modify weight vector

Randomly selects a weight and moves it to a new individual, if possible,
while preserving limitations.

## Usage

``` r
modify_weights(w, w_min = NULL, w_max = NULL, n_t)
```

## Arguments

- w:

  Integer vector of current weights.

- w_min:

  Optional integer vector of minimum allowable weights. May be:

  - NULL (defaults to 0).

  - a scalar applied to all individuals.

  - a vector of length n_g (number of individuals with trait data)
    indicating minimum weight for each individual.

- w_max:

  Optional integer vector of maximum allowable weights. May be:

  - NULL (unconstrained)

  - a scalar applied to all individuals.

  - a vector of length n_g (number of individuals with trait data)

- n_t:

  Integer. Total target population size.

## Value

Modified weight vector of same length as w.

# Internal helper for validating weight constraints.

Internal helper for validating weight constraints.

## Usage

``` r
.check_weights(n_g, n_t = NULL, w_min, w_max)
```

## Arguments

- n_g:

  Integer. Number of genotypes/individuals.

- n_t:

  Integer. Total target population size.

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

## Value

List of length 2 containing a vector of length n_g for both w_min and
w_max

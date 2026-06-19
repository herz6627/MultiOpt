# Generate initial weights

Generates initial weights for use before initiating simulated annealing.

## Usage

``` r
propose_initial_weights(n_g, n_t, w_max = NULL, w_min = NULL, verbose = FALSE)
```

## Arguments

- n_g:

  Integer. Number of genotypes/individuals.

- n_t:

  Integer. Total target population size.

- w_max:

  Optional integer vector of maximum allowable weights. May be:

  - NULL (unconstrained)

  - a scalar applied to all individuals.

  - a vector of length n_g (number of individuals with trait data)

- w_min:

  Optional integer vector of minimum allowable weights. May be:

  - NULL (defaults to 0).

  - a scalar applied to all individuals.

  - a vector of length n_g (number of individuals with trait data)
    indicating minimum weight for each individual.

- verbose:

  Logical. Print progress messages.

## Value

Integer vector of weights

# Determine annealing temperature

Determine annealing temperature

## Usage

``` r
temp_scheduler(s, max_steps, max_t, min_t = 0, nda = FALSE)
```

## Arguments

- s:

  Integer. Current time step.

- max_steps:

  Numeric scalar. Maximum allowed time steps.

- max_t:

  Numeric scalar. Maximum allowed temperature.

- min_t:

  Numeric scalar. Minimum allowed temperature. Only relevant when nda =
  TRUE.

- nda:

  Logical. If TRUE, indicates to hold temp at minimum temperature to
  allow exploration of Pareto front.

## Value

A numeric scalar of current annealing temperature.

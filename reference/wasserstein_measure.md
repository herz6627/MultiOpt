# Calculate the Wasserstein Distance for use in simulated annealing

Uses the
[`?wasserstein`](https://herz6627.github.io/MultiOpt/reference/wasserstein.md)
to calculate the Wasserstein distance (\\W_p\\) between the original
trait data distribution and the proposed subset of individuals. Assumes
that v provides the total population distribution of trait variation.

## Usage

``` r
wasserstein_measure(v, w, ...)
```

## Arguments

- v:

  single-column matrix of numeric values. Will be coerced to a vector.

- w:

  Numeric vector of individual weights with same length as \`v\`.

- ...:

  Optional additional arguments. See
  [`?wasserstein`](https://herz6627.github.io/MultiOpt/reference/wasserstein.md)
  for specific options.

## Value

A single numeric value representing the p-Wasserstein distance
multiplied by -1.

## Details

As \`MultiOpt\` aims to maximize measures but smaller \\W_p\\ are
better, we multiply the \\W_p\\ by -1 to allow for minimizing the
differences in distributions.

## Examples

``` r
v <- c(1, 2, 3)
w <- c(0, 1, 1)
wasserstein_measure(v, w)
#> Error in wasserstein_measure(v, w): v must be a matrix.
```

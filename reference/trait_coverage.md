# Calculate Trait Coverage

Calculates the proportion of trait-value bins represented in a selected
subset relative to the bins represented in the full dataset.

## Usage

``` r
trait_coverage(v, w, n_bins = 10)
```

## Arguments

- v:

  single-column matrix of numeric values. Will be coerced to a vector.

- w:

  Numeric vector of individual weights with same length as \`v\`. to the
  observations in v. Values with a count of zero are excluded from the
  selected subset.

- n_bins:

  The number of equally spaced bins used to divide the trait range.

## Value

A numeric value representing the proportion of trait-value bins covered
by the selected subset. A value of 1 indicates that the subset covers
all bins represented in the full dataset.

## Details

The range of trait values in v is divided into n_bins equally spaced
bins. The function counts how many bins contain observations in the full
dataset and how many contain observations in the selected subset,
defined by the selection counts in w. The returned value is the
proportion of occupied bins in the subset relative to the full dataset.

This function does not measure how well the proposed subset matches the
distribution of trait data. For this type of measurement,
[`wasserstein_measure()`](https://herz6627.github.io/MultiOpt/reference/wasserstein_measure.md)
would be a better fit.

## Examples

``` r
v <- as.matrix(c(1, 2, 3))
w <- c(0, 1, 2)

trait_coverage(v, w, n_bins = 10)
#> Warning: n_bins > length(v), check that this makes sense.
#> Warning: n_bins > sum(w), check that this makes sense,
#>                               as output will be artifically low due to forced blank bins.
#> [1] 0.6666667
```

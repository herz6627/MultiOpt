# \#' Compute weighted sum of squared deviations from a target value

Calculates the weighted sum of squared differences between a numeric
vector (\`v\`) and a displacement value (\`disp\`). Each squared
deviation is weighted by a corresponding value in \`w\`.This could be
used, e.g., to minimize the mean difference between temperature of
origin for each sample, and temperate of a site.

## Usage

``` r
sum_of_squared_difference(v, w, disp = 0, direction = 1)
```

## Arguments

- v:

  single-column matrix of numeric values. Will be coerced to a vector.

- w:

  Numeric vector of individual weights with same length as \`v\`.

- disp:

  Numeric scalar target value to compare the weighted mean against.
  Defaults to 0.

- direction:

  numeric scalar. Multiplier applied to the final metric value to
  control orientation. Use 1 for default direction, -1 to invert the
  sign.

## Value

A single numeric value representing the weighted sum of squared
deviations from \`disp\`.

## Details

The function computes: \$\$ \sum_i w_i (v_i - disp)^2 \$\$

Direction is applied after metric computation and does not alter the
underlying metric definition.

## Examples

``` r
v <- as.matrix(c(1, 2, 3))
w <- c(0.2, 0.3, 0.5)
sum_of_squared_difference(v, w, disp = 2)
#> [1] 0.7
```

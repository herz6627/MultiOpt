# Compute weighted mean absolute deviation from a target value

Calculates the weighted mean of absolute differences between a numeric
vector (\`v\`) and a displacement value (\`disp\`). Each absolute
deviation is weighted by a corresponding value in \`w\`, and the result
is normalized by the sum of weights.

## Usage

``` r
weighted_mean_of_absolute_difference(v, w, disp = 0, direction = 1)
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

A single numeric value representing the weighted mean absolute deviation
of \`v\` from \`disp\`.

## Details

This can be used as a loss function to quantify average absolute
departure from a reference value, for example deviation of sampled
environmental values from a target site condition.

The function computes: \$\$ \frac{\sum_i w_i \|v_i - disp\|}{\sum_i w_i}
\$\$

Direction is applied after metric computation and does not alter the
underlying metric definition.

## Examples

``` r
v <- as.matrix(c(1, 2, 3))
w <- c(0.2, 0.3, 0.5)
weighted_mean_of_absolute_difference(v, w, disp = 2)
#> [1] 0.7
```

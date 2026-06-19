# Compute weighted deviation of a vector from a target value

Calculates the absolute difference between a weighted mean of a numeric
vector (\`vs\`) and a specified target value (\`disp\`). Weights are
applied element-wise and normalized by their sum. This could be used to
minimize the use of individuals low levels of heterozygosity, or to
place arbitrary constraints on genotype composition.

## Usage

``` r
weighted_mean_of_vector(v, w, disp = 0, direction = 1)
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

A single numeric value representing the absolute deviation between the
weighted mean of \`v\` and \`disp\`.

## Details

The function computes: \$\$ \| ( \sum v_i w_i / \sum w_i ) - disp \|
\$\$

Direction is applied after metric computation and does not alter the
underlying metric definition.

## Examples

``` r
v <- as.matrix(c(1, 2, 3))
w <- c(0.2, 0.3, 0.5)
weighted_mean_of_vector(v, w, disp = 2)
#> [1] 0.3
```

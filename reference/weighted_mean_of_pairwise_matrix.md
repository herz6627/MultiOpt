# Compute weighted mean pairwise similarity (or distance)

Calculates the mean value of a pairwise similarity (or distance) matrix,
applying weights to individuals. The function supports two cases:

## Usage

``` r
weighted_mean_of_pairwise_matrix(v, w, direction = 1)
```

## Arguments

- v:

  Square pairwise matrix of similarities or distances between
  individuals. Must have dimensions n × n.

- w:

  Numeric vector of weights (length must match nrow(v)).

- direction:

  Numeric scalar. Multiplier applied to the final metric value to
  control orientation. Use 1 for default direction, -1 to invert the
  sign.

## Value

A single numeric value representing the weighted mean pairwise value.
Returns NA if fewer than two weighted individuals remain.

## Details

For binary weights (w == 1 or 0), the function reduces the matrix to
selected individuals and computes (to save time): \$\$ mean(sm\_{ij}) \\
for \\ i \< j \$\$

For general weights, it computes: \$\$ \frac{\sum\_{i \le j} w_i w_j
sm\_{ij}}{\sum_i w_i (\sum_i w_i - 1)/2} \$\$ with a correction term for
diagonal contributions due to repeated sampling.

Direction is applied after metric computation and does not alter the
underlying metric definition.

## Examples

``` r
v <- matrix(runif(16), 4, 4)
diag(v) <- 0
w <- c(1, 2, 0, 1)
weighted_mean_of_pairwise_matrix(v, w)
#> [1] 0.4025603
```

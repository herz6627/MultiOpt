# Calculate the Wasserstein Distance Between Two Empirical Distributions

Calculates the Wasserstein distance between two one-dimensional
empirical distributions. The Wasserstein distance quantifies the amount
of movement required to transform one distribution into another. Smaller
values indicate more similar distributions.

## Usage

``` r
wasserstein(
  observations_a,
  observations_b,
  p = 1,
  weights_a = NULL,
  weights_b = NULL
)
```

## Arguments

- observations_a:

  Numeric vector containing observations from the first distribution
  (e.g., source population trait values).

- observations_b:

  Numeric vector containing observations from the second distribution
  (e.g., selected collection trait values).

- p:

  Numeric value specifying the Wasserstein order. The default is
  `p = 1`, corresponding to the first Wasserstein distance. Higher
  values increasingly weight larger distributional differences.

- weights_a:

  Optional weights for observations in observations_a.

- weights_b:

  Optional weights for observations in observations_b.

## Value

Numeric value representing the p-Wasserstein distance. Smaller values
indicate distributions are more similar.

## Details

For one-dimensional distributions, the p-Wasserstein distance is
calculated as:

\$\$ W_p(F,G) = \left(\int_0^1 \|Q_F(u)-Q_G(u)\|^p du\right)^{1/p} \$\$

where \\Q_F\\ and \\Q_G\\ are the quantile functions of distributions
\\F\\ and \\G\\, respectively.

When `p = 1`, the metric represents the average distance that trait
values must be moved for one distribution to match the other. This makes
it particularly interpretable for evaluating phenotypic representation
in conservation collections.

## References

D. Schuhmacher, B. Bähre, N. Bonneel, C. Gottschlich, V. Hartmann, F.
Heinemann, B. Schmitzer and J. Schrieber (2024). transport: Computation
of Optimal Transport Plans and Wasserstein Distances. R package version
0.15-0. https://cran.r-project.org/package=transport

Vallender, S. S. (1974). Calculation of the Wasserstein Distance Between
Probability Distributions on the Line. Theory of Probability & Its
Applications. https://doi.org/10.1137/1118101

## Examples

``` r
# Compare a population trait distribution to a selected collection

population_trait <- rnorm(1000)
collection_trait <- sample(population_trait, 50)

wasserstein(
  population_trait,
  collection_trait
)
#> [1] 0.1117715
```

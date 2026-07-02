# Measures

``` r

library(MultiOpt)
```

## Overview

The optimization routines in **MultiOpt** are designed to maximize or
minimize one or more objectives of a collection. Each objective is
represented by a **measure function**, which evaluates a candidate
solution using a vector of individual selection weights.

Every measure function follows the same basic interface:

- **`v`** — a matrix containing data describing all candidate
  individuals.
- **`w`** — a vector of optimization weights, one for each individual.
- Additional optional arguments depending on the measure.

Current measures allow optimization of

- genetic diversity,
- adaptive allele frequencies,
- continuous trait distributions,
- environmental matching,
- pairwise similarity or distance, and
- arbitrary user-defined objectives.

------------------------------------------------------------------------

## Using `calculate_measure()`

Rather than calling each measure independently during optimization,
**MultiOpt** evaluates multiple objectives simultaneously using
[`calculate_measure()`](https://herz6627.github.io/MultiOpt/reference/calculate_measure.md).

``` r

trait_data <- list(
  genetics = geno_matrix,
  climate  = climate_vector,
  kinship  = kinship_matrix
)

measure_list <- list(
  genetics = nei_diversity,
  climate  = sum_of_squared_difference,
  kinship  = weighted_mean_of_pairwise_matrix
)

measure_args <- list(
  climate = list(disp = 15),
  kinship = list(direction = -1)
)

calculate_measure(
  list_of_trait_data = trait_data,
  list_of_measures = measure_list,
  list_of_args = measure_args,
  w = weights
)
```

## Summary of Available Measures

| Function | Required Data | Returns | Typical Objective |
|:---|:---|:---|:---|
| [`nei_diversity()`](https://herz6627.github.io/MultiOpt/reference/nei_diversity.md) | Genotype matrix | Mean expected heterozygosity | Maximize diversity |
| [`shannon_diversity()`](https://herz6627.github.io/MultiOpt/reference/shannon_diversity.md) | Genotype matrix | Hill diversity (q = 0, 1, 2) | Maximize diversity |
| [`allele_enrichment()`](https://herz6627.github.io/MultiOpt/reference/allele_enrichment.md) | Genotype matrix | Mean desirable allele frequency | Maximize adaptive alleles |
| [`weighted_mean_of_vector()`](https://herz6627.github.io/MultiOpt/reference/weighted_mean_of_vector.md) | Single continuous trait | Distance from target mean | Match target average |
| [`sum_of_squared_difference()`](https://herz6627.github.io/MultiOpt/reference/sum_of_squared_difference.md) | Single continuous trait | Weighted squared deviation | Penalize large deviations |
| [`weighted_mean_of_absolute_difference()`](https://herz6627.github.io/MultiOpt/reference/weighted_mean_of_absolute_difference.md) | Single continuous trait | Mean absolute deviation | Minimize average deviation |
| [`weighted_mean_of_pairwise_matrix()`](https://herz6627.github.io/MultiOpt/reference/weighted_mean_of_pairwise_matrix.md) | Pairwise similarity or distance matrix | Mean pairwise value | Maximize or minimize similarity |

------------------------------------------------------------------------

## Genetic Diversity Measures

### Nei Diversity

#### Description

[`nei_diversity()`](https://herz6627.github.io/MultiOpt/reference/nei_diversity.md)
estimates expected heterozygosity across loci.

#### Required data

- genotype matrix
- rows represent individuals
- columns represent loci
- genotypes coded as 0, 1, or 2 copies of the focal allele

#### Interpretation

Higher values indicate greater expected heterozygosity and therefore
greater representation of genetic diversity.

------------------------------------------------------------------------

### Shannon Diversity

#### Description

[`shannon_diversity()`](https://herz6627.github.io/MultiOpt/reference/shannon_diversity.md)
calculates Hill diversity numbers from allele frequencies.

Different values of `q` emphasize different components of diversity.

| q   | Metric            | Interpretation                   |
|:----|:------------------|:---------------------------------|
| 0   | Allelic richness  | Counts polymorphic loci equally  |
| 1   | Shannon diversity | Balances common and rare alleles |
| 2   | Simpson diversity | Emphasizes common alleles        |

#### Required data

The same genotype matrix used for
[`nei_diversity()`](https://herz6627.github.io/MultiOpt/reference/nei_diversity.md).

------------------------------------------------------------------------

### Allele Enrichment

#### Description

[`allele_enrichment()`](https://herz6627.github.io/MultiOpt/reference/allele_enrichment.md)
estimates the weighted frequency of desirable alleles rather than
overall diversity.

Unlike diversity measures, this function explicitly favors specific
alleles that are assumed to contribute positively to fitness or
management objectives.

The function assumes

    2 > 1 > 0

where genotype dosage reflects increasing benefit.

Optional locus weights allow some loci to contribute more strongly than
others.

#### Required data

- genotype matrix
- optional locus importance weights

#### Optional arguments

- `loc` — locus weights
- `rec = TRUE` — treat heterozygotes as non-carriers

------------------------------------------------------------------------

## Continuous Trait Measures

These functions operate on a **single continuous variable**, supplied as
a one-column matrix.

Examples include

- flowering date
- elevation
- annual precipitation
- heterozygosity
- seed mass
- climate of origin

------------------------------------------------------------------------

### Weighted Mean of a Vector

#### Description

[`weighted_mean_of_vector()`](https://herz6627.github.io/MultiOpt/reference/weighted_mean_of_vector.md)
measures how closely the weighted average of selected individuals
matches a target value.

The objective is

``` math
\left|
\frac{\sum w_i x_i}{\sum w_i}
-
target
\right|
```

#### Required data

- one-column matrix
- target value (`disp`)

#### Notes

Only the weighted average is considered, so large positive and negative
deviations may offset one another.

------------------------------------------------------------------------

### Sum of Squared Difference

#### Description

[`sum_of_squared_difference()`](https://herz6627.github.io/MultiOpt/reference/sum_of_squared_difference.md)
calculates the weighted sum of squared deviations from a target value.

Large deviations receive substantially larger penalties.

#### Required data

- one-column matrix
- target value (`disp`)

------------------------------------------------------------------------

### Weighted Mean Absolute Difference

#### Description

[`weighted_mean_of_absolute_difference()`](https://herz6627.github.io/MultiOpt/reference/weighted_mean_of_absolute_difference.md)
computes the weighted average absolute deviation from a target.

Unlike squared deviation, all deviations contribute linearly.

#### Required data

- one-column matrix
- target value (`disp`)

#### Advantages

Less sensitive to outliers than squared deviation.

------------------------------------------------------------------------

## Pairwise Matrix Measures

### Weighted Mean of Pairwise Matrix

#### Description

[`weighted_mean_of_pairwise_matrix()`](https://herz6627.github.io/MultiOpt/reference/weighted_mean_of_pairwise_matrix.md)
evaluates a candidate collection using an existing pairwise similarity
or distance matrix.

The matrix may describe

- genetic distance
- kinship
- geographic distance
- ecological similarity
- phylogenetic distance
- functional trait distance

or any other pairwise relationship.

#### Required data

A square symmetric matrix with rows and columns corresponding to
individuals.

#### Typical applications

- maximizing functional diversity
- minimizing relatedness
- maximizing geographic coverage
- ecological representation

------------------------------------------------------------------------

## Optimization Direction

All measure functions include a `direction` argument.

| Value | Optimization        |
|:------|:--------------------|
| `1`   | maximize the metric |
| `-1`  | minimize the metric |

The underlying metric itself is unchanged; only its optimization
direction is reversed. This avoids the need for separate maximizing and
minimizing implementations of the same statistic.

------------------------------------------------------------------------

## Combining Multiple Objectives

The primary strength of **MultiOpt** is the ability to optimize several
objectives simultaneously.

For example, a seed collection may seek to

- maximize neutral genetic diversity,
- maximize adaptive allele frequencies,
- minimize climatic mismatch to a restoration site, and
- maximize geographic representation.

Because every measure shares the same interface, additional measures can
easily be incorporated by writing a function that accepts a matrix `v`,
a vector of weights `w`, and returns a single numeric value.

This modular design makes the optimization framework flexible enough to
accommodate a wide variety of conservation, restoration, and breeding
applications.

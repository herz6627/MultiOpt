# Nei diversity

Calculates Nei genetic diversity from genotype data.

## Usage

``` r
nei_diversity(v, w = NULL, direction = 1)
```

## Arguments

- v:

  Genotype matrix (rows = individuals, columns = loci)

- w:

  Optional vector of weights

- direction:

  Numeric scalar. Must be 1 or -1. Applied as a multiplicative factor to
  the computed metric to control its optimisation direction.

## Value

Numeric scalar of calculated Nei diversity.

## Details

Direction is applied after metric computation and does not alter the
underlying metric definition.

## Note

Nei, M. (1973). Analysis of gene diversity in subdivided populations.
Proceedings of the National Academy of Sciences of the United States of
America, 70, 3321–3323.

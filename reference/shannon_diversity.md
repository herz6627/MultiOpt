# Shannon diversity

Shannon diversity

## Usage

``` r
shannon_diversity(v, w = NULL, q = 1, direction = 1)
```

## Arguments

- v:

  Numeric genotype matrix with individuals in rows and loci in columns.
  Entries are assumed to be allele counts (0, 1, 2).

- w:

  Optional numeric vector of individual weights (length must equal
  nrow(v)). If NULL, all individuals are treated equally.

- q:

  Diversity order:

  0

  :   Locus polymorphism indicator (monomorphic = 1, polymorphic = 2)

  1

  :   Shannon diversity (Hill number q = 1)

  2

  :   Simpson diversity (Hill number q = 2)

  @param direction numeric scalar. Multiplier applied to the final
  metric value to control orientation. Use 1 for default direction, -1
  to invert the sign.

## Value

A single numeric value representing the mean per-locus Hill diversity
across all loci.

## Details

Direction is applied after metric computation and does not alter the
underlying metric definition.

## Note

I have kept the function true to the version found in OptGenMix, but
have clarified what 'q' is actually doing. The original description in
OptGenMix was sparse, so I have had to do some extrapolation.

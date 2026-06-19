# Desirable allele enrichment

Calculates enrichment of desirable allele. Given genotypes (v), a vector
(w) of weight values (equal length to the number of individuals, a
vector (v) of values describing the importance of each locus, estimates
an index of enrichment of preferred alleles. Assumes genotype fitness 2
\> 1 \> 0.

## Usage

``` r
allele_enrichment(v, w = NULL, loc = NULL, rec = FALSE, direction = 1)
```

## Arguments

- v:

  Genotype matrix (individuals × loci), coded as allele dosage (0, 1,
  2).

- w:

  Optional numeric vector of individual weights (length = nrow(v)). If
  NULL, all individuals are weighted equally.

- loc:

  Optional numeric vector of locus weights (length = ncol(v)). If NULL,
  loci are weighted equally.

- rec:

  Logical. If TRUE, heterozygotes (1) are treated as 0, enforcing a
  recessive model where only homozygotes for the allele contribute.
  @param direction numeric scalar. Multiplier applied to the final
  metric value to control orientation. Use 1 for default direction, -1
  to invert the sign.

## Value

A single numeric value representing the weighted allele enrichment
index. Higher values indicate greater enrichment of the allele across
individuals and loci.

## Details

Direction is applied after metric computation and does not alter the
underlying metric definition.

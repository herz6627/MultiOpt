# Nei's Genetic Diversity

Calculates Nei's genetic diversity (expected heterozygosity) across loci
in a population. Allele frequencies are calculated from a genotype
matrix, with the option to incorporate individual-level weights.

## Usage

``` r
nei_diversity(v, w = NULL, opt = c("speed", "storage"), chunk_size = 10000)
```

## Arguments

- v:

  A genotype matrix with individuals in rows and loci in columns.
  Genotypes should be coded as the number of copies of an allele (e.g.,
  0, 1, or 2).

- w:

  An optional numeric vector of individual weights. If supplied,
  weighted allele frequencies are calculated.

- opt:

  A character string specifying the calculation method. "speed" uses a
  vectorized calculation, while "storage" calculates diversity
  iteratively. Defaults to "speed".

- chunk_size:

  Numeric value for how large (number of SNPs) to chunk large genotype
  datasets by. Only used if ncol(v) \> 10,000.

## Value

A numeric value representing the mean Nei's genetic diversity (expected
heterozygosity) across loci.

## Details

For each locus, the frequency of the focal allele is calculated as the
proportion of allele copies in the population. Nei's gene diversity is
then calculated as the expected heterozygosity.

For each locus, the allele frequency is calculated as:

\$\$p_i = \frac{\sum_j g\_{ji}}{2n}\$\$

where \\g\_{ji}\\ is the genotype of individual \\j\\ at locus \\i\\,
and \\n\\ is the number of individuals.

When weights are provided, the weighted allele frequency is:

\$\$p_i = \frac{\sum_j w_j g\_{ji}}{2\sum_j w_j}\$\$

Nei's gene diversity for each locus is calculated as:

\$\$H_i = 1 - p_i^2 - (1-p_i)^2 = 2p_i(1-p_i)\$\$

The final value is the mean diversity across all \\L\\ loci:

\$\$H = \frac{1}{L}\sum\_{i=1}^{L} H_i\$\$

Assumes loci are bi-allelic.

The "speed" and "storage" options calculate the same quantity using
different computational approaches.

## References

Nei, M. (1973). Analysis of gene diversity in subdivided populations.
Proceedings of the National Academy of Sciences of the United States of
America, 70, 3321–3323.

Based on code by Jason Bragg in the OptGenMix package
(jasongbragg/OptGenMix) as used in these publications:

Bragg, J. G., Cuneo, P., Sherieff, A., & Rossetto, M. (2020). Optimizing
the genetic composition of a translocation population: Incorporating
constraints and conflicting objectives. Molecular Ecology Resources,
20(1), 54–65. https://doi.org/10.1111/1755-0998.13074

Bragg, J. G., Yap, J.-Y. S., Wilson, T., Lee, E., & Rossetto, M. (2021).
Conserving the genetic diversity of condemned populations: Optimizing
collections and translocation. Evolutionary Applications, 14(5),
1225–1238. https://doi.org/10.1111/eva.13192

## Examples

``` r
# Create a small example genotype matrix.
# Rows represent individuals and columns represent loci.
# Genotypes are coded as the number of copies of the focal allele (0, 1, or 2).
geno <- matrix(
c(0, 1, 2,
1, 1, 0,
2, 0, 1,
1, 2, 1),
nrow = 4,
byrow = TRUE
)

# Calculate Nei's gene diversity.
nei_diversity(geno)
#> [1] 0.5

# Use the storage-oriented calculation.
nei_diversity(geno, opt = "storage")
#> [1] 0.5

# Calculate weighted Nei's gene diversity.
weights <- c(1, 0.5, 1.5, 1)
nei_diversity(geno, w = weights)
#> [1] 0.4921875
```

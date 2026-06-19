# Combine two archives and reduce if needed

Combines two archives from \`multiopt_sa\` and reduces to one archive,
filtering to only non-dominated values. Both archives should be produced
from the same simulated annealing arguments and trait data sets. Order
of which archive value is assigned to archive1 or archive2 is arbitrary
and won't make a difference.

## Usage

``` r
combine_archives(archive1, archive2)
```

## Arguments

- archive1:

  Archive from \`multiopt_sa\` (out1\$archive)

- archive2:

  Archive from \`multiopt_sa\` (out2\$archive)

## Value

a single archive list matching archive output style as multiopt_sa. A
list containing:

- archive_summary:

  Values for each non-dominated solution with a column for each trait.

- archive_weights:

  Weights for each non-dominated solution (row) with a column for each
  individual.

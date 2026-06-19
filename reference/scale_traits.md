# Min-max scale trait datasets to the range \[0, 1\]

Applies min-max scaling to each element in a list of trait datasets,
transforming values to the range \[0, 1\]. Scaling is performed
independently for each dataset using.

## Usage

``` r
scale_traits(trait_list)
```

## Arguments

- trait_list:

  A list of matrices (with a single column for a trait vector or a
  pairwise matrix), or other numeric objects to be scaled independently.

## Value

A list of scaled trait datasets with the same structure and names as
\`trait_list\`.

## Details

Each element of \`trait_list\` is processed independently using min-max
normalization:

\$\$ x\_{scaled} = \frac{x - \min(x)}{\max(x) - \min(x)} \$\$

The output preserves the original list structure.

If a dataset contains missing values, they are ignored during scaling
but retained in the output.

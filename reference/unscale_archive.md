# Back transform output to original scale

This function uses the original transformed trait data provided in
\`multiopt_sa\` or \`rand_multiopt\` and the corresponding output from
\`multiopt_sa\` or \`rand_multiopt\` to back transform archive values to
the original trait scale. This function assumes trait data was
transformed using \`scale_traits\`.

## Usage

``` r
unscale_archive(trait_list, archive_output)
```

## Arguments

- trait_list:

  Original list of trait data as used in \`multiopt_sa\`.

- archive_output:

  Archive output from \`multiopt_sa\` (out\$archive), \`rand_multiopt\`
  (out\$archive) or \`explore_pareto\`.

## Value

List with archive values and weights matching formatting of
\`archive_output\`.

## Details

Note that this function takes the absolute value of transformed trait
variables.

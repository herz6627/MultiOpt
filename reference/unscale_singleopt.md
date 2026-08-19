# Back transform output to original scale

This function uses the original transformed trait data provided to
\`singleopt_context\` and the corresponding output from
\`singleopt_context\` to back transform archive values to the original
trait scale. This function assumes trait data was transformed using
\`scale_traits\`.

## Usage

``` r
unscale_singleopt(trait_list, singleopt_output, skip_traits = NULL)
```

## Arguments

- trait_list:

  Original list of trait data as used in \`singleopt_context\`.

- singleopt_output:

  Archive output from \`singleopt_context\`.

- skip_traits:

  Optional characted vector specifing if any traits should be excluded
  from back transformation.

## Value

List with archive values and weights matching formatting of
\`singleopt_output\`.

## Details

Note that this function takes the absolute value of transformed trait
variables.

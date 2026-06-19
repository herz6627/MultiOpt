# Back transform output to original scale

This function uses the original transformed trait data provided in
\`multiopt_sa\` and the corresponding output from \`multiopt_sa\` to
back transform simulation output to the original trait scale. This
function assumes trait data was transformed using \`scale_traits\`.

## Usage

``` r
unscale_multiopt(trait_list, multiopt_output)
```

## Arguments

- trait_list:

  Original list of trait data as used in \`multiopt_sa\`.

- multiopt_output:

  Unmodified output from \`multiopt_sa\`.

## Value

Object of same dimensions and formatting as \`multiopt_output\`.

## Details

Note that this function takes the absolute value of transformed trait
variables.

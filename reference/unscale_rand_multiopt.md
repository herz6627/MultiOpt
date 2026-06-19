# Back transform output to original scale

This function uses the original transformed trait data provided in
\`rand_multiopt\` and the corresponding output from \`rand_multiopt\` to
back transform simulation output to the original trait scale. This
function assumes trait data was transformed using \`scale_traits\`.

## Usage

``` r
unscale_rand_multiopt(trait_list, rand_multiopt_output)
```

## Arguments

- trait_list:

  Original list of trait data as used in \`rand_multiopt\`.

- rand_multiopt_output:

  Unmodified output from \`rand_multiopt\`.

## Value

Object of same dimensions and formatting as \`rand_multiopt_output\`.

## Details

Note that this function takes the absolute value of transformed trait
variables.

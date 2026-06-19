# Reverse min-max scaling using a reference dataset

Transforms values from a min-max scaled scale back to their original
scale using the minimum and maximum of a reference unscaled dataset.

## Usage

``` r
min_max_unscale(x_scaled, x_unscaled)
```

## Arguments

- x_scaled:

  Numeric vector, matrix, or array containing values scaled to the range
  \[0, 1\].

- x_unscaled:

  Numeric vector, matrix, or array providing the reference scale used
  for back-transformation.

## Value

An object with the same structure as \`x_scaled\`, transformed back to
the scale of \`x_unscaled\`.

## Details

The transformation is:

\$\$ x = x\_{scaled} \times (\max(x\_{ref}) - \min(x\_{ref})) +
\min(x\_{ref}) \$\$

The minimum and maximum used for back-transformation are calculated from
\`x_unscaled\` with missing values ignored.

This function assumes that \`x_scaled\` was originally generated using
min-max scaling based on the same reference distribution.

# Min-max scale

Transforms numeric values using min-max normalization, rescaling values
to the interval \[0, 1\] according to:

## Usage

``` r
min_max_scale(x)
```

## Arguments

- x:

  A numeric vector, matrix, or array.

## Value

An object with the same structure as \`x\`, with values scaled to the
range \[0, 1\].

## Details

\$\$ x\_{scaled} = \frac{x - \min(x)}{\max(x) - \min(x)} \$\$

Missing values are ignored when calculating minima and maxima. If
missing values are present, a message is produced.

Scaling is performed using the global minimum and maximum of \`x\`.
Missing values (\`NA\`) are retained in the output.

If all non-missing values are identical, the function will return
\`NaN\` values due to division by zero.

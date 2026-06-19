# Visualize selection frequency across individuals and traits

Produces plots showing how frequently individuals are selected across
replicates. Supports both single-trait and multi-trait visualization
modes.

## Usage

``` r
plot_selection(individs_selected, trait_list)
```

## Arguments

- individs_selected:

  Matrix or data structure indicating selected individuals across
  replicates as produced by \`rand_multiopt\` (out\$individs_selected).

- trait_list:

  Named list of trait data frames, each containing a single column of
  trait values. Must not contain multi-column objects (e.g. pairwise
  matrices are not supported). This data does not need to be the same
  supplied to \`rand_multiopt\`, but does need to be in the same order
  for each individual.

## Value

A ggplot2 object (single trait case) or a patchwork object containing
pairwise trait scatterplots (multi-trait case).

## Details

For each individual, selection frequency is computed as the number of
replicates in which the individual was selected (ignoring weights).

The function operates in two modes:

\*\*Single trait case\*\*

- Aggregates trait values with selection frequency.

- Produces a bar plot ordered by trait value.

- Fill indicates number of replicates in which each individual was
  selected.

\*\*Multi-trait case\*\*

- Constructs a data frame of trait values.

- Computes selection frequency per individual.

- Generates all pairwise trait scatterplots.

- Uses fill color to indicate selection frequency.

Selection frequency is computed as the row-wise count of nonzero entries
in \`individs_selected\`.

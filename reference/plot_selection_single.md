# Visualize trait values of selected individuals

Produces plots showing how frequently individuals are selected across
optimization replicates in relation to their (and other) trait values.
This is similar to \`plot_selection\` except input values are formatted
for \`singleopt_context\` output.

## Usage

``` r
plot_selection_single(individs_selected, trait_list)
```

## Arguments

- individs_selected:

  A named list of matrices indicating selected individuals across
  replicates (e.g. \`out\$individs_selected\` returned by
  \`singleopt_context\`). Each element represents the results from a
  different optimization scenario.

- trait_list:

  Named list of trait data frames, each containing a single column of
  trait values. This data does not need to be the same data supplied to
  \`singleopt_context\`; in many cases it is preferable to use unscaled
  or raw trait values for visualization. Trait values must be in the
  same individual order as \`individs_selected\`. Pairwise matrices or
  multi-column objects are not supported.

## Value

A named list of plots, with one plot per optimization scenario. For a
single trait, each element is a ggplot2 bar plot. For multiple traits,
each element is a patchwork object containing all pairwise trait
scatterplots.

## Details

For each individual, selection frequency is computed as the number of
replicates in which the individual was selected (ignoring selection
weights).

The function operates in two modes:

\*\*Single trait case\*\*

- Combines trait values with selection frequency.

- Orders individuals by trait value.

- Produces a bar plot of trait values.

- Bar fill indicates the number of replicates in which each individual
  was selected.

\*\*Multi-trait case\*\*

- Combines all trait values into a single data frame.

- Computes selection frequency for each individual.

- Generates scatterplots for all pairwise combinations of traits.

- Point fill indicates the number of replicates in which each individual
  was selected.

- Pairwise plots are combined into a single figure using \`patchwork\`.

Selection frequency is computed as the number of replicates in which an
individual has a nonzero selection weight.

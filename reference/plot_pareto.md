# Visualize pairwise Pareto trade-offs

Generates a grid of pairwise scatterplots showing trade-offs among
traits in a Pareto archive. Optionally overlays results from
multi-objective and single-objective optimization runs for comparison.

## Usage

``` r
plot_pareto(
  archive_list,
  multi_list = NULL,
  single_list = NULL,
  measure_type = NULL
)
```

## Arguments

- archive_list:

  A list containing Pareto archive output, expected to include an
  element \`archive_summary\` with trait values. Either \$archive from
  \`multiopt_sa\` or \`rand_multiopt\` or output from \`explore_pareto\`

- multi_list:

  Optional list containing multi-objective optimization results from
  \`rand_multiopt\` or \`multiopt_sa\`. Must include either
  \`measure_summaries\` or \`final_selection\$measure_summary\`.

- single_list:

  Optional list of single-objective optimization results from
  \`singleopt_context\`.

- measure_type:

  Optional list of direction of optimization. List elements should have
  the same names as traits in archive_list and values should be
  \`"minimize"\`, \`"maximize"\`, or \`"diversify"\`.

## Value

A patchwork object consisting of multiple ggplot2 pairwise scatterplots
with shared legends.

## Details

Required packages: ggplot2, dplyr, patchwork.

Note that simulated annealing as used here is trying to maximize values.
Therefore, raw output from MultiOpt SA simulations will be maximizing.
If you back transform trait values to their original scale then you will
be possibly maximizing, minimizing, or diversifying.

## Examples

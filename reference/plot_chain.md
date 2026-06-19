# Plot trait value chains over iterations

Visualizes the temporal trajectory of trait values through the simulated
annealing procedure.

## Usage

``` r
plot_chain(chain_list)
```

## Arguments

- chain_list:

  Chain output from \`multiopt_sa\` (out\$chain).

## Value

A ggplot object showing trait trajectories over iteration steps, faceted
by trait.

## Details

This is primarily intended for diagnosing convergence, mixing, or
temporal dynamics in optimization or sampling chains.

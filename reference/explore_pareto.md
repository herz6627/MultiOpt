# Explore Pareto Front

Expands upon the archive values from an initial simulated annealing run
to expand the known coverage of the Pareto Front.

## Usage

``` r
explore_pareto(
  multiopt_archive_output,
  nd_samples = 100,
  max_steps = 5000,
  verbose = T,
  ...
)
```

## Arguments

- multiopt_archive_output:

  Unmodified archive output from \`multiopt_sa\`. out\$archive

- nd_samples:

  How many observations of the Pareto Front to explore. Will only impact
  settings for multiopt_sa internally, results may be larger than the
  provided \`nd_samples\` Should probably be set higher than what you
  initially set for \`multiopt_sa\` for maximum Pareto front
  exploration.

- max_steps:

  Maximum number of simulated annealing iterations. Since we are simply
  exploring the Pareto Front, not aiming for an optimum solution, this
  value does not need to be as large as when used for the full simulated
  annealing run.

- ...:

  The same arguments used to initially run \`multiopt_sa\`. If you
  supply \`initial_weights\` it will be overwritten.

## Value

Returns a list matching the formatting of \`multiopt_archive_output\`,
but with hopefully more observations of the front.

## Details

Can use \`unscale_archive\` to back transform data to original trait
scales.

## References

Bragg, J. G., van der Merwe, M., Yap, J.-Y. S., Borevitz, J., &
Rossetto, M. (2022). Plant collections for conservation and restoration:
Can they be adapted and adaptable? Molecular Ecology Resources, 22(6),
2171–2182. https://doi.org/10.1111/1755-0998.13605

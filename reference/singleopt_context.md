# Run single-objective simulated annealing in context of additional traits

This function aims to address "If I prioritize a single trait for
simulated annealing, what are the results of the other traits?." This is
complimentary analysis to \`multiopt_sa\` and \`rand_multiopt\`,
allowing users to explore what trade offs happen when incorporating
multiple objectives. No chains or archives are returned.

## Usage

``` r
singleopt_context(
  trait_list,
  measure_list,
  measure_args_list,
  n_t = NULL,
  weights_min = NULL,
  weights_max = NULL,
  initial_weights = NULL,
  max_steps = 10000,
  max_t = 1,
  p_depends_delta = F,
  acceptance_multipliers = 1,
  verbose = T,
  n_runs = 1,
  parallel = F
)
```

## Arguments

- trait_list:

  List of trait matrices. Each matrix should contain individuals in rows
  and trait values in columns.

- measure_list:

  List of objective functions corresponding to each element of
  \`trait_list\`. Possible functions include: \`nei_diversity\`,
  \`shannon_diversity\`, \`allele_enrichment\`,
  \`weighted_mean_of_vector\`, \`sum_of_squared_difference\`,
  \`weighted_mean_of_absolute_difference\`,
  \`weighted_mean_of_pairwise_matrix\`

- measure_args_list:

  List containing additional arguments supplied to each measure
  function.

- n_t:

  Integer giving the total number or total weight of individuals to
  select.

- weights_min:

  Optional minimum allowable value(s) for weights.

- weights_max:

  Optional maximum allowable value(s) for weights. Either a single
  numeric value or a vector with length equal to the number of
  individuals.

- initial_weights:

  Optional vector of starting weights for the optimization.

- max_steps:

  Maximum number of simulated annealing iterations.

- max_t:

  Initial annealing temperature.

- p_depends_delta:

  Logical; if \`TRUE\`, acceptance probabilities depend on the magnitude
  of objective change.

- acceptance_multipliers:

  Numeric scalar or vector of the same length as \`trait_list\`
  controlling acceptance probability for each objective.

- verbose:

  If TRUE, prints status updates in the console.

- n_runs:

  Integer. How many randomizations to run.

- parallel:

  Logical; whether to run in parallel. If TRUE, user will need to set up
  the session using \`future::plan()\`

## Value

A list of the same length of \`trait_list\` with an element for each
trait found in \`trait_list\`. Within each list element, a table with
measure results for each trait is returned, with the column with the
same name as the list element being the single objective used in the
simulated annealing run(s).

## Details

Compatible with parallel computing when n_reps \> 1. See
\`rand_multiopt()\` for more details.

Output value can be back transformed using \`unscale_singleopt\`.

## Examples

``` r

set.seed(12345)
n = 100
x = rnorm(n = n, mean = 120, sd = 2)
y = x * 3 + rnorm(n = n, mean = 0, sd = 20)
dat = data.frame(x = x, y = y)

trait_list = list(
x = as.matrix(dat$x),
y = as.matrix(dat$y)
)

measure_list = list(
 x = weighted_mean_of_vector,
 y = weighted_mean_of_vector
)

args_list = list(
 x = NULL,
 y = list(direction = -1) # to minimize the mean of this trait
)

trait_list_scaled <- scale_traits(trait_list)

set.seed(123)
singleopt_context(
 trait_list = trait_list_scaled,
 measure_list = measure_list,
 measure_args_list = args_list,
 n_t = 10,
 verbose = F
)
#> $x
#>          x          y
#> 1 1.061243 -0.8482948
#> 
#> $y
#>   x y
#> 1 0 0
#> 

```

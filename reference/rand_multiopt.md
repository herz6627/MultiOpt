# Multiple runs of simulated annealing

This is a wrapper function around \`multiopt_sa()\` for multiple
replicates of the simulated annealing process. Results are summarized
and condensed. Process can be run in parallel using the \`future\`
package.

## Usage

``` r
rand_multiopt(n_runs = 20, parallel = F, multiopt_args, verbose = T)
```

## Arguments

- n_runs:

  Integer. How many randomizations to run.

- parallel:

  Logical; whether to run in parallel. If TRUE, user will need to set up
  the session using \`future::plan()\`

- multiopt_args:

  List of arguments to supply to \`multiopt_sa\`

- verbose:

  Whether status updates should be printed to the console.

## Value

a list containing:

- measure_summaries:

  Table containing the final objective values.

- individs_selected:

  Matrix of weights for each individual (columns) for each randomization
  (row).

- archive:

  If \`nda = TRUE\`, a list containing archived non-dominated solutions
  (\`archive_summary\`) and associated weights (\`archive_weights\`);
  otherwise \`NULL\`. Archive is the combined archives over all
  randomizations.

## Examples

``` r

library(future)

future::plan(multisession, workers = 4) # multisession indicates this is run on this computer


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
 y = list(direction = -1)
)

trait_list_scaled <- scale_traits(trait_list)

sa_args = list(
trait_list = trait_list_scaled,
measure_list = measure_list,
measure_args_list = args_list,
n_t = 10,
weights_max = 1,
max_steps = 50000,
nda = T,
nd_samples = 500
)

test = rand_multiopt(n_runs = 5, multiopt_args = sa_args, parallel = T)
#> Starting simulation at 2026-06-19 15:52:35.870994 
#> Work completed in 0.72 minutes
str(test)
#> List of 3
#>  $ measure_summaries: num [1:5, 1:2] 0.755 0.794 0.764 0.744 0.594 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:2] "x" "y"
#>  $ individs_selected: num [1:5, 1:100] 0 0 0 0 0 1 0 1 0 0 ...
#>  $ archive          :List of 2
#>   ..$ archive_summary: num [1:77, 1:2] 0.866 0.86 0.849 0.842 0.78 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : NULL
#>   .. .. ..$ : chr [1:2] "x" "y"
#>   ..$ archive_weights: num [1:77, 1:100] 0 0 0 0 0 0 0 0 0 0 ...
#>   .. ..- attr(*, "dimnames")=List of 2
#>   .. .. ..$ : chr [1:77] "" "" "" "" ...
#>   .. .. ..$ : NULL
```

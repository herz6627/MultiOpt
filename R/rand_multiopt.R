#' Multiple runs of simulated annealing
#'
#' This is a wrapper function around `multiopt_sa()` for multiple replicates
#' of the simulated annealing process. Results are summarized and condensed.
#' Process can be run in parallel using the `future` package.
#'
#' @param n_runs Integer. How many randomizations to run.
#' @param parallel Logical; whether to run in parallel. If TRUE, user will need to set up the session using `future::plan()`
#' @param multiopt_args List of arguments to supply to `multiopt_sa`
#' @param verbose Whether status updates should be printed to the console.
#'
#' @returns a list containing:
#' \describe{
#'   \item{measure_summaries}{Table containing the final objective values.}
#'   \item{individs_selected}{Matrix of weights for each individual (columns) for each randomization (row).}
#'   \item{archive}{If `nda = TRUE`, a list containing archived
#'   non-dominated solutions (`archive_summary`) and associated weights (`archive_weights`); otherwise `NULL`.
#'   Archive is the combined archives over all randomizations.}
#' }
#'
#'
#' @export
#'
#' @examples
#'
#' library(future)
#'
#'future::plan(multisession, workers = 4) # multisession indicates this is run on this computer
#'
#'
#' set.seed(12345)
#'n = 100
#'x = rnorm(n = n, mean = 120, sd = 2)
#'y = x * 3 + rnorm(n = n, mean = 0, sd = 20)
#'dat = data.frame(x = x, y = y)
#'
#' trait_list = list(
#'x = as.matrix(dat$x),
#'y = as.matrix(dat$y)
#')
#'
#'measure_list = list(
#'  x = weighted_mean_of_vector,
#'  y = weighted_mean_of_vector
#')
#'
#'args_list = list(
#'  x = NULL,
#'  y = list(direction = -1)
#')
#'
#' trait_list_scaled <- scale_traits(trait_list)
#'
#' sa_args = list(
#'trait_list = trait_list_scaled,
#'measure_list = measure_list,
#'measure_args_list = args_list,
#'n_t = 10,
#'weights_max = 1,
#'max_steps = 50000,
#'nda = T,
#'nd_samples = 500
#')
#'
#'test = rand_multiopt(n_runs = 5, multiopt_args = sa_args, parallel = T)
#'str(test)
#'
rand_multiopt <- function(
    n_runs = 20,
    parallel = F,
    multiopt_args,
    verbose = T
){

  if(!n_runs%%1 == 0 | n_runs < 0) stop("`n_runs` must be positive and an integer.")

  # modify multiopt_sa args
  multiopt_args$save_chain = F # to save space
  multiopt_args$verbose = F

  # run simulation
  start_t = Sys.time()
  if(verbose) message("Starting simulation at ", start_t," \n")

  if(!parallel) {

    out <- replicate(
      n_runs,
      do.call(multiopt_sa, multiopt_args),
      simplify = F)

  } else {

    out <- future.apply::future_replicate(   # using future apply standardizes seeds, allowing for replication
      n_runs,
      do.call(multiopt_sa, multiopt_args),
      simplify = F
    )

  }

  # combine best results
  measure_summary <- do.call(
    rbind,
    lapply(out, function(x) x$final_selection$measure_summary)
  )

  # matrix of weights for each rep
  individs_selected <- do.call(
    rbind,
    lapply(out, function(x) x$final_selection$individs_selected)
  )

  # combine archives if needed
  if (!is.null(multiopt_args$nda) && isTRUE(multiopt_args$nda)) {

    all_archives <- lapply(out, `[[`, "archive")

    archive <- Reduce(combine_archives, all_archives)

    colnames(archive$archive_summary) = names(multiopt_args$trait_list)
  }

  if(verbose) message(
    sprintf(
      "Work completed in %.2f minutes",
      as.numeric(Sys.time() - start_t, units = "mins")
    )
  )

  out = list(
    measure_summaries = measure_summary,
    individs_selected = individs_selected,
    archive = if(!is.null(multiopt_args$nda) && isTRUE(multiopt_args$nda)) archive else NULL
  )

  return(out)
}

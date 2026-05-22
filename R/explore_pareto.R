#' Explore Pareto Front
#'
#' Expands upon the archive values from an initial simulated annealing run
#' to expand the known coverage of the Pareto Front.
#'
#'
#' Can use `unscale_archive` to back transform data to original trait scales.
#'
#' @param multiopt_archive_output Unmodified archive output from `multiopt_sa`. out$archive
#' @param nd_samples How many observations of the Pareto Front to explore.
#' Will only impact settings for multiopt_sa internally, results may be larger than the provided `nd_samples`
#' Should probably be set higher than what you initially set for `multiopt_sa` for maximum Pareto front exploration.
#' @param max_steps Maximum number of simulated annealing iterations.
#' Since we are simply exploring the Pareto Front, not aiming for an optimum solution,
#' this value does not need to be as large as when used for the full simulated annealing run.
#' @param ... The same arguments used to initially run `multiopt_sa`. If you supply `initial_weights` it will be overwritten.
#'
#' @returns Returns a list matching the formatting of `multiopt_archive_output`, but with hopefully more observations of the front.
#'
#' @references   Bragg, J. G., van der Merwe, M., Yap, J.-Y. S., Borevitz, J., & Rossetto, M. (2022). Plant collections for conservation and restoration: Can they be adapted and adaptable? Molecular Ecology Resources, 22(6), 2171–2182. https://doi.org/10.1111/1755-0998.13605
#'
#' @export
explore_pareto <- function(multiopt_archive_output, nd_samples = 100, max_steps = 5000, ...) {

  current_archive <- multiopt_archive_output

  nruns <- nrow(multiopt_archive_output$archive_summary)

  if(nruns == 1) {

    message("No reason to explore the Pareto Front with only 1 trait. Returning original archive.")
    return(multiopt_archive_output)

    }

  for (i in 1:nruns) {

    cat("\rRefining", i, "of", nruns, "archived solutions")

    initial_weights <- multiopt_archive_output$archive_weights[i,] # which individuals to start with

    # re-run SA starting from a previous non-dominated solution
    opt_tmp = multiopt_sa(
      initial_weights = initial_weights,
      nd_samples = nd_samples,
      max_steps = max_steps,
      nda = TRUE,
      verbose = F,
      ...)

    # now need to merge new archive with old
    current_archive <- combine_archives(current_archive, opt_tmp$archive)

    }

  return(current_archive)
}

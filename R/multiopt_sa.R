#' Multi-objective simulated annealing optimization
#'
#' Performs multi-objective (or single-objective) optimization using a simulated annealing
#' framework. Candidate solutions are represented as vectors of weights
#' describing the contribution or selection of individuals, and are
#' iteratively modified to improve one or more objective functions.
#'
#' Objective values are calculated by applying user-supplied measure
#' functions to corresponding trait datasets. Candidate solutions may be
#' accepted even when they worsen one or more objectives, with acceptance
#' probability determined by the current annealing temperature and the
#' specified acceptance parameters.
#'
#' Optionally, a non-dominated archive can be maintained during the search.
#'
#' @param trait_list List of trait matrices. Each matrix should contain
#'   individuals in rows and trait values in columns.
#' @param measure_list List of objective functions corresponding to each
#'   element of `trait_list`. Possible functions include: `nei_diversity`, `shannon_diversity`, `allele_enrichment`, `weighted_mean_of_vector`, `sum_of_squared_difference`, `weighted_mean_of_absolute_difference`, `weighted_mean_of_pairwise_matrix`
#' @param measure_args_list List containing additional arguments supplied
#'   to each measure function.
#' @param n_t Integer giving the total number or total weight of individuals
#'   to select.
#' @param weights_max Optional maximum allowable value(s) for weights.
#'   Either a single numeric value or a vector with length equal to the
#'   number of individuals.
#' @param weights_min Optional minimum allowable value(s) for weights.
#' @param initial_weights Optional vector of starting weights for the
#'   optimization.
#' @param max_steps Maximum number of simulated annealing iterations.
#' @param max_t Initial annealing temperature.
#' @param min_t Minimum allowable temperature when `nda = TRUE`.
#' @param p_depends_delta Logical; if `TRUE`, acceptance probabilities
#'   depend on the magnitude of objective change.
#' @param acceptance_multipliers Numeric scalar or vector of the same length as `trait_list` controlling
#'   acceptance probability for each objective.
#' @param acceptance_multiplier_all_worse Numeric scalar controlling
#'   acceptance probability when all objectives worsen simultaneously.
#' @param nda Logical; if `TRUE`, maintain a non-dominated archive of
#'   Pareto-optimal candidate solutions encountered during optimization.
#' @param nd_samples Maximum number of archived non-dominated solutions
#'   to retain when `nda = TRUE`. Once satisfied, simulated annealing run will stop.
#' @param save_chain Logical; Whether to save the chain of results. Helpful to set to FALSE if trying to save storage.
#' @param verbose If TRUE, prints status updates in the console.
#'
#' @return A list containing:
#' \describe{
#'   \item{final_selection}{List containing the final objective values
#'   (`measure_summary`) and selected weights (`individs_selected`).}
#'   \item{chain}{If `save_chain = TRUE`, a list containing optimization history, including proposed
#'   weights, objective values, and acceptance status for each iteration. Otherwise `NULL`}
#'   \item{archive}{If `nda = TRUE`, a list containing archived
#'   non-dominated solutions (`archive_summary`) and associated weights (`archive_weights`); otherwise `NULL`.}
#' }
#'
#' @export
multiopt_sa <- function(
    trait_list,
    measure_list,
    measure_args_list,
    n_t = NULL,
    weights_max = NULL,
    weights_min = NULL,
    initial_weights = NULL,
    max_steps = 10000,
    max_t = 1,
    min_t = 0,
    p_depends_delta = F,
    acceptance_multipliers = 1,
    acceptance_multiplier_all_worse = 1,
    nda = F,
    nd_samples = 100,
    save_chain = T,
    verbose = T
) {

  objectives = length(trait_list)

  if(verbose) message(paste("\nFound", objectives, "objective(s)."))

  # Run checks --------------------------------------------------------------

  if (!all(vapply(measure_list, is.function, logical(1)))) {
    stop("All supplied measures must be functions")
  }

  if(any(!sapply(trait_list, is.matrix))) stop("All trait data within 'trait_list' must be in matrix format.")

  check_similar_scale(trait_list)


  if (is.null(n_t) && is.null(initial_weights)) stop("Either 'n_t' or 'initial_weights' must be provided.")

  if(is.null(n_t) & !is.null(initial_weights)) n_t = sum(initial_weights)

  if ( !is.null(n_t) & !is.null(initial_weights)) {

    if (n_t != sum(initial_weights)) stop( "Conflict in number of target individuals and initial weights provided \n" )

  }

  if (!(length(trait_list) == length(measure_list) &&
        length(measure_list) == length(measure_args_list))) stop("`trait_list`, `measure_list`, and `measure_args_list`, do not have same lengths.")

  if (!(setequal(names(trait_list), names(measure_list)) &&
        setequal(names(measure_list), names(measure_args_list)))) stop("`trait_list`, `measure_list`, and `measure_args_list`, names do not match.")

  # Set up ------------------------------------------------------------------

  # number of individuals
  n_g <- nrow(trait_list[[1]])

  # if initial weights were not provided, need to assign
  if (is.null(initial_weights)) {

    # n_g = nrow(trait_list[[1]]) # how many individuals were provided in the trait dataset

    if (n_t > n_g) stop("'n_t' is larger than the number of individuals in trait data.")

    initial_weights <- propose_initial_weights(n_g = n_g, n_t = n_t, w_max = weights_max)

  }

  # get initial measure values
  measure_out = calculate_measure(trait_list, measure_list, measure_args_list, w = initial_weights)

  if (is.na(measure_out) || is.null(measure_out)) stop("Initial measure value is NA or NULL.")

  # begin simulated annealing ----------------------------------------------

  if (save_chain) {
    # allocate chain
    chain <- list(
      weight = matrix(NA_real_, max_steps, n_g),
      values = matrix(NA_real_, max_steps, objectives, dimnames = list(NULL, names(trait_list))),
      accept = rep(NA, max_steps) #logical(max_steps)
    )

    # add first observations to chain
    chain$weight[1,] <- initial_weights
    chain$values[1,] <- unlist(measure_out)
    # chain$accept[1] <- NA # not needed since it is already NA

  }

  # set up archive if needed
  if (nda) {

    archive = list(
      archive_summary = matrix(unlist(measure_out), nrow = 1, dimnames = list(NULL, names(trait_list))),
      archive_weights = matrix(initial_weights, nrow = 1)
    )

  }

  # set preliminary params

  s <- 2 # time step. Starting at 2 since we have already completed s = 1
  wts <- initial_weights

  if(verbose) message("Beginning simulation")

  while ( s <= max_steps ) {

    # get current temperature
    temp <- temp_scheduler(s, max_steps, max_t, nda)

    # make new mixture
    weights_mod <- modify_weights(wts, w_min = weights_min, w_max = weights_max, n_t = n_t)

    # and see what they measure
    measure_mod = calculate_measure(trait_list, measure_list, measure_args_list, w = weights_mod)

    if (is.na(measure_mod) || is.null(measure_mod)) stop("Proposed measure value is NA or NULL.")

    # test if we accept the new mix
    acceptance = accept_reject(
      summary = measure_out,
      proposal_summary = measure_mod,
      t = temp,
      p_depends_delta = p_depends_delta,
      c = acceptance_multipliers,
      c_all = acceptance_multiplier_all_worse
    )

    if(is.null(acceptance) || is.na(acceptance))

    if (nda) { # this deviates from OptGenMix: moved this outside of accept_proposal so that we can explore un-accepted front space

      archive <- clean_archive(measure_mod, weights_mod, archive)

      if (nrow(archive$archive_summary) >= nd_samples) {
        if(verbose) message("\nMaximum archive values reached. Stopping simulation")
        break # dont continue simulation if we have maxed out the archive
      }
    }

    # if we accept mix, store it for use in next round
    if (acceptance) {

      wts = weights_mod
      measure_out = measure_mod

    }

    # store chain information
    if (save_chain) {

    chain$weight[s,] <- wts
    chain$values[s,] <- unlist(measure_out, use.names = FALSE)
    chain$accept[s] <- acceptance

    }

    if(verbose) cat("\rFinished", s, "of", max_steps)

    s <- s + 1

  }

  final_selection = list(
    measure_summary = unlist(measure_out),
    individs_selected = wts
  )

  return(

    list(

      final_selection = final_selection,
      chain = if(save_chain) chain else NULL,
      archive = if (nda) archive else NULL

    )
  )
}




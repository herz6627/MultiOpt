



multiopt_sa <- function(
    trait_list, # list of trait data as vectors or matrices
    measure_list, # list of measure functions corresponding to trait order.
    measure_args_list, # list of any args needed for measures
    n_t = NULL, # how many total "individuals" to select
    weights_max = NULL, # single value or vector of length matching number of individuals in trait dataset
    weights_min = NULL,
    initial_weights = NULL,
    max_steps = 10000,
    max_t = 1,
    min_t = 0, # only needed for multi-optimization,
    p_depends_delta = F, # whether to account for the magnitude in the delta between each measure
    multiplier_list, # instead of c1, c2, etc
    nda = F, # whether to test for non-dominant archives. only needed for multi-opt
    nd_samples = 100 # if nda == T, how ,many samples to conduct
) {


  # Run checks --------------------------------------------------------------

  if (!all(vapply(measure_list, is.function, logical(1)))) {
    stop("All supplied measures must be functions")
  }


  if(any(!sapply(trait_list, is.matrix))) stop("All trait data within 'trait_list' must be in matrix format.")

  if (is.null(n_t) && is.null(initial_weights)) stop("Either 'n_t' or 'initial_weights' must be provided.")

  if ( !is.null(n_t) & !is.null(initial_weights)) {

    if (n_t != sum(initial_weights)) stop( "Conflict in number of target individuals and initial weights provided \n" )

  }


  # Set up ------------------------------------------------------------------

  # if initial weights were not provided, need to assign
  if (is.null(initial_weights)) {

    n_g = nrow(trait_list[[1]])

    if (n_t > n_g) stop("'n_t' is larger than the number of individuals in trait data.")

    initial_weights <- propose_initial_weights(n_g = n_g, n_t = n_t, w_max = weights_max)

  }

  # get initial measure values




  # set up non-dominant archive (nda)
  if (nda) {

    archive = list()

    nda_complete = FALSE

  } else  {

    nda_complete = FALSE

  }


  return(initial_weights)
}




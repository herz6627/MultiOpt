



multiopt_sa <- function(
    trait_list, # list of trait data as vectors or matrices
    measure_list, # list of measures corresponding to trait order. must match with supplied list.
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


  validate_measures <- function(measures) {
    allowed = c("nei", "shannon", "allele_enrichment", "vector_weighted_mean",
                "negative_vector_weighted_mean", "sum_squared_difference",
                "negative_sum_squared_difference", "vector_diff_weighted_mean",
                "negative_vector_diff_weighted_mean", "matrix_weighted_mean",
                "negative_matrix_weighted_mean")


    bad_idx <- which(!measures %in% allowed)

    if (length(bad_idx) > 0) {
      stop(
        paste0(
          "Invalid measure at position(s): ",
          paste(bad_idx, collapse = ", "),
          "\nValue(s): ",
          paste(measures[bad_idx], collapse = ", "),
          "\nAllowed: ",
          paste(allowed, collapse = ", ")
        ),
        call. = FALSE
      )
    }

    invisible(TRUE)
  }

  validate_measures(measure_list)

  if(any(!sapply(trait_list, is.matrix))) stop("All trait data within 'trait_list' must be in matrix format.")



  if (is.null(n_t) && is.null(initial_wights)) stop("Either 'n_t' or 'initial_weights' must be provided.")

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




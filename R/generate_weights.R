

#' Internal helper for validating weight constraints.
#'
#' @param n_g Integer. Number of genotypes/individuals.
#' @param n_t Integer. Total target population size.
#' @param w_min Optional integer vector of minimum allowable weights.
#'    May be:
#'    \itemize{
#'    \item NULL (defaults to 0).
#'    \item a scalar applied to all individuals.
#'    \item a vector of length n_g (number of individuals with trait data) indicating minimum weight for each individual.
#'    }
#' @param w_max Optional integer vector of maximum allowable weights.
#'  May be:
#'   \itemize{
#'     \item NULL (unconstrained)
#'     \item a scalar applied to all individuals.
#'     \item a vector of length n_g (number of individuals with trait data)
#'   }
#' @returns List of length 2 containing a vector of length n_g for both w_min and w_max
#' @export
.check_weights <- function(
    n_g,
    n_t = NULL,
    w_min,
    w_max
) {

  # maxima
  if (!is.null(w_max)) {

    if (length(w_max) == 1) w_max <- rep(w_max, n_g)     # single value: recycle

    if (length(w_max) != n_g) stop("w_max must have length 1 or n_g")

    if (any(w_max < 0)) stop("w_max cannot contain negative values")

    if (sum(w_max) < n_t) stop("Infeasible: sum(w_max) < n_t")

  } else {

    if (is.null(n_t)) {

      w_max <- rep(Inf, n_g)

    } else {

      w_max <- rep(n_t, n_g) # dont really need to limit by n_t but is a good catch if something goes wrong elsewhere

    }

  }

  # minima
  if (!is.null(w_min)) {

    if (length(w_min) == 1) w_min <- rep(w_min, n_g) # single value: recycle

    if (length(w_min) != n_g) stop("w_min must have length 1 or n_g")

    if (any(w_min < 0)) stop("w_min cannot contain negative values")

    if (sum(w_min) > n_t) stop("Infeasible: sum(w_min) > n_t")

    if (any(w_min > w_max)) stop("w_min cannot be > w_max")

  } else {

    w_min <- rep(0, n_g)

  }

  return(list(
    w_max = w_max,
    w_min = w_min
  ))

}



#' Generate initial weights
#'
#' Generates initial weights for use before initiating simulated annealing.
#'
#' @inheritParams .check_weights
#' @param verbose Logical. Print progress messages.
#'
#' @returns Integer vector of weights
#' @export
propose_initial_weights <- function(
    n_g,
    n_t,
    w_max = NULL, # same as weights_max
    w_min = NULL,
    verbose = FALSE
) {

  if (verbose) {
    cat("Generating constrained initial weights\n")
  }


  # unconstrained ------------------------------------------------------------
  # if no maximum weights are assigned, we'll just randomly pick n_t individuals however many times
  if (is.null(w_max) && is.null(w_min)) {

    return(
      as.vector(
        rmultinom(
          n = 1,
          size = n_t,
          prob = rep(1 / n_g, n_g)
        )
      )
    )

  }

  # check provided weights ----------------------------------------------------

  cleaned_weights = .check_weights(
    n_g = n_g,
    n_t = n_t,
    w_min = w_min,
    w_max = w_max
  )

  w_min <- cleaned_weights$w_min
  w_max <- cleaned_weights$w_max

  # constrained -------------------------------------------------------------


  # assign minima
  w <- w_min

  # already complete
  if ((n_t - sum(w)) == 0) {
    return(w)
  }

  # remaining available weights
  remaining <- w_max - w

  # assign remaining weights

  for (k in seq_len(n_t - sum(w))) {

    # indices still able to receive weight
    available <- which(remaining > 0)

    if (length(available) == 0) {
      stop("Unexpected allocation")
    }

    # sample one feasible individual
    i <- sample(available, size = 1)

    w[i] <- w[i] + 1
    remaining[i] <- remaining[i] - 1
  }

  return(w)
}


#' Modify weight vector
#'
#' Randomly selects a weight and moves it to a new individual, if possible, while preserving limitations.
#'
#' @inheritParams .check_weights
#' @param w Integer vector of current weights.
#'
#' @returns Modified weight vector of same length as w.
#' @export
modify_weights <- function(
    w,
    w_min = NULL,
    w_max = NULL
) {

  n_g <- length(w)

  # run checks
  cleaned_weights = .check_weights(
    n_g = n_g,
    n_t = n_t,
    w_min = w_min,
    w_max = w_max
  )

  w_min <- cleaned_weights$w_min
  w_max <- cleaned_weights$w_max

  # options to modify
  dec_candidates <- which(w > w_min)   # can decrease
  inc_candidates <- which(w < w_max)   # can increase

  if (length(dec_candidates) == 0) stop("No weights can be decreased")
  if (length(inc_candidates) == 0) stop("No weights can be increased")

  # pick which individuals to change
  i_dec <- sample(dec_candidates, size = 1) # who to decrease

  # prevent null move unless unavoidable
  inc_candidates <- setdiff(inc_candidates, i_dec)

  if (length(inc_candidates) == 0) {
    inc_candidates <- i_dec

    warning("Same individual has been decreased/increased")
  }

  i_inc <- sample(inc_candidates, size = 1) # who to increase

  # apply change
  w_new <- w

  w_new[i_dec] <- w_new[i_dec] - 1
  w_new[i_inc] <- w_new[i_inc] + 1

  return(w_new)
}





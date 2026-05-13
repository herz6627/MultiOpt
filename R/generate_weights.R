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

  # maxima
  if (!is.null(w_max)) {

    if (length(w_max) == 1) w_max <- rep(w_max, n_g)     # single value: recycle

    if (length(w_max) != n_g) stop("w_max must have length 1 or n_g")

    if (any(w_max < 0)) stop("w_max cannot contain negative values")

    if (sum(w_max) < n_t) stop("Infeasible: sum(w_max) < n_t")

  }

  # minima
  if (!is.null(w_min)) {

    if (length(w_min) == 1) w_min <- rep(w_min, n_g) # single value: recycle

    if (length(w_min) != n_g) stop("w_min must have length 1 or n_g")

    if (any(w_min < 0)) stop("w_min cannot contain negative values")

    if (sum(w_min) > n_t) stop("Infeasible: sum(w_min) > n_t")

    if (w_min > w_max) stop("w_min cannot be > w_max")

  } else {

    w_min <- rep(0, n_g)

  }


  # constrained -------------------------------------------------------------


  # assign minima
  w <- w_min

  # already complete
  if ((n_t - sum(w)) == 0) {
    return(w)
  }

  # remaining
  remaining <- w_max - w


  # assign remaining weights

  for (k in seq_len(n_t)) {

    # indices still able to receive weight
    available <- which(remaining > 0)

    if (length(available) == 0) {
      stop("Unexpected allocation")
    }

    # sample one feasible individual
    i <- sample(available, size = 1)

    w[i] <- w[i] + 1
  }

  return(w)
}

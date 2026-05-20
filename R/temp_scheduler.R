#' Determine annealing temperature
#'
#'
#' @param s Integer. Current time step.
#' @param max_steps Numeric scalar. Maximum allowed time steps.
#' @param max_t Numeric scalar. Maximum allowed temperature.
#' @param min_t Numeric scalar. Minimum allowed temperature. Only relevant when nda = TRUE.
#' @param nda Logical. If TRUE, indicates to hold temp at minimum temperature to allow exploration of Pareto front.
#'
#' @returns A numeric scalar of current annealing temperature.
#' @export
temp_scheduler <- function(s, max_steps, max_t, min_t = 0, nda = FALSE) {

  proportion_ahead = 1 - s/max_steps

  t <- proportion_ahead * max_t

  # adjust temp if we are exploring nd solutions
  if (nda) {

    if (t < min_t) {

      t <- min_t

    }
  }

  return(t)
}

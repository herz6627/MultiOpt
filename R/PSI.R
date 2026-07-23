
# calculate the shape of the Pareto using the Pareto Shape Index (PSI)
# as described in Unal, M., Warn, G. P., & Simpson, T. W. (2017). Quantifying the Shape of Pareto Fronts During Multi-Objective Trade Space Exploration. Journal of Mechanical Design, 140(021402). https://doi.org/10.1115/1.4038005



#' Pareto Shape Index (PSI)
#'
#' The Pareto Shape Index (PSI; \eqn{\gamma}) describes the shape of the Pareto front
#' (ranging from convex, linear, or concave) for
#' a two-objective comparison of traits. The PSI quantifies
#' the geometric shape of a Pareto front by averaging the vertical distances of
#' pairwise solution-line intersections from the center of the normalized
#' objective space.
#'
#' This index ranges from -1 to 1,
#' where -1 indicates a concave Pareto front (weak trade off),
#' 0 indicates a linear front (strong trade off),
#' and 1 indicates a convex front (which is not desirable when we are
#' maximizing our objectives as this is a very harsh trade off).
#' When there is no trade off between objectives, then the Pareto front is
#' formed by a single solution; therefore, the value of Ncr
#' in Equation (6) in Unal et al (2017) becomes zero, leading the Pareto Shape Index to have a value of infinity.
#'
#' Trait measures must be in the same direction (minimizing or maximizing).
#' This should be automatically accounted for when using the MultiOpt measure functions
#' as MultiOpt simulated annealing attempts to maximize objectives
#' (measure functions will add a negative for minimizations).
#' Thus it is important to not do any transformations to the Pareto outputs before calculating PSI.
#'
#' @param pareto_front `archive_summary` output (non-dominated archive) from a
#' simulated annealing run where rows are non-dominated archive values and a column for each trait.
#' Only two traits at a time are possible.
#'
#' @references Unal, M., Warn, G. P., & Simpson, T. W. (2017).
#' Quantifying the Shape of Pareto Fronts During Multi-Objective Trade Space Exploration.
#' Journal of Mechanical Design, 140(021402). https://doi.org/10.1115/1.4038005
#'
#' @returns Single numeric value of the Pareto Shape Index (PSI).
#'
#' @details
#' The Pareto Shape Index is calculated following Equations (4)--(6) of
#' Unal et al. (2017). For each pair of Pareto-optimal solutions
#' \eqn{i} and \eqn{j}, the intersection of the corresponding solution lines
#' between objectives \eqn{k} and \eqn{l} is calculated as
#'
#' \deqn{
#' Y_{ij,kl} =
#' \frac{
#' (S_{i,k}S_{j,l}) - (S_{i,l}S_{j,k})
#' }{
#' S_{i,k} - S_{i,l} - S_{j,k} + S_{j,l}
#' }
#' }{
#' Yij = (Si,k * Sj,l - Si,l * Sj,k) /
#'       (Si,k - Si,l - Sj,k + Sj,l)
#' }
#'
#' where \eqn{S_{i,k}} and \eqn{S_{i,l}} are the normalized objective values
#' for solution \eqn{i}. Trait values should be normalized to the interval
#' [0, 1] and transformed to a common optimization direction (minimizing or maximizing) prior to
#' computing the index.
#'
#' The signed distance of each intersection from the center of the normalized
#' objective space is then
#'
#' \deqn{
#' P_{ij,kl} = Y_{ij,kl} - 0.5
#' }{
#' Pij = Yij - 0.5
#' }
#'
#' Finally, the Pareto Shape Index is calculated as the mean signed distance
#' across all pairwise solution intersections:
#'
#' \deqn{
#' \gamma_{kl} =
#' \frac{1}{N_{cr}}
#' \sum_{i=1}^{N_p-1}
#' \sum_{j=i+1}^{N_p}
#' P_{ij,kl}
#' }{
#' gamma = mean(Pij)
#' }
#'
#' where \eqn{N_p} is the number of Pareto-optimal solutions and
#' \eqn{N_{cr} = N_p(N_p-1)/2} is the number of pairwise solution
#' intersections. The resulting index ranges from -1 to 1, with
#' negative values indicating concave Pareto fronts, values near zero
#' indicating approximately linear tradeoffs, and positive values
#' indicating convex Pareto fronts.
#'
#' @export
pareto_shape_index <- function(pareto_front) {
  # pareto_front = $archive_summary
  # rows = Pareto solutions
  # columns = two objectives
  # objectives are already converted so that both are in the same direction (e.g., both maximized or both minimized)
  # objective values are normalized to 0–1



  pareto_front <- as.matrix(pareto_front)

  if(any(pareto_front > 1) | any(pareto_front < 0)) stop("Pareto values should be using scaled [0-1] trait data.")

  if (ncol(pareto_front) != 2) {
    stop("Pareto Shape Index is defined for two objectives.")
  }

  n <- nrow(pareto_front)

  if (n < 3) {
    warning("Very few Pareto (<3) solutions provided, this will likely cause bias.")
  }

  # number of crossings
  n_cross <- choose(n, 2)

  if (n_cross == 0) {
    return(Inf)
  }

  # store crossing deviations
  P <- numeric(n_cross)

  counter <- 1

  for (i in 1:(n-1)) {

    for (j in (i+1):n) {

      Si_k <- pareto_front[i,1]
      Si_l <- pareto_front[i,2]

      Sj_k <- pareto_front[j,1]
      Sj_l <- pareto_front[j,2]

      denominator <- Si_k - Si_l - Sj_k + Sj_l

      # parallel lines do not have an intersection
      if (abs(denominator) < 1e-12) {
        next
      }

      # Equation 4
      Yij_kl <- ((Si_k * Sj_l) - (Si_l * Sj_k)) / denominator

      # Equation 5
      P[counter] <- 0.5 - Yij_kl

      counter <- counter + 1
    }
  }

  # remove unused entries (parallel lines)
  P <- P[!is.na(P)]

  if(length(P) == 0) {
    return(NA_real_)
  }

  # Equation 6
  gamma <- mean(P)

  if (gamma > 1 | gamma < -1) stop("PSI should be between 0 and 1. There is an error somewhere.")

  return(gamma)
}







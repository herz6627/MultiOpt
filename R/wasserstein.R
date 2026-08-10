#' Calculate the Wasserstein Distance Between Two Empirical Distributions
#'
#' Calculates the Wasserstein distance between two one-dimensional empirical
#' distributions. The Wasserstein distance quantifies the amount of movement
#' required to transform one distribution into another. Smaller values indicate
#' more similar distributions.
#'
#' @param observations_a Numeric vector containing observations from the first
#' distribution (e.g., source population trait values).
#'
#' @param observations_b Numeric vector containing observations from the second
#' distribution (e.g., selected collection trait values).
#'
#' @param p Numeric value specifying the Wasserstein order. The default is
#' \code{p = 1}, corresponding to the first Wasserstein distance. Higher values
#' increasingly weight larger distributional differences.
#'
#' @param weights_a Optional weights for observations in observations_a.
#'
#' @param weights_b Optional weights for observations in observations_b.
#'
#' @return Numeric value representing the p-Wasserstein distance. Smaller values
#' indicate distributions are more similar.
#'
#' @references D. Schuhmacher, B. Bähre, N. Bonneel, C. Gottschlich, V. Hartmann, F. Heinemann, B. Schmitzer and J. Schrieber (2024). transport: Computation of Optimal Transport Plans and Wasserstein Distances. R package version 0.15-0. https://cran.r-project.org/package=transport
#'
#' Vallender, S. S. (1974). Calculation of the Wasserstein Distance Between Probability Distributions on the Line. Theory of Probability & Its Applications. https://doi.org/10.1137/1118101
#'
#'
#' @details
#' For one-dimensional distributions (\code{p = 1}), the p-Wasserstein distance is calculated
#' from the quantile functions of the two distributions:
#'
#' \deqn{
#' W_p(F,G) =
#' \left(\int_0^1 |Q_F(u)-Q_G(u)|^p du\right)^{1/p}
#' }
#'
#' where \eqn{Q_F} and \eqn{Q_G} are the quantile functions of distributions
#' \eqn{F} and \eqn{G}, respectively.
#'
#'
#' For equally weighted distributions with the same number of observations,
#' this simplifies to the mean p-th powered difference between corresponding
#' ordered observations:
#'
#' \deqn{
#' W_p =
#' \left[
#' \frac{1}{n}\sum_{i=1}^{n}
#' |x_{(i)}-y_{(i)}|^p
#' \right]^{1/p}.
#' }
#'
#'
#' When \code{p = 1}, the metric represents the average distance that trait
#' values must be moved for one distribution to match the other.
#' The equation also simplifies to:
#'
#' \deqn{
#' W_p =
#' \left[
#' \frac{1}{n}\sum_{i=1}^{n}
#' |x_{(i)}-y_{(i)}|
#' \right].
#' }
#'
#'
#' This metric can be sensitive to outliers because it tracks raw displacements
#' rather than overlap.
#'
#' @examples
#' # Compare a population trait distribution to a selected collection
#'
#' population_trait <- rnorm(1000)
#' collection_trait <- sample(population_trait, 50)
#'
#' wasserstein(
#'   population_trait,
#'   collection_trait
#' )
#'
#' @export
wasserstein <- function(observations_a,
                        observations_b,
                        p = 1,
                        weights_a = NULL,
                        weights_b = NULL) {


  # ---- Input checks ----

  stopifnot(
    is.numeric(observations_a),
    is.numeric(observations_b),
    length(observations_a) > 0,
    length(observations_b) > 0,
    length(p) == 1,
    p >= 1
  )


  # ---- Remove missing observations ----

  if (anyNA(observations_a)) {

    keep <- !is.na(observations_a)

    observations_a <- observations_a[keep]

    if (!is.null(weights_a)) {
      weights_a <- weights_a[keep]
    }
  }


  if (anyNA(observations_b)) {

    keep <- !is.na(observations_b)

    observations_b <- observations_b[keep]

    if (!is.null(weights_b)) {
      weights_b <- weights_b[keep]
    }
  }


  # ---- Fast calculation for equal-sized unweighted samples ----
  #
  # When both distributions contain the same number of observations and
  # observations have equal weight:
  #
  # Wp = [ mean(|x(i)-y(i)|^p) ]^(1/p)
  #
  # where x(i) and y(i) are sorted observations.
  #
  # When p = 1 the equation simplifies to
  # Wp = [ mean(|x(i)-y(i)| ]


  if (length(observations_a) == length(observations_b) &&
      is.null(weights_a) &&
      is.null(weights_b)) {

    return(
      mean(
        abs(
          sort(observations_a) -
            sort(observations_b)
        )^p
      )^(1/p)
    )
  }


  # ---- Prepare observation weights ----
  #
  # Convert missing weights to equal weighting.
  # Remove observations with zero contribution.

  prepare_distribution <- function(values, weights) {

    if (is.null(weights)) {

      return(
        list(
          values = values,
          weights = rep(1, length(values))
        )
      )
    }

    stopifnot(length(values) == length(weights))

    keep <- weights > 0

    list(
      values = values[keep],
      weights = weights[keep]
    )
  }


  observations_a_data <- prepare_distribution(
    observations_a,
    weights_a
  )

  observations_b_data <- prepare_distribution(
    observations_b,
    weights_b
  )


  observations_a <- observations_a_data$values
  weights_a <- observations_a_data$weights

  observations_b <- observations_b_data$values
  weights_b <- observations_b_data$weights


  n_a <- length(observations_a)
  n_b <- length(observations_b)


  # ---- Sort observations and weights ----

  order_a <- order(observations_a)
  order_b <- order(observations_b)

  observations_a <- observations_a[order_a]
  weights_a <- weights_a[order_a]

  observations_b <- observations_b[order_b]
  weights_b <- weights_b[order_b]


  # ---- Calculate cumulative distribution probabilities ----
  #
  # These represent the jump locations of the empirical cumulative
  # distribution functions.

  cumulative_a <- cumsum(
    weights_a / sum(weights_a)
  )[-n_a]

  cumulative_b <- cumsum(
    weights_b / sum(weights_b)
  )[-n_b]


  # ---- Determine quantile overlap ----
  #
  # Each observation contributes over an interval of cumulative probability.
  # Determines how long each interval of the empirical CDF lasts.

  observations_a_repeats <- hist(
    cumulative_b,
    breaks = c(-Inf, cumulative_a, Inf),
    plot = FALSE
  )$counts + 1

  observations_b_repeats <- hist(
    cumulative_a,
    breaks = c(-Inf, cumulative_b, Inf),
    plot = FALSE
  )$counts + 1


  # ---- Expand observations over quantile intervals ----
  # so that the quantile values match

  observations_a_expanded <- rep(
    observations_a,
    times = observations_a_repeats
  )

  observations_b_expanded <- rep(
    observations_b,
    times = observations_b_repeats
  )


  # ---- Integrate distance between quantile functions ----
  #
  # Wp = [ integral_0^1 |Qa(u)-Qb(u)|^p du ]^(1/p)

  probability_breaks <- sort(
    c(
      cumulative_a,
      cumulative_b
    )
  )

  probability_lower <- c(
    0,
    probability_breaks
  )

  probability_upper <- c(
    probability_breaks,
    1
  )

  wasserstein_distance <- sum(
    (probability_upper - probability_lower) *
      abs(
        observations_a_expanded - observations_b_expanded
      )^p
  )^(1/p)

  return(wasserstein_distance)
}

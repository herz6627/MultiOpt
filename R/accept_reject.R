#' Determines acceptance or rejection of a proposal
#'
#' @param summary           List of current objective values.
#' @param proposal_summary  List of proposed objective values.
#' @param t                 Current temperature.
#' @param p_depends_delta   If TRUE, acceptance depends on magnitude of worsening.
#' @param c                 Penalty multipliers. Single value of vector of same length as `summary`.
#' @param c_all             Multiplier when all objectives worsen.
#'
#' @return Logical scalar indicating whether proposal is accepted.
#' @export

accept_reject <- function(
    summary,
    proposal_summary,
    t,
    p_depends_delta = FALSE,
    c = 1,
    c_all = 1
) {

  # convert lists to numeric vectors
  current  <- unlist(summary, use.names = FALSE)
  proposal <- unlist(proposal_summary, use.names = FALSE)

  # checks
  if (length(current) != length(proposal)) {
    stop("'summary' and 'proposal_summary' must have equal length.")
  }

  if (anyNA(current) || anyNA(proposal)) {
    stop("NA detected in current or proposal measure values.")
  }

  if (!is.numeric(current) || !is.numeric(proposal)) {
    stop("All summary values must be numeric")
  }

  if(c_all != 1 && !is.numeric(c_all)) stop("c_all must be a single numeric value.")

  # recycle scalar penalties
  if (length(c) == 1) {
    c <- rep(c, length(current))
  }

  if (length(c) != length(current)) {
    stop("'c' must have length 1 or number of objectives.")
  }

  # check improvement
  better <- proposal > current
  worse  <- proposal < current
  equal  <- proposal == current

  n_better <- sum(better)
  n_worse  <- sum(worse)

  # all improvement
  # proposal >= current in all objectives and strictly better in at least one
  if (all(!worse) && any(better)) {
    return(TRUE)
  }

  # how much worse
  delta <- pmax(current - proposal, 0) # forces improvements to 0. This logic matches OptGenMix.

  # all worse
  if (all(worse)) {

    # penalty numerator
    if (p_depends_delta) {

      numerator <- sum(delta * c)

    } else {

      numerator <- sum(c)

    }

    numerator <- numerator * c_all

  } else {

    # some worse

    if (p_depends_delta) {

      numerator <- sum(delta * c)

    } else {

      # penalty depends only on number of worse objectives
      numerator <- sum(worse * c)

    }


  }

  # acceptance probability
  p_accept_worse <- exp(-numerator / t)

  out <- runif(1) < p_accept_worse # stochastic acceptance. Standard to Metropolis algorithm method.

  if (is.null(out) || is.na(out)) stop("Acceptance value is NA or NULL.")

}

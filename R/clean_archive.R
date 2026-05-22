#' Reconcile a candidate solution with a non-dominated (Pareto) archive
#'
#' Updates a Pareto archive by comparing a new mixture solution against
#' existing archived solutions across multiple objectives. The function
#' removes any archived solutions dominated by the new mixture and rejects
#' the candidate if it is dominated by any existing archive member.
#'
#' Dominance is defined in the Pareto sense:
#' a solution A dominates B if A is greater than or equal to B in all
#' objectives and strictly greater in at least one objective.
#'
#' Exact duplicates (within objective space) are treated as redundant and are
#' not added to the archive.
#'
#' @param summary Named list or numeric vector of objective values for the
#'   candidate solution. If provided as a list, it will be coerced to a numeric
#'   vector via `unlist()`.
#'
#' @param wts Numeric vector decision variables associated with
#'   the candidate solution.
#'
#' @param archive List of existing non-dominated solutions. Each element is a
#'   list containing:
#' \describe{
#'   \item{archive_summary}{Numeric matrix of objective values, where each row
#'   corresponds to a solution and each column to an objective.}
#'
#'   \item{archive_weights}{Numeric matrix of decision variables, where each row
#'   corresponds to a solution and each column an individual.}
#' }
#'
#' @return Updated archive list containing only non-dominated solutions after
#'   reconciliation with the candidate.
#'
#' @details
#' The function performs three operations:
#' \enumerate{
#'   \item Checks whether the candidate is dominated by any archived solution.
#'   \item Identifies and removes archived solutions dominated by the candidate.
#'   \item Inserts the candidate if it is not dominated.
#' }
#'
#'
#' @examples
#' archive <- list(
#'   archive_summary = matrix(c(1, 2, 2, 1), ncol = 2, byrow = TRUE),
#'   archive_weights = matrix(c(0.1, 0.2, 0.3, 0.4), ncol = 2, byrow = TRUE)
#' )
#'
#' candidate_summary <- list(x = 2, y = 2)
#' candidate_weights <- c(0.5, 0.6)
#'
#' archive <- reconcile_sample_nondominated_archive(
#'   summary = candidate_summary,
#'   wts = candidate_weights,
#'   archive = archive
#' )
#' @export
clean_archive <- function(
    summary,  # current measure summary list
    wts, # current vector of weights
    archive # current archive
    ) {

  S <- archive$archive_summary
  W <- archive$archive_weights

  rownames(S) <- NULL
  rownames(W) <- NULL


  if (!is.matrix(S) & !is.matrix(W)) {
    stop("`archive$archive_summary` and `archive$weights` need to be matrices.")
  }

  # ensure candidate is numeric vector
  new_vec <- unlist(summary, use.names = FALSE)

  # handle empty archive case
  if (is.null(S) || nrow(S) == 0) {

    archive$archive_summary <- matrix(new_vec, nrow = 1)
    archive$archive_weights <- matrix(wts, nrow = 1)

    return(archive)
  }

  # check for exact duplicates
  if (nrow(S) > 0) {
    is_duplicate <- apply(S, 1, function(old) {
      all(old == new_vec)
    })

    if (any(is_duplicate)) {
      return(archive)
    }
  }

  # Check if new is dominated
  dominates_new <- apply(S, 1, function(old) {
    all(old >= new_vec) && any(old > new_vec)
  })

  if (any(dominates_new)) {
    return(archive)  # reject new mixture immediately
  }

  # Find which old solutions are dominated by new
  dominated_old <- apply(S, 1, function(old) {
    all(new_vec >= old) && any(new_vec > old)
  })

  # Filter archive
  S_new <- S[!dominated_old, , drop = FALSE]
  W_new <- W[!dominated_old, , drop = FALSE]

  # Append new solution
  archive$archive_summary <- rbind(S_new, unname(new_vec))
  archive$archive_weights <- rbind(W_new, wts)

  return(archive)
}






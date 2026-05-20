#' Combine two archives and reduce if needed
#'
#' Combines two archives from `multiopt_sa` and reduces to one archive,
#' filtering to only non-dominated values. Both archives should be produced
#' from the same simulated annealing arguments and trait data sets.
#' Order of which archive value is assigned to archive1 or archive2
#' is arbitrary and won't make a difference.
#'
#' @param archive1 Archive from `multiopt_sa` (out1$archive)
#' @param archive2 Archive from `multiopt_sa` (out2$archive)
#'
#' @returns a single archive list matching archive output style as multiopt_sa.
#' A list containing:
#' \describe{
#'   \item{archive_summary}{Values for each non-dominated solution with a column for each trait.}
#'   \item{archive_weights}{Weights for each non-dominated solution (row) with a column for each individual.}
#' }
#'
#'
#' @export
combine_archives <- function(archive1, archive2) {

  anew   <- archive1

  anew$archive_summary <- rbind(archive1$archive_summary, archive2$archive_summary)
  anew$archive_weights <- rbind(archive1$archive_weights, archive2$archive_weights)

  anew_len <- nrow(anew$archive_summary)

  i_dom <- NULL

  for (i in 1:anew_len) {

    iv1 <- anew$archive_summary[i,1]
    iv2 <- anew$archive_summary[i,2]


    for (j in 1:anew_len) {

      iv1_smaller <- FALSE
      iv2_smaller <- FALSE

      if ( i != j ) {
        jv1 <- anew$archive_summary[j,1]
        jv2 <- anew$archive_summary[j,2]

        if (iv1 < jv1) {iv1_smaller <- TRUE}
        if (iv2 < jv2) {iv2_smaller <- TRUE}

        if ( iv1_smaller & iv2_smaller ) {
          i_dom <- c(i_dom, i)
          i_dom <- unique(i_dom)
        }

      }

    }

  }

  if (length(i_dom) > 0) {
    anew$archive_summary  <- anew$archive_summary[-i_dom,]
    anew$archive_weights <- anew$archive_weights[-i_dom,]
  }

  return(anew)

}

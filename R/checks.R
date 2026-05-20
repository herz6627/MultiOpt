#' Check trait data scaling
#'
#' @param x List of trait values
#' @param tolerance How much variation in trait values are allowed.
#' Compared to the ratio of max(difference in trait values) / min(difference in trait values)
#'
#' @returns An error or nothing.
#' @export
check_similar_scale <- function(x, # trait list
                                tolerance = 2) {

  on.exit(options("warn"))
  options(warn=1) # immediately print warnings

  # get size of range for each trait
  ranges <- lapply(x, range)
  diffs <- lapply(ranges, diff)

  # get ratio of differences
  ratio <- max(unlist(diffs)) / min(unlist(diffs))

  if (ratio > tolerance) {
    warning( # could also make this `stop`
      "\n\nTraits appear to be on very different scales. This will likely bias simulation results.\nConsider using `scale_traits()`.\nIgnore if using a genotype matrix.\n",
      "Range widths: ",
      paste(names(diffs), ":",
            signif(unlist(diffs), 3),
            collapse = ", "),
      "\n",
      "Ratio = ", signif(ratio, 3), "\n\n"
    )
  } else {

    # even if the scales are not on very different scales, there could be some minor issues.
    if (any(

      abs(do.call(rbind, ranges)[,1] - 0) > 1e-8 | # adding a slight tolerance instead of using just 0
      abs(do.call(rbind, ranges)[,2] - 1) > 1e-8

    )) {

      warning("\n\nTraits do not appear to be scaled to [0,1]. This could bias simulation results.\n\n")
    }
  }


  invisible(TRUE)
}

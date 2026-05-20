#' Min-max scale trait datasets to the range [0, 1]
#'
#' Applies min-max scaling to each element in a list of trait datasets,
#' transforming values to the range [0, 1]. Scaling is performed
#' independently for each dataset using.
#'
#' @param trait_list A list of matrices (with a single column for a trait vector or a pairwise matrix), or other
#'   numeric objects to be scaled independently.
#'
#' @return A list of scaled trait datasets with the same structure
#'   and names as `trait_list`.
#'
#' @details
#' Each element of `trait_list` is processed independently using
#' min-max normalization:
#'
#' \deqn{
#' x_{scaled} = \frac{x - \min(x)}{\max(x) - \min(x)}
#' }
#'
#' The output preserves the original list structure.
#'
#' If a dataset contains missing values, they are ignored during
#' scaling but retained in the output.
#'
#' @export
scale_traits <- function(trait_list){

  min_max_scale = function(x) {
    if (any(is.na(x))) message("Found NAs. Ignoring.")

    (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))
  }

  trait_list_scaled <- lapply(trait_list, min_max_scale)

  return(trait_list_scaled)
}

#' Back transform output to original scale
#'
#' This function uses the original transformed trait data provided in `multiopt_sa` and the corresponding output
#' from `multiopt_sa` to back transform simulation output to the original trait scale.
#' This function assumes trait data was transformed using `scale_traits`.
#'
#' @param trait_list Original list of trait data as used in `multiopt_sa`.
#' @param multiopt_output Unmodified output from `multiopt_sa`.
#'
#' @details
#' Note that this function takes the absolute value of transformed trait variables.
#'
#'
#' @returns Object of same dimensions and formatting as `multiopt_output`.
#' @export
unscale_output <- function(trait_list, multiopt_output){

  min_max_unscale <- function(x_scaled, x_unscaled) {
    x_scaled * (max(x_unscaled, na.rm = T) - min(x_unscaled, na.rm = T)) + min(x_unscaled, na.rm = T)
  }

  n_traits <- length(trait_list)

  # check everything lines up (this check is not comprehensive)
  if(!all(names(trait_list) %in% names(multiopt_output$final_selection$measure_summary))) stop("Original trait names don't line up with multiopt output.")

  # make output
  new_out = multiopt_output

  for (i in seq(n_traits)) { # for each trait

    # transform final_selection
    new_out$final_selection$measure_summary[i] <-
      min_max_unscale(
        x_scaled = abs(multiopt_output$final_selection$measure_summary[i]),
        x_unscaled = trait_list[[i]]
      )

    # transform chain
    new_out$chain$values[,i] <-
      min_max_unscale(
        x_scaled = abs(multiopt_output$chain$values[,i]),
        x_unscaled = trait_list[[i]]
      )

    # transform archive (if needed)
    if (!is.null(new_out$archive)) {

      new_out$archive$archive_summary[,i] <-
        min_max_unscale(
          x_scaled = abs(multiopt_output$archive$archive_summary[,i]),
          x_unscaled = trait_list[[i]]
        )

    }

  }

  return(new_out)
}

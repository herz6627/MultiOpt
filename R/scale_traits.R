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

  trait_list_scaled <- lapply(trait_list, min_max_scale)

  return(trait_list_scaled)
}

#' Min-max scale
#'
#' Transforms numeric values using min-max normalization, rescaling
#' values to the interval [0, 1] according to:
#'
#' \deqn{
#' x_{scaled} = \frac{x - \min(x)}{\max(x) - \min(x)}
#' }
#'
#' Missing values are ignored when calculating minima and maxima.
#' If missing values are present, a message is produced.
#'
#' @param x A numeric vector, matrix, or array.
#'
#' @return An object with the same structure as `x`, with values
#'   scaled to the range [0, 1].
#'
#' @details
#' Scaling is performed using the global minimum and maximum of `x`.
#' Missing values (`NA`) are retained in the output.
#'
#' If all non-missing values are identical, the function will return
#' `NaN` values due to division by zero.
#'
#' @export
min_max_scale = function(x) {
  if (any(is.na(x))) message("Found NAs. Ignoring.")

  (x - min(x, na.rm = T)) / (max(x, na.rm = T) - min(x, na.rm = T))
}

#' Reverse min-max scaling using a reference dataset
#'
#' Transforms values from a min-max scaled scale back to their
#' original scale using the minimum and maximum of a reference
#' unscaled dataset.
#'
#' The transformation is:
#'
#' \deqn{
#' x = x_{scaled} \times (\max(x_{ref}) - \min(x_{ref})) + \min(x_{ref})
#' }
#'
#' @param x_scaled Numeric vector, matrix, or array containing values
#'   scaled to the range [0, 1].
#' @param x_unscaled Numeric vector, matrix, or array providing the
#'   reference scale used for back-transformation.
#'
#' @return An object with the same structure as `x_scaled`, transformed
#'   back to the scale of `x_unscaled`.
#'
#' @details
#' The minimum and maximum used for back-transformation are calculated
#' from `x_unscaled` with missing values ignored.
#'
#' This function assumes that `x_scaled` was originally generated using
#' min-max scaling based on the same reference distribution.
#'
#' @export
min_max_unscale <- function(x_scaled, x_unscaled) {
  x_scaled * (max(x_unscaled, na.rm = T) - min(x_unscaled, na.rm = T)) + min(x_unscaled, na.rm = T)
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
unscale_multiopt <- function(trait_list, multiopt_output){

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

      new_out$archive <- unscale_archive(trait_list, multiopt_output$archive)

    }

  }

  return(new_out)
}

#' Back transform output to original scale
#'
#' This function uses the original transformed trait data provided in `rand_multiopt` and the corresponding output
#' from `rand_multiopt` to back transform simulation output to the original trait scale.
#' This function assumes trait data was transformed using `scale_traits`.
#'
#' @param trait_list Original list of trait data as used in `rand_multiopt`.
#' @param rand_multiopt_output Unmodified output from `rand_multiopt`.
#'
#' @details
#' Note that this function takes the absolute value of transformed trait variables.
#'
#'
#' @returns Object of same dimensions and formatting as `rand_multiopt_output`.
#' @export
unscale_rand_multiopt <- function(trait_list, rand_multiopt_output){

  n_traits <- length(trait_list)

  # check everything lines up (this check is not comprehensive)
  if(!all(names(trait_list) %in% colnames(rand_multiopt_output$measure_summaries))) stop("Original trait names don't line up with multiopt output.")

  # make output
  new_out = rand_multiopt_output

  for (i in seq(n_traits)) { # for each trait

    # transform final_selection
    new_out$measure_summaries[,i] <-
      min_max_unscale(
        x_scaled = abs(rand_multiopt_output$measure_summaries[,i]),
        x_unscaled = trait_list[[i]]
      )

    # transform archive (if needed)
    if (!is.null(new_out$archive)) {

      new_out$archive <- unscale_archive(trait_list, rand_multiopt_output$archive)

    }

  }

  return(new_out)
}

#' Back transform output to original scale
#'
#' This function uses the original transformed trait data provided in `multiopt_sa` or `rand_multiopt` and the corresponding output
#' from `multiopt_sa` or `rand_multiopt` to back transform archive values to the original trait scale.
#' This function assumes trait data was transformed using `scale_traits`.
#'
#' @param trait_list Original list of trait data as used in `multiopt_sa`.
#' @param archive_output Archive output from `multiopt_sa` (out$archive), `rand_multiopt` (out$archive) or `explore_pareto`.
#'
#' @details
#' Note that this function takes the absolute value of transformed trait variables.
#'
#' @returns List with archive values and weights matching formatting of `archive_output`.
#' @export
unscale_archive <- function(trait_list, archive_output){

  n_traits <- length(trait_list)

  # check everything lines up (this check is not comprehensive)
  if(!all(names(trait_list) %in% colnames(archive_output$archive_summary))) stop("Original trait names don't line up with multiopt output.")

  if (is.null(archive_output)) stop("Archive is NULL")

  # make output
  new_out = archive_output

  # transform archive
  for (i in seq(n_traits)) { # for each trait

    new_out$archive_summary[,i] <-
      min_max_unscale(
        x_scaled = abs(archive_output$archive_summary[,i]),
        x_unscaled = trait_list[[i]]
      )
  }

  return(new_out)
}


#' Back transform output to original scale
#'
#' This function uses the original transformed trait data provided to `singleopt_context` and the corresponding output
#' from `singleopt_context` to back transform archive values to the original trait scale.
#' This function assumes trait data was transformed using `scale_traits`.
#'
#' @param trait_list Original list of trait data as used in `singleopt_context`.
#' @param singleopt_output Archive output from `singleopt_context`.
#'
#' @details
#' Note that this function takes the absolute value of transformed trait variables.
#'
#' @returns List with archive values and weights matching formatting of `singleopt_output`.
#' @export
unscale_singleopt <- function(trait_list, singleopt_output){

  n_traits <- length(trait_list)

  # check everything lines up (this check is not comprehensive)
  if(!all(names(trait_list) %in% names(singleopt_output$measure_summaries))) stop("Original trait names don't line up with multiopt output.")

  if (is.null(singleopt_output)) stop("Archive is NULL")

  # make output
  new_out = singleopt_output

  # transform values
  for (i in seq(n_traits)) { # for each trait element
    for (z in seq(n_traits)) { # for each trait col within the element

      new_out$measure_summaries[[i]][,z] <-
        min_max_unscale(
          x_scaled = abs(singleopt_output$measure_summaries[[i]][,z]),
          x_unscaled = trait_list[[z]]
        )
    }
  }
  return(new_out)
}

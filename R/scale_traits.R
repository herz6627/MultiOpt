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
#' @param skip_traits Character vector of trait names that should not be
#' back transformed.
#'
#' @details
#' Note that this function takes the absolute value of transformed trait variables.
#'
#'
#' @returns Object of same dimensions and formatting as `multiopt_output`.
#' @export
unscale_multiopt <- function(trait_list, multiopt_output, skip_traits = NULL){

  n_traits <- length(trait_list)
  trait_names <- names(trait_list)

  # check everything lines up (this check is not comprehensive)
  if(!all(names(trait_list) %in% names(multiopt_output$final_selection$measure_summary))) stop("Original trait names don't line up with multiopt output.")

  # check skip_traits are valid
  if (!is.null(skip_traits) &&
      !all(skip_traits %in% trait_names)) {
    warning("One or more skip_traits are not present in trait_list.")
  }

  # traits to back transform
  unscale_traits <- setdiff(trait_names, skip_traits)

  # make output
  new_out = multiopt_output

  # back transform final_selection and chain
  for (trait_name in unscale_traits) {

    # find corresponding column in chain
    col <- match(
      trait_name,
      colnames(multiopt_output$chain$values)
    )

    # back transform final_selection
    new_out$final_selection$measure_summary[trait_name] <-
      min_max_unscale(
        x_scaled = abs(
          multiopt_output$final_selection$measure_summary[trait_name]
        ),
        x_unscaled = trait_list[[trait_name]]
      )

    # back transform chain
    new_out$chain$values[, col] <-
      min_max_unscale(
        x_scaled = abs(
          multiopt_output$chain$values[, col]
        ),
        x_unscaled = trait_list[[trait_name]]
      )
  }

  # back transform archive, if present
  if (!is.null(multiopt_output$archive)) {

    new_out$archive <-
      unscale_archive(
        trait_list = trait_list,
        archive_output = multiopt_output$archive,
        skip_traits = skip_traits
      )
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
#' @param skip_traits Character vector of trait names that should not be
#' back transformed.
#'
#' @details
#' Note that this function takes the absolute value of transformed trait variables.
#'
#'
#' @returns Object of same dimensions and formatting as `rand_multiopt_output`.
#' @export
unscale_rand_multiopt <- function(trait_list, rand_multiopt_output, skip_traits = NULL){

  n_traits <- length(trait_list)
  trait_names <- names(trait_list)

  # check everything lines up (this check is not comprehensive)
  if(!all(trait_names %in% colnames(rand_multiopt_output$measure_summaries))) stop("Original trait names don't line up with multiopt output.")

  # check skip_traits are valid
  if (!is.null(skip_traits) &&
      !all(skip_traits %in% trait_names)) {
    warning("One or more skip_traits are not present in trait_list.")
  }

  # traits to back transform
  unscale_traits <- setdiff(trait_names, skip_traits)

  # make output
  new_out = rand_multiopt_output

  # back transform final measure summaries
  for (trait_name in unscale_traits) {

    # find corresponding column
    col <- match(
      trait_name,
      colnames(rand_multiopt_output$measure_summaries)
    )

    new_out$measure_summaries[, col] <-
      min_max_unscale(
        x_scaled = abs(
          rand_multiopt_output$measure_summaries[, col]
        ),
        x_unscaled = trait_list[[trait_name]]
      )
  }

  # back transform archive, if present
  if (!is.null(rand_multiopt_output$archive)) {

    new_out$archive <-
      unscale_archive(
        trait_list = trait_list,
        archive_output = rand_multiopt_output$archive,
        skip_traits = skip_traits
      )
  }

  return(new_out)
}

#' Back transform output to original scale
#'
#' This function uses the original transformed trait data provided in
#' `multiopt_sa` or `rand_multiopt` and the corresponding output
#' from `multiopt_sa` or `rand_multiopt` to back transform archive values to
#' the original trait scale.
#'
#' This function assumes trait data was transformed using `scale_traits()`.
#'
#' @param trait_list Original list of trait data as used in `multiopt_sa`.
#' @param archive_output Archive output from `multiopt_sa` (out$archive),
#' `rand_multiopt` (out$archive) or `explore_pareto`.
#' @param skip_traits Character vector of trait names that should not be
#' back transformed.
#'
#' @details
#' Note that this function takes the absolute value of transformed trait variables.
#'
#' @returns List with archive values and weights matching formatting of `archive_output`.
#' @export
unscale_archive <- function(trait_list, archive_output, skip_traits = NULL){

  if (is.null(archive_output)) stop("Archive is NULL")

  trait_names <- names(trait_list)

  # check everything lines up
  if(!all(trait_names %in% colnames(archive_output$archive_summary))) {
    stop("Original trait names don't line up with multiopt output.")
  }

  # make output
  new_out <- archive_output

  # traits to back transform
  unscale_traits <- setdiff(trait_names, skip_traits)

  # back transform each trait
  for (trait_name in unscale_traits) {

    # find corresponding archive column
    col <- match(
      trait_name,
      colnames(archive_output$archive_summary)
    )

    new_out$archive_summary[, col] <-
      min_max_unscale(
        x_scaled = abs(archive_output$archive_summary[, col]),
        x_unscaled = trait_list[[trait_name]]
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
#' @param skip_traits Optional characted vector specifing if any traits should
#' be excluded from back transformation.
#'
#' @details
#' Note that this function takes the absolute value of transformed trait variables.
#'
#' @returns List with archive values and weights matching formatting of `singleopt_output`.
#' @export
unscale_singleopt <- function(trait_list, singleopt_output, skip_traits = NULL){

  n_traits <- length(trait_list)

  trait_names <- names(trait_list)

  # check that all original traits are present in output
  if(!all(trait_names %in% names(singleopt_output$measure_summaries))) {
    stop("Original trait names don't line up with multiopt output.")
  }

  if (is.null(singleopt_output)) stop("Archive is NULL")

  # traits to unscale
  if(!is.null(skip_traits)) {

    unscale_traits <- setdiff(trait_names, skip_traits)

  } else {

    unscale_traits <- names(trait_list)

  }

  # make output
  new_out = singleopt_output

  # transform values
  for (trait_name in unscale_traits) {

    # Unscale each measure summary
    for (i in seq_along(singleopt_output$measure_summaries)) {

      # Find column based on name
      col <- match(
        trait_name,
        colnames(singleopt_output$measure_summaries[[i]])
      )

      if (is.na(col)) {

        stop(
          "Trait '", trait_name,
          "' was not found in measure_summaries[[", i, "]]."
        )

      }

      # unscale column
      new_out$measure_summaries[[i]][, col] <-
        min_max_unscale(
          x_scaled = abs(
            singleopt_output$measure_summaries[[i]][, col]
          ),
          x_unscaled = trait_list[[trait_name]]
        )
    }
  }

  return(new_out)
}

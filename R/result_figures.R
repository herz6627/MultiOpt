
#' Plot trait value chains over iterations
#'
#' Visualizes the temporal trajectory of trait values through the simulated annealing
#' procedure.
#'
#' @param chain_list Chain output from `multiopt_sa` (out$chain).
#'
#' @return A ggplot object showing trait trajectories over iteration steps,
#'   faceted by trait.
#'
#' @details
#'
#' This is primarily intended for diagnosing convergence, mixing, or
#' temporal dynamics in optimization or sampling chains.
#'
#'
#' @import ggplot2
#' @import dplyr
#' @import patchwork
#' @export
plot_chain = function(
    chain_list
){

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for this function.",
      call. = FALSE
    )
  }

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop(
      "Package 'dplyr' is required for this function.",
      call. = FALSE
    )
  }

  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for this function.",
      call. = FALSE
    )
  }

  chain_dat <- as.data.frame(chain_list$values)

  chain_dat |>
    dplyr::mutate(t = dplyr::row_number()) |>
    tidyr::pivot_longer(cols = -t, names_to = "trait", values_to = "val") |>
    ggplot2::ggplot(ggplot2::aes(y = val, x = t)) +
    ggplot2::geom_line() +
    ggplot2::facet_wrap(~trait, scales = "free_y") +
    ggplot2::labs(
      x = "Step",
      y =  "Trait value",
      title = "Chain") +
    ggplot2::theme_bw()

}


#' Visualize selection frequency across individuals and traits
#'
#' Produces plots showing how frequently individuals are selected across
#' replicates. Supports both single-trait and multi-trait visualization modes.
#'
#' For each individual, selection frequency is computed as the number of
#' replicates in which the individual was selected (ignoring weights).
#'
#' @param individs_selected Matrix or data structure indicating selected
#'   individuals across replicates as produced by `rand_multiopt` (out$individs_selected).
#' @param trait_list Named list of trait data frames, each containing a
#'   single column of trait values. Must not contain multi-column objects
#'   (e.g. pairwise matrices are not supported). This data does not need to be the same supplied to `rand_multiopt`,
#'   but does need to be in the same order for each individual.
#'
#' @return A ggplot2 object (single trait case) or a patchwork object
#'   containing pairwise trait scatterplots (multi-trait case).
#'
#' @details
#' The function operates in two modes:
#'
#' **Single trait case**
#' \itemize{
#'   \item Aggregates trait values with selection frequency.
#'   \item Produces a bar plot ordered by trait value.
#'   \item Fill indicates number of replicates in which each individual
#'         was selected.
#' }
#'
#' **Multi-trait case**
#' \itemize{
#'   \item Constructs a data frame of trait values.
#'   \item Computes selection frequency per individual.
#'   \item Generates all pairwise trait scatterplots.
#'   \item Uses fill color to indicate selection frequency.
#' }
#'
#' Selection frequency is computed as the row-wise count of nonzero
#' entries in `individs_selected`.
#'
#' @import ggplot2
#' @import dplyr
#' @import patchwork
#' @export
plot_selection <- function(
    individs_selected,
    trait_list # a  list of trait data. Does not have to be the same as used in multiopt, in fact, in some cases it is better to not use the same data, as we are more interested in raw trait data than in pairwise data (genomic or otherwise) or scaled variables.
){

  # checks -----------------------------------------------------------

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for this function.",
      call. = FALSE
    )
  }

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop(
      "Package 'dplyr' is required for this function.",
      call. = FALSE
    )
  }

  if (!requireNamespace("patchwork", version = , quietly = TRUE)) {
    stop(
      "Package 'patchwork' is required for this function.",
      call. = FALSE
    )
  }

  stopifnot(utils::packageVersion("patchwork") >= "1.3.0")

  if(any(lapply(trait_list, ncol) != 1)) stop("Trait data should only have 1 column for each trait. If you are attempting to use a pairwise matrix, this is not supported.")

  # loop through traits if needed -------------------------------------------
  dat <- as.data.frame(individs_selected)

  if(length(trait_list) == 1) { # only 1 trait

    # get data table put together
    trait_dat = trait_list[[1]]

    trait_dat = cbind(
      trait_dat,
      colSums(dat > 0) # how many times an individual was selected across replicates (ignores weights)
    )

    colnames(trait_dat) <- c("trait", "n_selected")

    # put plot together
    trait_dat |>
      as.data.frame() |>
      dplyr::mutate(id = dplyr::row_number()) |>
      ggplot2::ggplot(ggplot2::aes(x = reorder(id, -trait),
                                   y = trait,
                                   fill = n_selected)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::theme_classic() +
      ggplot2::scale_fill_gradientn(
        colours = rev(
          grDevices::hcl.colors(
            300,
            palette = "Purple-Yellow"
          )
        )
      ) +
      ggplot2::labs(
        x = "Individual",
        y = names(trait_list),
        fill = "Numer of times\nselected in\nreplicates") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank()
      ) +
      ggplot2::scale_y_continuous(expand = c(0, 0))

  } else {

    # multiple traits

    # format data
    trait_names = names(trait_list)

    trait_dat = as.data.frame(trait_list)

    trait_dat = cbind(
      trait_dat,
      n_selected = rowSums(dat > 0) # how many times an individual was selected across replicates (ignores weights)
    )


    # All pairwise combinations
    trait_pairs <- combn(trait_names, 2, simplify = FALSE)

    plot_list <- lapply(trait_pairs, function(vars){

      trait_dat |>
        ggplot2::ggplot(
          ggplot2::aes(
            x = .data[[vars[1]]],
            y = .data[[vars[2]]],
            fill = n_selected
          )) +
        ggplot2::geom_point(alpha = 0.5, size = 3, shape = 21, color = "gray30") +
        ggplot2::theme_bw() +
        ggplot2::scale_fill_gradientn(
          colours = rev(
            grDevices::hcl.colors(
              300,
              palette = "Purple-Yellow"
            )
          )
        ) +          ggplot2::labs(
          fill = "Numer of times\nselected in\nreplicates"
        )


    }
    )

    patchwork::wrap_plots(plot_list) +
      patchwork::plot_layout(guides = "collect") &
      ggplot2::theme(legend.position = "bottom")

  }

}



#' Visualize trait values of selected individuals
#'
#' Produces plots showing how frequently individuals are selected across
#' optimization replicates in relation to their (and other) trait values.
#' This is similar to `plot_selection` except input values are formatted
#' for `singleopt_context` output.
#'
#' For each individual, selection frequency is computed as the number of
#' replicates in which the individual was selected (ignoring selection
#' weights).
#'
#' @param individs_selected A named list of matrices indicating selected
#'   individuals across replicates (e.g. `out$individs_selected` returned by
#'   `singleopt_context`). Each element represents the results from a different
#'   optimization scenario.
#' @param trait_list Named list of trait data frames, each containing a
#'   single column of trait values. This data does not need to be the same
#'   data supplied to `singleopt_context`; in many cases it is preferable to use
#'   unscaled or raw trait values for visualization. Trait values must be in
#'   the same individual order as `individs_selected`. Pairwise matrices or
#'   multi-column objects are not supported.
#'
#' @return A named list of plots, with one plot per optimization scenario.
#'   For a single trait, each element is a ggplot2 bar plot. For multiple
#'   traits, each element is a patchwork object containing all pairwise
#'   trait scatterplots.
#'
#' @details
#' The function operates in two modes:
#'
#' **Single trait case**
#' \itemize{
#'   \item Combines trait values with selection frequency.
#'   \item Orders individuals by trait value.
#'   \item Produces a bar plot of trait values.
#'   \item Bar fill indicates the number of replicates in which each
#'         individual was selected.
#' }
#'
#' **Multi-trait case**
#' \itemize{
#'   \item Combines all trait values into a single data frame.
#'   \item Computes selection frequency for each individual.
#'   \item Generates scatterplots for all pairwise combinations of traits.
#'   \item Point fill indicates the number of replicates in which each
#'         individual was selected.
#'   \item Pairwise plots are combined into a single figure using
#'         `patchwork`.
#' }
#'
#' Selection frequency is computed as the number of replicates in which an
#' individual has a nonzero selection weight.
#'
#' @import ggplot2
#' @import dplyr
#' @import patchwork
#' @export
plot_selection_single <- function(
    individs_selected,
    trait_list # a  list of trait data. Does not have to be the same as used in multiopt, in fact, in some cases it is better to not use the same data, as we are more interested in raw trait data than in pairwise data (genomic or otherwise) or scaled variables.
){

  # checks -----------------------------------------------------------

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for this function.",
      call. = FALSE
    )
  }

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop(
      "Package 'dplyr' is required for this function.",
      call. = FALSE
    )
  }

  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop(
      "Package 'patchwork' is required for this function.",
      call. = FALSE
    )
  }

  if(any(lapply(trait_list, ncol) != 1)) stop("Trait data should only have 1 column for each trait. If you are attempting to use a pairwise matrix, this is not supported.")

  # loop through traits if needed -------------------------------------------

  if(length(trait_list) == 1) { # only 1 trait

    # get data table put together
    # for each set of weights provided in individs selected, make a table of trait value and weights
    trait_dat = lapply(individs_selected, function(x){

      temp = cbind(
        trait_list[[1]],
        colSums(x > 0) # how many times an individual was selected across replicates (ignores weights)
      )

      colnames(temp) <- c("trait", "n_selected")

      return(temp)
    })


    # put plot together
    plot_out <- lapply(trait_dat, function(x){

      x |>
        as.data.frame() |>
        dplyr::mutate(id = dplyr::row_number()) |>
        ggplot2::ggplot(ggplot2::aes(
          x = stats::reorder(id, -trait),
          y = trait,
          fill = n_selected
        )) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::theme_classic() +
        ggplot2::scale_fill_gradientn(
          colours = rev(
            grDevices::hcl.colors(
              300,
              palette = "Purple-Yellow"
            )
          )
        ) +
        ggplot2::labs(
          x = "Individual",
          y = names(trait_list),
          fill = "Numer of times\nselected in\nreplicates") +
        ggplot2::theme(
          axis.text.x = ggplot2::element_blank(),
          axis.ticks = ggplot2::element_blank()
        ) +
        ggplot2::scale_y_continuous(expand = c(0, 0))

    })


  } else {

    # multiple traits

    # format data
    trait_names = names(trait_list)

    trait_dat = lapply(individs_selected, function(x){

      trait_df = as.data.frame(trait_list)

      temp = cbind(
        trait_df,
        colSums(x > 0) # how many times an individual was selected across replicates (ignores weights)
      )

      colnames(temp) <- c(colnames(trait_df), "n_selected")

      return(temp)

    })

    # All pairwise combinations
    trait_pairs <- combn(trait_names, 2, simplify = FALSE)

    plot_out <- lapply(seq_along(individs_selected), function(x){

      plot_list <- lapply(trait_pairs, function(vars){

        trait_dat[[x]] |>
          ggplot2::ggplot(
            ggplot2::aes(
              x = .data[[vars[[1]]]],
              y = .data[[vars[[2]]]],
              fill = n_selected
            )) +
          ggplot2::geom_point(alpha = 0.5, size = 3, shape = 21, color = "gray30") +
          ggplot2::theme_bw() +
          ggplot2::scale_fill_gradientn(
            colours = rev(
              grDevices::hcl.colors(
                300,
                palette = "Purple-Yellow"
              )
            )
          ) +
          ggplot2::labs(
            fill = "Numer of times\nselected in\nreplicates"
            # title = paste("Optimization for", names(individs_selected)[x])
          )

      })

      patchwork::wrap_plots(plot_list) +
        patchwork::plot_layout(guides = "collect") +
        patchwork::plot_annotation(
          title = paste("Optimization for", names(individs_selected)[x]),
          theme = ggplot2::theme(legend.position = "bottom")
        ) #&
        # ggplot2::theme(legend.position = "bottom")


    })

    names(plot_out) <- names(individs_selected)

  }

  return(plot_out)

}


#' Visualize pairwise Pareto trade-offs
#'
#' Generates a grid of pairwise scatterplots showing trade-offs among
#' traits in a Pareto archive. Optionally overlays results from
#' multi-objective and single-objective optimization runs for comparison.
#'
#' @param archive_list A list containing Pareto archive output,
#'   expected to include an element `archive_summary` with trait values.
#'   Either  $archive from `multiopt_sa` or `rand_multiopt` or output from `explore_pareto`
#' @param multi_list Optional list containing multi-objective optimization
#'   results from `rand_multiopt` or `multiopt_sa`. Must include either `measure_summaries` or
#'   `final_selection$measure_summary`.
#' @param single_list Optional list of single-objective optimization
#'   results from `singleopt_context`.
#' @param measure_type Optional list of direction of optimization.
#'   List elements should have the same names as traits in archive_list and
#'   values should be `"minimize"`, `"maximize"`, or `"diversify"`.
#'
#' @return A patchwork object consisting of multiple ggplot2
#'   pairwise scatterplots with shared legends.
#'
#' @details
#'
#' Required packages: ggplot2, dplyr, patchwork.
#'
#'Note that simulated annealing as used here is trying to maximize values.
#'Therefore, raw output from MultiOpt SA simulations will be maximizing.
#'If you back transform trait values to their original scale then you will
#'be possibly maximizing, minimizing, or diversifying.
#'
#' @import ggplot2
#' @import dplyr
#' @import patchwork
#'
#' @examples
#'
#' # simulate data
#' set.seed(12345)
#' n = 100
#' x = rnorm(n = n, mean = 120, sd = 2)
#' y = x * 3 + rnorm(n = n, mean = 0, sd = 20)
#' dat = data.frame(x = x, y = y)
#'
#' trait_list = list(
#'x = as.matrix(dat$x),
#'y = as.matrix(dat$y)
#')
#'
#'measure_list = list(
#'  x = weighted_mean_of_vector,
#'  y = weighted_mean_of_vector
#')
#'
#'args_list = list(
#'  x = NULL,
#'  y = list(direction = -1)
#')
#'
#'# scale data
#' trait_list_scaled <- scale_traits(trait_list)
#'
#' sa_args = list(
#'trait_list = trait_list_scaled,
#'measure_list = measure_list,
#'measure_args_list = args_list,
#'n_t = 10, # select 10 individuals
#'weights_max = 1, # individuals can only be selected once
#'nda = T
#')
#'
#'test = rand_multiopt(n_runs = 5, multiopt_args = sa_args)
#'
#' # plot the pareto front
#' plot_pareto(archive_list = test$archive)
#'
#'# You can also add multi- and single- objective outputs
#'
#' # single objective simulation
#'single_out <-
#'  singleopt_context(
#'    trait_list = trait_list_scaled,
#'    measure_list = measure_list,
#'    measure_args_list = args_list,
#'    n_t = n_t,
#'    n_runs = 10,
#'    verbose = F
#'  )
#'
#'  plot_pareto(archive_list = test$archive,
#'   multi_list = test,
#'   single_list = single_out
#'  )
#'
#'# and you can additionally add some arrows which helps viewer's remember what
#'direction we are trying to optimize.
#'
#'  plot_pareto(archive_list = test$archive,
#'   multi_list = test,
#'   single_list = single_out,
#'   measure_type = list(x = "maximize", y = "maximize")
#'  )
#'
#'
#' @export
plot_pareto <- function(
    archive_list,
    multi_list = NULL, # optional for multi-objective; rand_multiopt or multiopt_sa
    single_list = NULL, # optional for single-objective output from singleopt_context
    measure_type = NULL # optional list of what direction SA was optimizing
) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for this function.",
      call. = FALSE
    )
  }

  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop(
      "Package 'dplyr' is required for this function.",
      call. = FALSE
    )
  }


  if (is.null(archive_list)) stop("Archive is NULL.")



  # get ready ---------------------------------------------------------------
  archive_dat = as.data.frame(archive_list$archive_summary)

  # get trait info
  trait_names = names(archive_dat)

  # all pairwise combinations
  trait_pairs <- combn(trait_names, 2, simplify = FALSE)

  archive_dat$group <- "archive" # for legend

  ## If multi-objective results are provided ---------------------------------
  if (!is.null(multi_list)) {

    # try and assign
    multi_dat <- multi_list$measure_summaries

    # if it fails, try the other version
    if (is.null(multi_dat)) multi_dat <- multi_list$final_selection$measure_summary

    # if that fails, we have an issue
    if (is.null(multi_dat)) stop("There is an issue with the multi_list formatting.")

    multi_dat = as.data.frame(multi_dat)

    multi_dat$group <- "multi" # for legend

  }


  ## If single-objective outputs are provided --------------------------------
  if (!is.null(single_list)) {

    single_dat <-
      single_list$measure_summaries |>
      dplyr::bind_rows(.id = "trait") |>
      dplyr::as_tibble() |>
      dplyr::mutate(
        group = trait # for legend
      )

  }

  ## if measure information is provided --------------------------------------
  if (!is.null(measure_type)) {

    # measure_type should specify the direction of each trait
    # check inputs

    if(!all(measure_type %in% c("maximize", "minimize", "diversify"))) {
      stop("measure_type list must contain only values of `maximize`, `minimize`, or `diversify`")
    }

    if(!setequal(names(measure_type), colnames(archive_list$archive_summary))) {
      stop("List names for measure_type do not match names in archive_list")
    }

    # set arrow direction. This will be the argument for arrow(type = )
    arrow_directions = lapply(measure_type, function(z) {

      lookup <- c(maximize = "last",
                  minimize = "first",
                  diversify = "both")

      temp <- unname(lookup[z])

      return(temp)
    })

  }


  # make plot -------------------------------------------------------

  plot_list <- lapply(trait_pairs, function(vars){

    p <- archive_dat |>
      ggplot2::ggplot(
        ggplot2::aes(
          x = .data[[vars[1]]],
          y = .data[[vars[2]]],
          color = group
        )
      ) +
      ggplot2::geom_point(
        alpha = 0.5,
        size = 3
      ) +
      ggplot2::theme_bw()

    # add multi-objective if needed
    if (!is.null(multi_list)) {

      p <- p +
        ggplot2::geom_point(data = multi_dat,
                            shape = 3,
                            alpha = 0.75,
                            size = 3,
                            ggplot2::aes(
                              color = group # need to specify in aes() so it appears on the legend
                            )
        )

    }

    # add single-objective if needed
    if (!is.null(single_list)) {

      # make the color and shape palettes
      cols = grDevices::hcl.colors(
        n = length(trait_names),
        palette = "Purple-Yellow"
      )

      shapes = rep(
        c(21, 22, 23, 24, 25), # only want the ones with outlines
        length.out = length(trait_names)
      )

      # add to plot
      p <- p +
        ggplot2::geom_point(data = single_dat,
                            ggplot2::aes(
                              fill = trait,
                              shape = trait
                            ),
                            alpha = 0.75,
                            size = 3,
                            color = "gray40"
        ) +
        ggplot2::scale_fill_manual(values = cols) +
        ggplot2::scale_shape_manual(values = shapes)+
        ggplot2::labs(
          shape = "Single-objective\nresults",
          fill = "Single-objective\nresults",
          color = NULL
        )


    }

    # add objective direction if specified
    if(!is.null(measure_type)) {

      # get information about the plot so far
      p_info <- ggplot2::ggplot_build(p)

      # X-axis limits
      x_lims = p_info$layout$panel_params[[1]]$x$get_limits()

      # Y-axis limits
      y_lims = p_info$layout$panel_params[[1]]$y$get_limits()

      # desired length of arrow as a % of axis length
      arrow_length <- 0.15

      # Build arrow along x-axis
      x_arrow <- data.frame(
        x = x_lims[1] + 0.01 * diff(x_lims),
        xend = x_lims[1] + (0.01 + arrow_length) * diff(x_lims),
        y = y_lims[1],
        yend =  y_lims[1]
      )

      # Build arrow along y-axis
      y_arrow <- data.frame(
        x = x_lims[1],
        xend = x_lims[1],
        y = y_lims[1] + 0.01 * diff(y_lims),
        yend = y_lims[1] + (0.01 + arrow_length) * diff(y_lims)
      )

      # add arrows to plot
      p <- p +
        geom_segment(
          data = x_arrow,
          aes(x = x, xend = xend, y = y, yend = yend),
          inherit.aes = FALSE,
          linewidth = 1,
          lineend = "butt", # this makes the arrow sharp rather than smoothed
          linejoin = "mitre", # this makes the arrow sharp rather than smoothed
          arrow = arrow(
            ends = arrow_directions[[vars[1]]],
            type = "closed",
            angle = 25, # width of arrow
            length = unit(0.2, "cm")
          )
        ) +

        geom_segment(
          data = y_arrow,
          aes(x = x, xend = xend, y = y, yend = yend),
          inherit.aes = FALSE,
          linewidth = 1,
          lineend = "butt",
          linejoin = "mitre",
          arrow = arrow(
            ends = arrow_directions[[vars[2]]],
            type = "closed",
            angle = 25,
            length = unit(0.2, "cm")
          )
        ) +
        coord_cartesian(clip = "off") +
        theme(aspect.ratio = 1) # Forces the plot canvas to be square which helps make the arrows the same length

    }


    # scale front and final results
    p <-  p +
      ggplot2::scale_color_manual(
        values = c(
          "black", # Pareto front
          "hotpink" # multi-objective results
        ),
        labels = c(
          "Pareto front",
          "Multi-objective\nresults"
        )) +
      ggplot2::labs(
        color = NULL
      )

    return(p)
  })

  final_plot <- patchwork::wrap_plots(plot_list) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")

  return(final_plot)

}



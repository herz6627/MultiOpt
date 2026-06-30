#' Run single-objective simulated annealing in context of additional traits
#'
#' This function aims to address "If I prioritize a single trait for simulated annealing,
#' what are the results of the other traits?." This is complimentary analysis
#' to `multiopt_sa` and `rand_multiopt`, allowing users to explore what trade offs
#' happen when incorporating multiple objectives. No chains or archives are returned.
#'
#'
#'Compatible with parallel computing when n_reps > 1. See `rand_multiopt()` for more details.
#'
#'Output value can be back transformed using `unscale_singleopt`.
#'
#' @inheritParams multiopt_sa
#' @inheritParams rand_multiopt
#'
#' @returns A list of the same length of `trait_list` with an element for each
#' trait found in `trait_list`. Within each list element, a table with measure
#' results for each trait is returned, with the column with the same name as the
#' list element being the single objective used in the simulated annealing run(s).
#' @export
#'
#' @examples
#'
#' set.seed(12345)
#' n = 100
#'x = rnorm(n = n, mean = 120, sd = 2)
#'y = x * 3 + rnorm(n = n, mean = 0, sd = 20)
#'dat = data.frame(x = x, y = y)
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
#'  y = list(direction = -1) # to minimize the mean of this trait
#')
#'
#' trait_list_scaled <- scale_traits(trait_list)
#'
#'set.seed(123)
#'singleopt_context(
#'  trait_list = trait_list_scaled,
#'  measure_list = measure_list,
#'  measure_args_list = args_list,
#'  n_t = 10,
#'  verbose = F
#')
#'
#'

singleopt_context <- function(
    trait_list,
    measure_list,
    measure_args_list,
    n_t = NULL,
    weights_min = NULL,
    weights_max = NULL,
    initial_weights = NULL,
    max_steps = 10000,
    max_t = 1,
    p_depends_delta = F,
    acceptance_multipliers = 1,
    verbose = T,
    n_runs = 1,
    parallel = F
){

  # checks
  if(n_runs <= 0 || n_runs %% 1 != 0) stop("`n_runs` must be positive and an integer.")

  # run multi_opt for a single trait ----------------------------------------

  if (n_runs == 1) {

    # single run --------------------------------------------------------------

    # loop over each list element for a single-objective simulated annealing
    sim_out <- Map(
      f = function(trait, measure, measure_args, nm) {

        multiopt_sa(
          trait_list = stats::setNames(list(trait), nm),
          measure_list = stats::setNames(list(measure), nm),
          measure_args_list = stats::setNames(list(measure_args), nm),
          n_t = n_t,
          weights_max = weights_max,
          weights_min = weights_min,
          initial_weights = initial_weights,
          max_t = max_t,
          max_steps = max_steps,
          p_depends_delta = p_depends_delta,
          acceptance_multipliers = acceptance_multipliers,
          nda = F, # cant have an nda with only 1 trait
          save_chain = F, # not storing chain information
          verbose = verbose
        )

      },
      trait_list,
      measure_list,
      measure_args_list,
      names(trait_list) # helps with re-naming the list
    )


    ## get measures for non-target trait ---------------------------------------

    # loop over every trait summary and condense to a table
    sum_out = lapply(
      sim_out,
      function(x){
        z = unlist(calculate_measure(
          trait_list,
          measure_list,
          list_of_args = measure_args_list,
          w = x$final_selection$individs_selected
        ))

        z = data.frame(t(as.matrix(z)))
      })

    # quick check that everything makes sense
    # first column of the tables we just made should match the original measure_summaries output (for the first trait)
    if(sum_out[[1]][1] != sim_out[[1]]$final_selection$measure_summary[1]) stop("Something went wrong. Outputs do not match.")

  } else {

    # multiple runs -----------------------------------------------------------

    # loop over each list element for a single-objective simulated annealing
    sim_out <- Map(

      f = function(trait, measure, measure_args, nm) {

        rand_args <- list(
          trait_list = stats::setNames(list(trait), nm),
          measure_list = stats::setNames(list(measure), nm),
          measure_args_list = stats::setNames(list(measure_args), nm),
          n_t = n_t,
          weights_max = weights_max,
          weights_min = weights_min,
          initial_weights = initial_weights,
          max_t = max_t,
          max_steps = max_steps,
          p_depends_delta = p_depends_delta,
          acceptance_multipliers = acceptance_multipliers,
          nda = FALSE,
          save_chain = FALSE,
          verbose = verbose
        )

        rand_multiopt(
          n_runs = n_runs,
          parallel = parallel,
          multiopt_args = rand_args,
          verbose = verbose
        )

      },

      trait_list,
      measure_list,
      measure_args_list,
      names(trait_list) # helps with re-naming the list

    )

    ## get measures for non-target trait ---------------------------------------

    # use the weights to calculate the measures for the other traits (individs_selected)
    # summarise each row
    calc_one_row <- function(x, i, trait_list, measure_list, measure_args_list){

      w <- x$individs_selected[i, ]

      res <- calculate_measure(
        trait_list,
        measure_list,
        list_of_args = measure_args_list,
        w = w
      )

      unlist(res)
    }

    # iterate over each trait (element in list)
    calc_one_sim <- function(x, trait_list, measure_list, measure_args_list){

      n <- nrow(x$individs_selected)

      rows <- lapply(
        seq_len(n),
        calc_one_row,
        x = x,
        trait_list = trait_list,
        measure_list = measure_list,
        measure_args_list = measure_args_list
      )

      as.data.frame(do.call(rbind, rows))
    }

    # iterate over each row and element
    sum_out <- lapply(
      sim_out,
      calc_one_sim,
      trait_list = trait_list,
      measure_list = measure_list,
      measure_args_list = measure_args_list
    )

    # quick check that everything makes sense
    # first column of the tables we just made should match the original measure_summaries output (for the first trait)
    if(!all(sum_out[[1]][,1] == sim_out[[1]]$measure_summaries[,1])) stop("Something went wrong. Outputs do not match.")

  }

  return(sum_out)
}



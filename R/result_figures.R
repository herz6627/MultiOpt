



sa_plots_multi = function(sa_loop_out, max_steps = NULL, trait1_name, trait2_name) {

  archive_dat = as_tibble(sa_loop_out$archive$archive_values)
  trait_dat = sa_loop_out$trait_means
  chain_dat = sa_loop_out$sample_chain
  individ_dat = sa_loop_out$individs

  # pareto front and final results for all replicates
  p1 = ggplot() +
    geom_point(data = archive_dat,
               aes(x = value_1, y = value_2), size = 3) +                         # pareto front
    geom_point(data =trait_dat,                                                   # final result from each replicate
               aes(x = v1, y = v2), color = "red", alpha = 0.5) +
    labs(
      x = trait1_name,
      y  = trait2_name,
      title = "Pareto front") +
    theme_bw()

  # chain from one of the runs, with the final result in red
  p2 = ggplot() +
    geom_point(data = chain_dat, # one of the runs
               aes(x = value_1, y = value_2), color = "gray", alpha = .2) +
    geom_point(data = chain_dat[max_steps,], # final result
               aes(x = value_1, y = value_2), color = "red", alpha = .5) +
    labs(
      x = trait1_name,
      y = trait2_name,
      title = "SA sample") +
    theme_bw()

  # chain and pareto together
  p3 = ggplot() +

    # chain
    geom_point(data = chain_dat, # one of the runs
               aes(x = value_1, y = value_2), color = "gray", alpha = .1) +

    # pareto front
    geom_point(data = archive_dat,
               aes(x = value_1, y = value_2), size = 3) +                         # pareto front
    geom_point(data =trait_dat,                                                   # final result from each replicate
               aes(x = v1, y = v2), color = "red", alpha = 0.5) +

    labs(
      x = trait1_name,
      y  = trait2_name,
      title = "Pareto front and sample run") +
    theme_bw()

  # number of times an individual was selected
  p4 = individ_dat %>%
    ggplot(aes(x = .data[[trait1_name]], y = .data[[trait2_name]], fill = value)) +
    geom_point(alpha = 0.5, size = 3, shape = 21, color = "gray30") +
    theme_bw() +
    paletteer::scale_fill_paletteer_c("grDevices::Purple-Yellow", direction = -1)+
    labs(
      # x = "Stem length",
      # y = "Leaf length",
      fill = "Numer of times\nselected in\nreplicates")

  return(list(parato = p1, chain = p2, pareto_chain = p3, individ_selec = p4))

}


sa_plots_single = function(sa_loop_out, max_steps) {

  traits = names(sa_loop_out)

  for (i in seq(length(traits))) {

    trait_dat = sa_loop_out[[traits[[i]]]]$trait_means
    chain_dat = sa_loop_out[[traits[[i]]]]$sample_chain
    individ_dat = sa_loop_out[[traits[[i]]]]$individs


    # final results for all replicates
    p1 = ggplot(data = trait_dat, aes(x = v1)) +                                                   # final result from each replicate
      geom_histogram(bins = 50) +
      labs(
        x = traits[[i]],
        title = "Run results") +
      theme_bw()

    # chain from one of the runs, with the final result in red
    p2 = chain_dat %>%
      mutate(t = row_number(),
      ) %>%
      ggplot(aes(y = value, x = t)) +
      geom_line() +
      labs(
        x = "Step",
        y =  traits[[i]],
        title = "SA sample") +
      theme_bw()


    # number of times an individual was selected
    p3 = individ_dat %>%
      ggplot(aes(x = reorder(bolt_id, -.data[[traits[[i]]]]), y = .data[[traits[[i]]]], fill = value)) +
      geom_bar(stat = "identity") +
      theme_classic() +
      paletteer::scale_fill_paletteer_c("grDevices::Purple-Yellow", direction = -1)+
      labs(
        x = "Individual",
        y = traits[[i]],
        fill = "Numer of times\nselected in\nreplicates") +
      theme(
        axis.text.x = element_blank(),
        axis.ticks = element_blank()
      ) +
      scale_y_continuous(limits = c(0,1), expand = c(0, 0))
  }

  return(list(results = p1, chain = p2, individ_selec = p3))

}

sa_plots = function(
    sa_multi_loop_out = NULL,
    sa_single_loop_out = NULL,
    max_steps,
    trait_names,
    return_all = T, # if T and both multi and single results are supplied, all possible figures are returned. If F, returns only figures using both datasets.
    show_multi = T # If T, all multi objective results will be shown. If F, just pareto front and single objective results will be shown.
) {

  out = list()

  # multi objective ------------------------------------------------
  if (!is.null(sa_multi_loop_out) & return_all) {

    out <- append(out, sa_plots_multi(sa_multi_loop_out, trait1_name = trait_names[[1]], trait2_name = trait_names[[2]], max_steps = max_steps))

  } # end of multi

  # single objective ------------------------------------------------
  if (!is.null(sa_single_loop_out) & return_all) {

    out <- append(out, sa_plots_single(sa_single_loop_out, max_steps))

  } # end of single


  # both single and multi --------------------------------------------
  if (!is.null(sa_single_loop_out) & !is.null(sa_multi_loop_out)) {

    if(any(!trait_names == names(sa_single_loop_out))) stop("`trait_names` do not match names in `sa_single_loop_out`.")
    if (is.null(sa_single_loop_out[[trait_names[[1]]]]$sup_trait_means) | is.null(sa_single_loop_out[[trait_names[[2]]]]$sup_trait_means)) {
      stop("`sup_trait_means` is missing from `sa_single_loop_out`")
    }

    # prep single-objective
    single_dat1 = sa_single_loop_out[[trait_names[[1]]]]$sup_trait_means %>%
      filter(trait == trait_names[[2]]) %>%
      rename(v2 = trait.mean) %>%
      mutate(v1 = pull(sa_single_loop_out[[trait_names[[1]]]]$trait_means))

    single_dat2 = sa_single_loop_out[[trait_names[[2]]]]$sup_trait_means %>%
      filter(trait == trait_names[[1]]) %>%
      rename(v1 = trait.mean) %>%
      mutate(v2 = pull(sa_single_loop_out[[trait_names[[2]]]]$trait_means))


    # pareto front and final results for all replicates
    lab1 <- paste("Single objective:\n", trait_names[[1]])
    lab2 <- paste("Single objective:\n", trait_names[[2]])

    p = ggplot() +

      # pareto front
      geom_point(data = as_tibble(sa_multi_loop_out$archive$archive_values),
                 aes(x = value_1, y = value_2, color = "Pareto front"), size = 5) +

      # final result to trait 1 single
      geom_point(data =  single_dat1,
                 aes(x = v1, y = v2, color = lab1), alpha = 0.75, size = 4) +

      # final result to trait 2 single
      geom_point(data =  single_dat2,
                 aes(x = v1, y = v2, color = lab2), alpha = 0.75, size = 4) +
      scale_color_manual(
        values = setNames(
          c("black", "#5495CFFF", "#DB4743FF", "#7C873EFF"),
          c("Pareto front", "Multi-objective", lab1, lab2)
        ),
        name = "Legend"
      ) +
      labs(
        x = trait_names[[1]],
        y  = trait_names[[2]],
        title = "Pareto front and all results") +
      theme_bw()


    if(show_multi) {
      p = p +
        # final result from each multi replicate
        geom_point(data = sa_multi_loop_out$trait_means,
                   aes(x = v1, y = v2, color = "Multi-objective"), alpha = 0.75, size = 4)
    } else p = p

    out <- append(out, p)
  }

  return(out)
}

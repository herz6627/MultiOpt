#' Apply measure functions across multiple trait datasets
#'
#' Evaluates a set of user-supplied measure functions across multiple
#' trait datasets using a common vector of individual weights.
#' Each trait dataset is paired with a corresponding measure function
#' and optional additional arguments.
#'
#' The function iterates over the named elements of
#' `list_of_trait_data`, applies the corresponding function from
#' `list_of_measures`, and returns the resulting values as a named list.
#'
#' Measure functions must accept the trait data through an argument
#' named `v`, and weights through an argument named `w`.
#'
#' @param list_of_trait_data Named list of trait datasets. Each element
#'   is passed to its corresponding measure function as argument `v`. Trait data must be matrices.
#' @param list_of_measures Named list of functions corresponding to
#'   `list_of_trait_data`. Each function is applied to the matching
#'   trait dataset.
#' @param list_of_args Optional named list of additional argument lists
#'   supplied to each measure function. Names should match
#'   `list_of_trait_data`.
#' @param w Numeric vector of weights applied across all measure functions.
#'
#' @return A named list containing the output of each measure function.
#'
#' @details
#' For each trait:
#' \enumerate{
#'   \item The corresponding trait dataset is extracted.
#'   \item The associated measure function is retrieved and validated.
#'   \item Additional user-supplied arguments are combined with the
#'         common weight vector `w`.
#'   \item The function is evaluated using `do.call()`.
#' }
#'
#' The names of `list_of_trait_data`, `list_of_measures`, and
#' `list_of_args` are expected to align.
#'
#' @export
calculate_measure <- function(
    list_of_trait_data,
    list_of_measures,
    list_of_args = list(),
    w
){

  out <- vector("list", length(list_of_trait_data))
  names(out) <- names(list_of_trait_data)

  for (trait in names(list_of_trait_data)) {

    # get the trait data and what function will be used
    x <- list_of_trait_data[[trait]]
    fun <- list_of_measures[[trait]]

    if (is.null(fun)) {
      stop("No measure supplied for: ", trait)
    }

    if (!is.function(fun)) {
      stop("Measure for ", trait, " must be a function")
    }

    measure_args <- list_of_args[[trait]]

    if (is.null(measure_args)) measure_args <- list()

    # weights are an argument for most measure functions
    measure_args[["w"]] <- w

    out[[trait]] <- do.call(fun, c(list(v = x), measure_args))

  }

  out
}


#' Nei's Genetic Diversity
#'
#' Calculates Nei's genetic diversity (expected heterozygosity) across loci
#' in a population. Allele frequencies are calculated from a genotype matrix,
#' with the option to incorporate individual-level weights.
#'
#' @param v A genotype matrix with individuals in rows and loci in columns.
#' Genotypes should be coded as the number of copies of an allele (e.g.,
#' 0, 1, or 2).
#' @param w An optional numeric vector of individual weights. If supplied,
#' weighted allele frequencies are calculated.
#' @param chunk_size Numeric value for how large (number of SNPs) to chunk large
#' genotype datasets by. If value is larger than ncol(v), data is not chunked.
#'
#' @return A numeric value representing the mean Nei's genetic diversity
#' (expected heterozygosity) across loci.
#'
#' @details
#' For each locus, the frequency of the focal allele is calculated as the
#' proportion of allele copies in the population. Nei's gene diversity is
#' then calculated as the expected heterozygosity.
#'
#' For each locus, the allele frequency is calculated as:
#'
#' \deqn{p_i = \frac{\sum_j g_{ji}}{2n}}
#'
#' where \eqn{g_{ji}} is the genotype of individual \eqn{j} at locus \eqn{i},
#' and \eqn{n} is the number of individuals.
#'
#' When weights are provided, the weighted allele frequency is:
#'
#' \deqn{p_i = \frac{\sum_j w_j g_{ji}}{2\sum_j w_j}}
#'
#' Nei's gene diversity for each locus is calculated as:
#'
#' \deqn{H_i = 1 - p_i^2 - (1-p_i)^2 = 2p_i(1-p_i)}
#'
#' The final value is the mean diversity across all \eqn{L} loci:
#'
#' \deqn{H = \frac{1}{L}\sum_{i=1}^{L} H_i}
#'
#' Assumes loci are bi-allelic.
#'
#' For extra large matrices, this function accepts BEDmatrix formats
#' (BEDMatrix package). There are probably other options that will work. If the
#' genotype matrix object runs with \code{colSums()}, it will probably work here.
#'
#' @examples
#' # Create a small example genotype matrix.
#' # Rows represent individuals and columns represent loci.
#' # Genotypes are coded as the number of copies of the focal allele (0, 1, or 2).
#' geno <- matrix(
#'   c(0, 1, 2,
#'     1, 1, 0,
#'     2, 0, 1,
#'     1, 2, 1),
#'    nrow = 4,
#'    byrow = TRUE
#')
#'
#' # Calculate Nei's gene diversity.
#' nei_diversity(geno)
#'
#'
#' # Calculate weighted Nei's gene diversity.
#' weights <- c(1, 0.5, 1.5, 1)
#' nei_diversity(geno, w = weights)
#'
#' @references
#' Nei, M. (1973). Analysis of gene diversity in subdivided populations. Proceedings of the National Academy of Sciences of the United States of America, 70, 3321–3323.
#'
#' Based on code by Jason Bragg in the OptGenMix package (jasongbragg/OptGenMix)
#' as used in these publications:
#'
#' Bragg, J. G., Cuneo, P., Sherieff, A., & Rossetto, M. (2020). Optimizing the genetic composition of a translocation population: Incorporating constraints and conflicting objectives. Molecular Ecology Resources, 20(1), 54–65. https://doi.org/10.1111/1755-0998.13074
#'
#' Bragg, J. G., Yap, J.-Y. S., Wilson, T., Lee, E., & Rossetto, M. (2021). Conserving the genetic diversity of condemned populations: Optimizing collections and translocation. Evolutionary Applications, 14(5), 1225–1238. https://doi.org/10.1111/eva.13192
#'
#' @export
nei_diversity <- function(v,
                          w = NULL,
                          chunk_size = 5000
){

  # get info
  n_snps <- ncol(v) # number of snps
  n_ind <- nrow(v) # number of individuals with snp data

  if (!is.null(w) && n_ind != length(w)) stop("\nnrow(v) does not equal length(w).")

  # --- calculate frequency of alleles across population (p) ---

  total <- 0
  # to allow for large matrices without making huge vectors, we are chunking
  # up the matrix and slowly adding up the p values, rather than making a large
  # p vector which we would then sum.

  for (start in seq(1, n_snps, by = chunk_size)) {

    end <- min(start + chunk_size - 1, n_snps)

    v_chunk = as.matrix(v[, start:end])

    if (is.null(w)) {

      p <- colSums(v_chunk) / (2 * n_ind)

    } else {

      p <- colSums(v_chunk * w) / (2 * sum(w)) # since our weights indicate if an individual is selected more than once, we use the weights, not number of individuals here.
    }

    total <- total + sum(2 * p * (1 - p))

  }

  # --- calculate Nei ---
  nei <- total / n_snps

  return(nei)
}

#' Shannon diversity
#'
#' @param v
#'   Numeric genotype matrix with individuals in rows and loci in columns.
#'   Entries are assumed to be allele counts (0, 1, 2).
#' @param w
#'   Optional numeric vector of individual weights (length must equal
#'   nrow(v)). If NULL, all individuals are treated equally.
#' @param q
#'   Diversity order:
#'   \describe{
#'     \item{0}{Locus polymorphism indicator (monomorphic = 1, polymorphic = 2)}
#'     \item{1}{Shannon diversity (Hill number q = 1)}
#'     \item{2}{Simpson diversity (Hill number q = 2)}
#'   }
#'  @param direction numeric scalar. Multiplier applied to the final metric value
#'   to control orientation. Use 1 for default direction, -1 to invert the sign.
#'
#' @details
#' Direction is applied after metric computation and does not alter the
#' underlying metric definition.
#'
#' @returns
#' A single numeric value representing the mean per-locus Hill diversity
#' across all loci.
#'
#' @note
#' I have kept the function true to the version found in OptGenMix, but have clarified what 'q' is actually doing.
#' The original description in OptGenMix was sparse, so I have had to do some extrapolation.
#'
#' @export
shannon_diversity <- function(v, w=NULL, q=1, direction = 1) {

  # run checks
  if (!is.matrix(v)) stop("v must be a matrix.")

  if(!direction %in% c(-1, 1)) stop("`direction` must be -1 or 1")


  # calculate allele frequencies

  if (is.null(w)) {

    p  <- colSums(v) / (nrow(v) *2)

  } else {

    wm <- matrix(rep( w, ncol(v)), ncol=ncol(v), byrow = FALSE)
    p  <- colSums(v*wm) / (nrow(v)*2) / mean(w)

  }

  # q=0, allele count
  if (q == 0) {

    a   <- p
    a[ which(p == 1) ] <- 1
    a[ which(p == 0) ] <- 1
    a[ which(p > 0 & p < 1) ] <- 2

    D_0_mean <- mean(a)

    return(direction * D_0_mean)
  }

  # q=1, Shannon diversity
  if (q == 1) {

    H_1 <- rep(0,ncol(v))

    for (i in 1:ncol(v)) {

      p_i <- p[i]
      q_i <- 1 - p_i

      if (p_i == 0 | p_i == 1) {
        # 0 times log 0 has limit of 0
        # but R will return NaN
        H_1[i] <- - log(1)

      } else {
        H_1[i] <- -1 * ( p_i * log(p_i) + q_i * log(q_i) )
      }

    }

    D_1 <- exp(H_1)

    D_1_mean <- mean(D_1)

    return(direction * D_1_mean)

  }

  # q=2, Simpson diversity
  if (q == 2) {

    H_2 <- rep(0,ncol(v))

    for (i in 1:ncol(v)) {

      p_i <- p[i]
      q_i <- 1 - p_i

      H_2[i] <- 1 - ( p_i*p_i + q_i*q_i) # Nei diversity

    }

    D_2 <- 1 / ( 1 - H_2 )

    D_2_mean <- mean(D_2)

    return(direction * D_2_mean)
  }
}

#' Desirable allele enrichment
#'
#' Calculates enrichment of desirable allele. Given genotypes (v),
#' a vector (w) of weight values (equal length to the
#' number of individuals, a vector (v) of values describing
#' the importance of each locus, estimates an index of enrichment of
#' preferred alleles. Assumes genotype fitness 2 > 1 > 0.
#'
#' @param v Genotype matrix (individuals × loci), coded as allele dosage (0, 1, 2).
#' @param w Optional numeric vector of individual weights (length = nrow(v)).
#'   If NULL, all individuals are weighted equally.
#' @param loc Optional numeric vector of locus weights (length = ncol(v)).
#'   If NULL, loci are weighted equally.
#' @param rec Logical. If TRUE, heterozygotes (1) are treated as 0,
#'   enforcing a recessive model where only homozygotes for the allele contribute.
#'  @param direction numeric scalar. Multiplier applied to the final metric value
#'   to control orientation. Use 1 for default direction, -1 to invert the sign.
#'
#' @details
#' Direction is applied after metric computation and does not alter the
#' underlying metric definition.
#'
#' @return A single numeric value representing the weighted allele enrichment index.
#'   Higher values indicate greater enrichment of the allele across individuals and loci.
#' @export
allele_enrichment <- function(v, w = NULL, loc = NULL, rec = FALSE, direction = 1) {

  # run checks
  if (!is.matrix(v)) stop("v must be a matrix.")

  if(!direction %in% c(-1, 1)) stop("`direction` must be -1 or 1")

  # calculate
  if (rec) {
    v[ v == 1 ] <- 0
  }

  if (is.null(w)) {
    p  <- colSums(v) / (nrow(v) *2)
  } else {
    wm <- matrix(rep(w,ncol(v)),ncol=ncol(v),byrow=FALSE)
    p  <- colSums(v*wm) / (nrow(v)*2) / mean(w)
  }

  if (is.null(loc)) {
    a  <- sum(p) / length(p)
  } else {
    a  <- sum(p*loc) / sum(loc)
  }

  return(direction * a)
}

#' Compute weighted deviation of a vector from a target value
#'
#' Calculates the absolute difference between a weighted mean of a numeric vector (`vs`)
#' and a specified target value (`disp`). Weights are applied element-wise and
#' normalized by their sum. This could be used to minimize the use of individuals
#' low levels of heterozygosity, or to place arbitrary
#' constraints on genotype composition.
#'
#' @param v single-column matrix of numeric values. Will be coerced to a vector.
#' @param w Numeric vector of individual weights with same length as `v`.
#' @param disp Numeric scalar target value to compare the weighted mean against. Defaults to 0.
#' @param direction numeric scalar. Multiplier applied to the final metric value
#'   to control orientation for simulated annealing. Use 1 for default direction
#'   (maximize trait value), -1 to invert the sign (minimize trait value).
#'
#' @return A single numeric value representing the absolute deviation between
#'         the weighted mean of `v` and `disp`.
#'
#' @details The function computes:
#' \deqn{ | ( \sum v_i w_i / \sum w_i ) - disp | }
#'
#' Direction is applied after metric computation and does not alter the
#' underlying metric definition.
#' @examples
#' v <- as.matrix(c(1, 2, 3))
#' w <- c(0.2, 0.3, 0.5)
#' weighted_mean_of_vector(v, w, disp = 2)
#'
#' @export
weighted_mean_of_vector <- function(v, w, disp = 0, direction = 1) {

  # run checks
  if (!is.matrix(v)) stop("v must be a matrix.")

  if (dim(v)[2] != 1) stop("v should have only 1 column.")

  if(!direction %in% c(-1, 1)) stop("`direction` must be -1 or 1")

  # force into vector
  vs_vec <- as.numeric(v)

  # check lengths
  if (length(vs_vec) != length(w)) {
    stop("`v` and `w` must have the same length.")
  }

  ss <- abs((sum(v*w) / sum(w)) - disp)

  return(direction * ss)

}


#' #' Compute weighted sum of squared deviations from a target value
#'
#' Calculates the weighted sum of squared differences between a numeric vector
#' (`v`) and a displacement value (`disp`). Each squared deviation is weighted
#' by a corresponding value in `w`.This could be used, e.g., to minimize
#' the mean difference between temperature of origin for each sample,
#' and temperate of a site.
#'
#' @inheritParams weighted_mean_of_vector
#' @return A single numeric value representing the weighted sum of squared
#'         deviations from `disp`.
#' @details The function computes:
#' \deqn{ \sum_i w_i (v_i - disp)^2 }
#'
#' @details
#' Direction is applied after metric computation and does not alter the
#' underlying metric definition.
#'
#' @examples
#' v <- as.matrix(c(1, 2, 3))
#' w <- c(0.2, 0.3, 0.5)
#' sum_of_squared_difference(v, w, disp = 2)
#'
#' @export
sum_of_squared_difference <- function(v, w, disp=0, direction = 1) {

  # run checks
  if (!is.matrix(v)) stop("v must be a matrix.")

  if (dim(v)[2] != 1) stop("v should have only 1 column.")

  if(!direction %in% c(-1, 1)) stop("`direction` must be -1 or 1")

  # force into vector
  vs_vec <- as.numeric(v)

  # check lengths
  if (length(vs_vec) != length(w)) {
    stop("`v` and `w` must have the same length.")
  }

  ss <- sum(((vs_vec - disp)^2) * w)

  return(direction * ss)
}

#' Compute weighted mean absolute deviation from a target value
#'
#' Calculates the weighted mean of absolute differences between a numeric vector
#' (`v`) and a displacement value (`disp`). Each absolute deviation is weighted
#' by a corresponding value in `w`, and the result is normalized by the sum of
#' weights.
#'
#' This can be used as a loss function to quantify average absolute departure
#' from a reference value, for example deviation of sampled environmental values
#' from a target site condition.
#'
#' @inheritParams weighted_mean_of_vector
#'
#' @return A single numeric value representing the weighted mean absolute
#'         deviation of `v` from `disp`.
#'
#' @details The function computes:
#' \deqn{ \frac{\sum_i w_i |v_i - disp|}{\sum_i w_i} }
#'
#' Direction is applied after metric computation and does not alter the
#' underlying metric definition.
#'
#' @examples
#' v <- as.matrix(c(1, 2, 3))
#' w <- c(0.2, 0.3, 0.5)
#' weighted_mean_of_absolute_difference(v, w, disp = 2)
#'
#' @export
weighted_mean_of_absolute_difference <- function(v, w, disp=0, direction = 1) {

  # run checks
  if (!is.matrix(v)) stop("v must be a matrix.")

  if (dim(v)[2] != 1) stop("v should have only 1 column.")

  if(!direction %in% c(-1, 1)) stop("`direction` must be -1 or 1")

  # force into vector
  vs_vec <- as.numeric(v)

  # check lengths
  if (length(vs_vec) != length(w)) {
    stop("`v` and `w` must have the same length.")
  }

  ss <- sum(abs(vs_vec - disp) * w) / sum(w)

  return(direction * ss)
}


#' Compute weighted mean pairwise similarity (or distance)
#'
#' Calculates the mean value of a pairwise similarity (or distance) matrix,
#' applying weights to individuals. The function supports two cases:
#'
#' @param v Square pairwise matrix of similarities or distances between individuals.
#'        Must have dimensions n × n.
#' @param w Numeric vector of weights (length must match nrow(v)).
#' @param direction Numeric scalar. Multiplier applied to the final metric value
#'   to control orientation. Use 1 for default direction, -1 to invert the sign.
#'
#' @return A single numeric value representing the weighted mean pairwise value.
#'         Returns NA if fewer than two weighted individuals remain.
#'
#' @details
#' For binary weights (w == 1 or 0), the function reduces the matrix to selected
#' individuals and computes (to save time):
#' \deqn{ mean(sm_{ij}) \; for \; i < j }
#'
#' For general weights, it computes:
#' \deqn{
#' \frac{\sum_{i \le j} w_i w_j sm_{ij}}{\sum_i w_i (\sum_i w_i - 1)/2}
#' }
#' with a correction term for diagonal contributions due to repeated sampling.
#'
#' Direction is applied after metric computation and does not alter the
#' underlying metric definition.
#' @examples
#' v <- matrix(runif(16), 4, 4)
#' diag(v) <- 0
#' w <- c(1, 2, 0, 1)
#' weighted_mean_of_pairwise_matrix(v, w)
#'
#' @export
weighted_mean_of_pairwise_matrix <- function(v, w, direction = 1) {

  if(!is.matrix(v) | dim(v)[1] != dim(v)[2]) stop("v must be a pairwise matrix.")

  if(!direction %in% c(-1, 1)) stop("`direction` must be -1 or 1")

  wlen <- length(w)
  wsum <- sum(w)

  if (length(w) != nrow(v)) {stop("Length of w must match number of rows in v.")}

  # if weights are 0 or 1, the calculation is simple (and fast)
  if (all(w %in% 0:1)) {

    sub_sm <- v[w == 1, w == 1] # select only weighted individuals
    return(direction * mean(sub_sm[upper.tri(sub_sm, diag = F)])) # dont need to wory about the diagonal when an individual is only selected max once

  } else{

    # if weights are more complicated we need to adjust methods (slower)

    # subset matrix by removing unselected individuals (speeds up the function)
    idx <- which(w > 0)

    if (length(idx) < 2) return(NA_real_)

    sub_sm <- v[idx, idx]
    sub_w  <- w[idx]

    W <- outer(sub_w, sub_w) # get matrix of weights
    weighted_mat <- sub_sm * W # multiply matrix by weights

    diag_term <- sum(diag(sub_sm) * (sub_w - 1) * sub_w / 2) # diagonal contribution (repeated sampling of an individual)

    off_diag_term <- sum(weighted_mat[upper.tri(weighted_mat)]) # off-diagonal (i < j)

    numerator <- off_diag_term + diag_term # sum of all trait values

    denominator <- sum(sub_w) * (sum(sub_w) - 1) / 2 # number of all pairwise combinations

    return(direction * (numerator / denominator)) # get mean

  }
}


#' Calculate the Wasserstein Distance for use in simulated annealing
#'
#' Uses the \code{wasserstein()} to calculate the Wasserstein distance (\eqn{W_p}) between
#' the original trait data distribution and the proposed subset of individuals.
#' Assumes that v provides the total population distribution of trait variation.
#'
#' @param v single-column matrix of numeric values. Will be coerced to a vector.
#' @param w Numeric vector of individual weights with same length as `v`.
#' @param ... Optional additional arguments. See \code{?wasserstein} for specific
#' options.
#'
#' @details
#' As `MultiOpt` aims to maximize measures, but smaller \eqn{W_p} are
#' better, we multiply the \eqn{W_p} by -1 to allow for minimizing
#' the differences in distributions.
#'
#'
#' This function examines how well distributions match. If you are more interested
#' in checking the full range of trait data is captured in the subset,
#' \code{trait_coverage()} would be a better fit.
#'
#'
#' @returns A single numeric value representing the p-Wasserstein distance multiplied
#' by -1.
#' @export
#'
#' @examples
#' v <- as.matrix(c(1, 2, 3))
#' w <- c(0, 1, 1)
#' wasserstein_measure(v, w)
wasserstein_measure <- function(v, w, ...) {

  # run checks
  if (!is.matrix(v)) stop("v must be a matrix.")

  if (dim(v)[2] != 1) stop("v should have only 1 column.")

  # if(!direction %in% c(-1, 1)) stop("`direction` must be -1 or 1")

  if (any(w < 0) || any(w != floor(w))) {
    stop("`w` must contain non-negative integer selection counts.")
  }

  # force into vector
  vs_vec <- as.numeric(v)

  # check lengths
  if (length(vs_vec) != length(w)) {
    stop("`v` and `w` must have the same length.")
  }

  # replicate each observed value w times (weights of 0 will get dropped).
  v_sub <- rep(vs_vec, w)

  # Get W
  wass <- wasserstein(
    observations_a = vs_vec,
    observations_b = v_sub,
    ...
  )

  # since we are trying to maximize values in SA but smaller Wasserstein values
  # are better, we'll just multiply the distance value by -1
  return(wass * -1)

}


#' Calculate Trait Coverage
#'
#' Calculates the proportion of trait-value bins represented in a selected
#' subset relative to the bins represented in the full dataset.
#'
#' The range of trait values in v is divided into n_bins equally spaced
#' bins. The function counts how many bins contain observations in the full
#' dataset and how many contain observations in the selected subset, defined
#' by the selection counts in w. The returned value is the proportion of
#' occupied bins in the subset relative to the full dataset.
#'
#' @param v single-column matrix of numeric values. Will be coerced to a vector.
#' @param w Numeric vector of individual weights with same length as `v`.
#' to the observations in v. Values with a count of zero are excluded
#' from the selected subset.
#' @param n_bins The number of equally spaced bins used to divide the trait
#' range.
#'
#' @return A numeric value representing the proportion of trait-value bins
#' covered by the selected subset. A value of 1 indicates that the subset
#' covers all bins represented in the full dataset.
#'
#' @details
#' This function does not measure how well the proposed subset matches the
#' distribution of trait data. For this type of measurement,
#' \code{wasserstein_measure()} would be a better fit.
#'
#' @examples
#' v <- as.matrix(c(1, 2, 3))
#' w <- c(0, 1, 2)
#'
#' trait_coverage(v, w, n_bins = 10)
trait_coverage <- function(v, w, n_bins = 10) {

  # run checks
  if (!is.matrix(v)) stop("v must be a matrix.")

  if (dim(v)[2] != 1) stop("v should have only 1 column.")

  if (any(w < 0) || any(w != floor(w))) {
    stop("`w` must contain non-negative integer selection counts.")
  }

  if(n_bins > length(v)) warning("n_bins > length(v), check that this makes sense.")

  if(n_bins > sum(w)) warning("n_bins > sum(w), check that this makes sense,
                              as output will be artifically low due to forced blank bins.")

  if (length(v) != length(w)) stop("length(v) != length(w)")

  # force into vector
  v_vec <- as.numeric(v)

  # for proposed solution:
  # replicate each observed value w times (weights of 0 will get dropped).
  v_sub <- rep(v_vec, w)

  # define breakpoints using the full dataset
  breakpoints <- seq(

    from = min(v_vec),          # should be 0 if data is scaled
    to = max(v_vec),            # should be 1 if data is scaled
    length.out = n_bins + 1     # length.out will be how many edges we have
    # (so breakpoints-1 is how many bins we will have)
  )

  # categorize full dataset into intervals
  intervals_all <- table(cut(v_vec, breaks = breakpoints, include.lowest = TRUE)) # table() makes the frequency table
  # barplot(intervals)

  freqs_all <- as.data.frame(intervals_all) # format to df

  # now categorize for the subset data
  intervals_sub <- table(cut(v_sub, breaks = breakpoints, include.lowest = TRUE))

  freqs_sub <- as.data.frame(intervals_sub)

  # total bins filled in the full dataset (if we have some outliers we could have a few blanks)
  filled_bins_all <- sum(freqs_all$Freq > 0)
  filled_bins_sub <- sum(freqs_sub$Freq > 0)

  # proportion bins filled.
  # making a proportion fits in with the scales of the SA algorithm
  # and also adjusts for missing bins in the _all dataset
  # (and is more interpretable across bin sizes)
  prop_filled <- filled_bins_sub / filled_bins_all

  if(prop_filled > 1) warning("Proportion is >1. This may be an issue.")

  return(prop_filled)
}



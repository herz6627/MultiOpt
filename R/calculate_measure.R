



calculate_measure <- function(
    list_of_trait_data,
    list_of_measures,
    list_of_args = list(), # arguments needed for measure function
    w # weight vector
){

  out <- vector("list", length(list_of_trait_data))
  names(out) <- names(list_of_trait_data)

  for (trait in names(list_of_trait_data)) {

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


#' Nei diversity
#'
#' Calculates Nei genetic diversity from genotype data.
#'
#' @param v Genotype matrix (rows = individuals, columns = loci)
#' @param w Optional vector of weights
#' @param direction numeric scalar. Must be 1 or -1. Applied as a multiplicative
#'   factor to the computed metric to control its optimisation direction.
#'
#' @details
#' Direction is applied after metric computation and does not alter the
#' underlying metric definition.
#'
#' @returns Numeric scalar of calculated Nei diversity.
#' @note
#' Nei, M. (1973). Analysis of gene diversity in subdivided populations. Proceedings of the National Academy of Sciences of the United States of America, 70, 3321–3323.
#' @export
nei_diversity <- function(v, w=NULL, direction = 1) {

  # run checks
  if (!is.matrix(v)) stop("v must be a matrix.")
  if (nrow(v) <= 1 || ncol(v) <= 1) stop("v should have more than 1 row/col.")
  if(!direction %in% c(-1, 1)) stop("`direction` must be -1 or 1")

  # estimate allele frequencies
  if (is.null(w)) {

    p  <- colSums(v) / (nrow(v) *2)

  } else {

    if (nrow(v) == length(w)) {

      wm <- matrix(rep(w, ncol(v)), ncol = ncol(v), byrow=FALSE)
      p  <- colSums(v*wm) / (nrow(v)*2) / mean(w)

    } else {

      stop("weights supplied: must have length equal to number of rows in v.")

    }
  }

  # calculate diversity from frequencies
  Hs <- 0

  for (i in 1:ncol(v)) {

    Hs <- Hs + ( 1-(p[i])^2 - (1-p[i])^2  )

  }

  # finish and return
  nei <- 1 / ncol(v) * Hs

  return(direction * nei)
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
  if (nrow(v) <= 1 || ncol(v) <= 1) stop("v should have more than 1 row/col.")
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
  if (nrow(v) <= 1 || ncol(v) <= 1) stop("v should have more than 1 row/col.")
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
#'   to control orientation. Use 1 for default direction, -1 to invert the sign.
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
  if (nrow(v) <= 1) stop("v should have only 1 row.")
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
  if (nrow(v) <= 1) stop("v should have only 1 row.")
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
  if (nrow(v) <= 1) stop("v should have only 1 row.")
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

  if(!is.matrix(v) | nrow(v) != ncol(v)) stop("v must be a pairwise matrix.")
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

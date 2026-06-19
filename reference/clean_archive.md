# Reconcile a candidate solution with a non-dominated (Pareto) archive

Updates a Pareto archive by comparing a new mixture solution against
existing archived solutions across multiple objectives. The function
removes any archived solutions dominated by the new mixture and rejects
the candidate if it is dominated by any existing archive member.

## Usage

``` r
clean_archive(summary, wts, archive)
```

## Arguments

- summary:

  Named list or numeric vector of objective values for the candidate
  solution. If provided as a list, it will be coerced to a numeric
  vector via \`unlist()\`.

- wts:

  Numeric vector decision variables associated with the candidate
  solution.

- archive:

  List of existing non-dominated solutions. Each element is a list
  containing:

  archive_summary

  :   Numeric matrix of objective values, where each row corresponds to
      a solution and each column to an objective.

  archive_weights

  :   Numeric matrix of decision variables, where each row corresponds
      to a solution and each column an individual.

## Value

Updated archive list containing only non-dominated solutions after
reconciliation with the candidate.

## Details

Dominance is defined in the Pareto sense: a solution A dominates B if A
is greater than or equal to B in all objectives and strictly greater in
at least one objective.

Exact duplicates (within objective space) are treated as redundant and
are not added to the archive.

The function performs three operations:

1.  Checks whether the candidate is dominated by any archived solution.

2.  Identifies and removes archived solutions dominated by the
    candidate.

3.  Inserts the candidate if it is not dominated.

## Examples

``` r
archive <- list(
  archive_summary = matrix(c(1, 2, 2, 1), ncol = 2, byrow = TRUE),
  archive_weights = matrix(c(0.1, 0.2, 0.3, 0.4), ncol = 2, byrow = TRUE)
)

candidate_summary <- list(x = 2, y = 2)
candidate_weights <- c(0.5, 0.6)

archive <- reconcile_sample_nondominated_archive(
  summary = candidate_summary,
  wts = candidate_weights,
  archive = archive
)
#> Error in reconcile_sample_nondominated_archive(summary = candidate_summary,     wts = candidate_weights, archive = archive): could not find function "reconcile_sample_nondominated_archive"
```

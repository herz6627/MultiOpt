# Pareto Shape Index (PSI)

The Pareto Shape Index (PSI; \\\gamma\\) describes the shape of the
Pareto front (ranging from convex, linear, or concave) for a
two-objective comparison of traits. The PSI quantifies the geometric
shape of a Pareto front by averaging the vertical distances of pairwise
solution-line intersections from the center of the normalized objective
space.

## Usage

``` r
pareto_shape_index(pareto_front)
```

## Arguments

- pareto_front:

  \`archive_summary\` output (non-dominated archive) from a simulated
  annealing run where rows are non-dominated archive values and a column
  for each trait. Only two traits at a time are possible.

## Value

Single numeric value of the Pareto Shape Index (PSI).

## Details

This index ranges from -1 to 1, where -1 indicates a concave Pareto
front (weak trade off), 0 indicates a linear front (strong trade off),
and 1 indicates a convex front (which is not desirable when we are
maximizing our objectives as this is a very harsh trade off). When there
is no trade off between objectives, then the Pareto front is formed by a
single solution; therefore, the value of Ncr in Equation (6) in Unal et
al (2017) becomes zero, leading the Pareto Shape Index to have a value
of infinity.

Trait measures must be in the same direction (minimizing or maximizing).
This should be automatically accounted for when using the MultiOpt
measure functions as MultiOpt simulated annealing attempts to maximize
objectives (measure functions will add a negative for minimizations).
Thus it is important to not do any transformations to the Pareto outputs
before calculating PSI.

The Pareto Shape Index is calculated following Equations (4)–(6) of Unal
et al. (2017). For each pair of Pareto-optimal solutions \\i\\ and
\\j\\, the intersection of the corresponding solution lines between
objectives \\k\\ and \\l\\ is calculated as

\$\$ Y\_{ij,kl} = \frac{ (S\_{i,k}S\_{j,l}) - (S\_{i,l}S\_{j,k}) }{
S\_{i,k} - S\_{i,l} - S\_{j,k} + S\_{j,l} } \$\$

where \\S\_{i,k}\\ and \\S\_{i,l}\\ are the normalized objective values
for solution \\i\\. Objectives should be normalized to the interval \[0,
1\] and transformed to a common optimization direction prior to
computing the index.

The signed distance of each intersection from the center of the
normalized objective space is then

\$\$ P\_{ij,kl} = Y\_{ij,kl} - 0.5 \$\$

Finally, the Pareto Shape Index is calculated as the mean signed
distance across all pairwise solution intersections:

\$\$ \gamma\_{kl} = \frac{1}{N\_{cr}} \sum\_{i=1}^{N_p-1}
\sum\_{j=i+1}^{N_p} P\_{ij,kl} \$\$

where \\N_p\\ is the number of Pareto-optimal solutions and \\N\_{cr} =
N_p(N_p-1)/2\\ is the number of pairwise solution intersections. The
resulting index ranges from -1 to 1, with negative values indicating
concave Pareto fronts, values near zero indicating approximately linear
tradeoffs, and positive values indicating convex Pareto fronts.

## References

Unal, M., Warn, G. P., & Simpson, T. W. (2017). Quantifying the Shape of
Pareto Fronts During Multi-Objective Trade Space Exploration. Journal of
Mechanical Design, 140(021402). https://doi.org/10.1115/1.4038005

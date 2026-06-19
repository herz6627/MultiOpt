# Determines acceptance or rejection of a proposal

Determines acceptance or rejection of a proposal

## Usage

``` r
accept_reject(
  summary,
  proposal_summary,
  t,
  p_depends_delta = FALSE,
  c = 1,
  c_all = 1
)
```

## Arguments

- summary:

  List of current objective values.

- proposal_summary:

  List of proposed objective values.

- t:

  Current temperature.

- p_depends_delta:

  If TRUE, acceptance depends on magnitude of worsening.

- c:

  Penalty multipliers. Single value of vector of same length as
  \`summary\`.

- c_all:

  Multiplier when all objectives worsen.

## Value

Logical scalar indicating whether proposal is accepted.

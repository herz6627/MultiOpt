# Apply measure functions across multiple trait datasets

Evaluates a set of user-supplied measure functions across multiple trait
datasets using a common vector of individual weights. Each trait dataset
is paired with a corresponding measure function and optional additional
arguments.

## Usage

``` r
calculate_measure(
  list_of_trait_data,
  list_of_measures,
  list_of_args = list(),
  w
)
```

## Arguments

- list_of_trait_data:

  Named list of trait datasets. Each element is passed to its
  corresponding measure function as argument \`v\`. Trait data must be
  matrices.

- list_of_measures:

  Named list of functions corresponding to \`list_of_trait_data\`. Each
  function is applied to the matching trait dataset.

- list_of_args:

  Optional named list of additional argument lists supplied to each
  measure function. Names should match \`list_of_trait_data\`.

- w:

  Numeric vector of weights applied across all measure functions.

## Value

A named list containing the output of each measure function.

## Details

The function iterates over the named elements of \`list_of_trait_data\`,
applies the corresponding function from \`list_of_measures\`, and
returns the resulting values as a named list.

Measure functions must accept the trait data through an argument named
\`v\`, and weights through an argument named \`w\`.

For each trait:

1.  The corresponding trait dataset is extracted.

2.  The associated measure function is retrieved and validated.

3.  Additional user-supplied arguments are combined with the common
    weight vector \`w\`.

4.  The function is evaluated using \`do.call()\`.

The names of \`list_of_trait_data\`, \`list_of_measures\`, and
\`list_of_args\` are expected to align.

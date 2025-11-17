# Generate markdown reports

Generate a markdown report for outcomes from a single operating model
projection

## Usage

``` r
# S4 method for class 'SMSE'
report(
  object,
  name = object@Name,
  filename = "SMSE",
  dir = tempdir(),
  open_file = TRUE,
  render_args = list(),
  ...
)
```

## Arguments

- object:

  [SMSE](https://docs.salmonmse.com/reference/SMSE-class.md) object

- name:

  Character string for the model name to include in the report, e.g.,
  model run number.

- filename:

  Character string for the name of the markdown and HTML files.

- dir:

  The directory in which the markdown and HTML files will be saved.

- open_file:

  Logical, whether the HTML document is opened after it is rendered.

- render_args:

  List of arguments to pass to
  [`rmarkdown::render()`](https://pkgs.rstudio.com/rmarkdown/reference/render.html).

- ...:

  Additional arguments (not used)

## Value

Returns invisibly the output of
[`rmarkdown::render()`](https://pkgs.rstudio.com/rmarkdown/reference/render.html),
typically the path of the output file

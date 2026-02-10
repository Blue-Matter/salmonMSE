# Compare scenarios in markdown

Generate a markdown report for multiple model runs to compare scenarios.
3-5 is likely the ideal number of scenarios for comparison.

## Usage

``` r
compare(
  SMSE_list,
  names,
  col_vec,
  filename = "SMSEcompare",
  dir = tempdir(),
  open_file = TRUE,
  render_args = list(),
  ...
)
```

## Arguments

- SMSE_list:

  List of [SMSE](https://docs.salmonmse.com/reference/SMSE-class.md)
  objects

- names:

  Character vector `length(SMSE_list)` to label individual model runs

- col_vec:

  Character vector `length(SMSE_list)` for custom colour schemes for
  comparing across model scenarios in figures

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

## See also

[`report()`](https://docs.salmonmse.com/reference/report.md)

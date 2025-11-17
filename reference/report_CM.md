# Conditioning model markdown report

Generate a markdown report to plot time series and MCMC posteriors of
estimates from the conditioning model. See
[`get_report()`](https://docs.salmonmse.com/reference/CMfigures.md) for
the various plotting functions used in the report.

## Usage

``` r
report_CM(
  stanfit,
  year,
  cov1_names,
  cov_names,
  rs_names,
  name,
  filename = "CM",
  dir = tempdir(),
  open_file = TRUE,
  render_args = list(),
  ...
)
```

## Arguments

- stanfit:

  Output from
  [`sample_CM()`](https://docs.salmonmse.com/reference/fit_CM.md)

- year:

  Optional vector of calendar years

- cov1_names:

  Optional character vector for names of covariates that predict age-1
  natural mortality

- cov_names:

  Optional character vector for names of covariates that predict age-2+
  natural mortality

- rs_names:

  Optional character vector for names of hatchery release strategies

- name:

  Optional character string for the model name to include in the report,
  e.g., model run number

- filename:

  Character string for the name of the markdown and HTML files

- dir:

  The directory in which the markdown and HTML files will be saved.

- open_file:

  Logical, whether the HTML document is opened after it is rendered

- render_args:

  List of arguments to pass to
  [`rmarkdown::render()`](https://pkgs.rstudio.com/rmarkdown/reference/render.html)

- ...:

  Additional arguments (not used)

## Value

Returns invisibly the output of
[`rmarkdown::render()`](https://pkgs.rstudio.com/rmarkdown/reference/render.html),
typically the path of the output file

## Details

Report excludes MCMC values from warmup iterations

## See also

[`fit_CM()`](https://docs.salmonmse.com/reference/fit_CM.md)
[`get_report()`](https://docs.salmonmse.com/reference/CMfigures.md)

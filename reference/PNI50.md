# Example performance metrics

Functions that evaluate return probabilities of outcomes from the
simulations.

## Usage

``` r
PNI50(SMSE, Ref = 0.5, Yrs = NULL)

PNI80(SMSE, Ref = 0.8, Yrs = NULL)

WILD50(SMSE, Ref = 0.5, Yrs = NULL)

SMSY85(SMSE, Ref = 0.85, Yrs = NULL)

Sgen100(SMSE, Ref = 1, Yrs = NULL)
```

## Arguments

- SMSE:

  SMSE object returned by
  [`salmonMSE()`](https://docs.salmonmse.com/reference/salmonMSE.md)

- Ref:

  Threshold for the performance metric, used to calculate the
  probability that the metric exceeds this value

- Yrs:

  Numeric vector of length 2 to indicate the year range over which to
  summarize performance. If NULL, the performance is summarized over all
  projection years.

## Value

A vector of probabilities corresponding to population

## Details

- `PNI50` calculates the probability that PNI exceeds 0.50 (threshold
  for an integrated-transition population, Withler et al. 2018)

- `PNI80` calculates the probability that PNI exceeds 0.80 (threshold
  for an integrated-wild population, Withler et al. 2018)

- `WILD50` calculates the probability that at least 50 percent of
  natural spawners are wild

- `SMSY85` calculates the probability that NOS/SMSY exceeds 0.85

- `Sgen100` calculates the probability that NOS/Sgen exceeds 1

## References

Withler et al. 2018. Genetically Based Targets for Enhanced
Contributions to Canadian Pacific Chinook Salmon Populations. DFO Can.
Sci. Advis. Sec. Res. Doc. 2018/019. xii + 88 p.

# Reference points

## Introduction

Reference points for salmon typically use the Ricker stock function that
predict the return as a density-dependent relationship with the number
of spawners, typically of the previous generation for salmon with single
brood year returns.

The Ricker function is

$$R = a \times S\exp( - bS)$$

where $R$ is the return and $S$ is the spawners. Parameters $a$ is in
units of returns/spawner and $b$ is in units of reciprocal spawners.

Let $S_{\text{MSY}}$ be the spawner abundance that maximizes catch
$R - S$, which can be obtained by maximizing the function:

$$R - S = a \times S\exp( - bS) - S$$

The implicit solution, where $\alpha = \log(a)$, satisfies

$$\left( 1 - bS_{\text{MSY}} \right)\exp\left( \alpha - bS_{\text{MSY}} \right) = 1$$

and can be solved numerically ([Scheuerell
2016](https://doi.org/10.7717/peerj.1623)).

$S_{\text{MSY}}$ is also the spawner abundance that maximizes excess
returns above the one-to-one line.

The harvest rate corresponding to $S_{\text{MSY}}$ is

$$U_{\text{MSY}} = bS_{\text{MSY}}$$

The quantity $S_{\text{gen}}$ is the spawner abundance that would reach
the spawner abundance at MSY after one generation without fishing ([Holt
et
al. 2009](https://www.dfo-mpo.gc.ca/csas-sccs/Publications/ResDocs-DocRech/2009/2009_058-eng.htm)),
and satisfies the equation:

$$S_{\text{MSY}} = a \times S_{\text{gen}}\exp\left( - b \times S_{\text{gen}} \right)$$

There are several implicit assumptions frequently associated with these
reference points:

1.  The stock-recruit relationship predicts the number of returns, not
    juveniles, from the number of spawners
2.  All returns are the same age
3.  All exploitation occurs on mature individuals
4.  All returns are equally vulnerable to the terminal fishery
5.  All spawners are equally fecund

## Reference points in salmonMSE

salmonMSE is an age-structured model that can relax the assumptions from
the typical spawner-return model:

1.  The stock-recruit relationship predicts the number of outmigrating
    juveniles, not the return, from the egg production (not the number
    of spawners)
2.  The calendar-year return can comprise of multiple brood years
3.  Exploitation can occur on both immature and mature individuals
4.  Fish vulnerability can vary by age
5.  Older, larger spawners have larger contributions to the egg
    production because they are more fecund than younger fish

Accordingly, reference points are calculated with more structural
complexity compared to the simple approach. However, comparisons below
indicate the many situations where the salmonMSE reference points
overlap with the typical salmon approach.

### Calculations for MSY

Given exploitation in both the immature and mature component of the
population, the juvenile survival of an outmigrating individual
(“smolt”) to age $a$ is

$$\ell_{a} = \prod\limits_{j = 1}^{a - 1}\exp\left( - \left\lbrack M_{j} + v_{j}^{\text{PT}}F^{\text{PT}} \right\rbrack \right)\left( 1 - p_{j} \right)$$
where $M$ is natural mortality, $v^{\text{PT}}$ is vulnerability at age
in the preterminal fishery, $F^{\text{PT}}$ is fishing mortality in the
preterminal fishery, and $p$ is maturity.

The corresponding return per juvenile
$r_{a} = \ell_{a}\exp\left( - v_{a}^{\text{PT}}F^{\text{PT}}\rbrack \right)\left( 1 - p_{a} \right)$.

The corresponding escapement per juvenile is
$s_{a} = r_{a}\exp\left( - v_{a}^{\text{T}}F^{\text{T}}\rbrack \right)$,
with $v^{\text{T}}$ is vulnerability at age in the fishery and
$F^{\text{T}}$ is fishing mortality in the terminal fishery.

The corresponding egg production per juvenile
$\phi = \sum_{a}s_{a}f_{a}p^{\text{female}}$.

The equilibrium juvenile production $J$ is calculated from $\phi$ as:

$J = \log(a\prime\phi)/(b\prime\phi)$

where $a\prime$ is in units of juvenile per egg and $b\prime$ is in
units of per egg.

The total return is $\sum_{a}J \times r_{a}$ and total spawners is
$\sum_{a}J \times s_{a}$.

The maximum sustainable yield state variables are values obtained from
maximizing the total catch:

$$C = \sum\limits_{a}J \times \ell_{a} \times \left( 1 - \exp\left( v_{a}^{\text{PT}}F^{\text{PT}} \right) \right) + \sum\limits_{a}J \times r_{a} \times \left( 1 - \exp\left( v_{a}^{\text{T}}F^{\text{T}} \right) \right)$$

The optimization calculates the corresponding fishery effort
$E_{\text{MSY}}$. The user needs to specify the relative effort to the
preterminal and terminal fisheries ($e_{1}$ and $e_{2}$), with
$F^{\text{PT,MSY}} = e_{1} \times E_{\text{MSY}}$ and
$F^{\text{T,MSY}} = e_{2} \times E_{\text{MSY}}$.

An alternative optimization scheme maximizes excess recruitment
($R - S$) but should be equivalent to maximizing terminal fishery catch
(no preterminal fishery catch):

$$R - S = \sum\limits_{a}J \times r_{a} - \sum\limits_{a}J \times s_{a}$$

### Calculations for Sgen

Sgen is calculated by a subsequent optimization algorithm to search for
the fishery effort $E_{\text{gen}}$ that generates an equilibrium
distribution at age, then projects that equilibrium population for
$n_{\text{gen}}$ years with no fishing, such that the spawning abundance
reaches $S_{\text{MSY}}$ after $n_{\text{gen}}$ years.

If there is partial maturity at several age classes, the number of years
for the projection is unclear. The population can recover to
$S_{\text{MSY}}$ with either increased survival, increased juvenile
production, or both. Currently, salmonMSE will use the first age of
maturity, e.g., typically two years for Chinook salmon, for the
projection duration. Comparisons in the next section indicate that Sgen
between the two methods sometimes, but not always, agree.

## Comparison between methods

Currently, salmonMSE parameterizes the operating model with $a$
(returns/spawner), $\phi$ (juveniles/egg) and $Smax = 1/b\prime$ (eggs).

The following function facilitates comparison of the salmonMSE and
Ricker SRR reference points in the various examples. Some conversion is
also needed to convert between $b\prime$ (1/eggs) to $b$ (1/spawner).

In summary, if the age-structured reference point calculation only
includes terminal exploitation, then the reference points are equivalent
with values from the Ricker model unless both fecundity and
vulnerability vary with age.

``` r
compare_ref <- function(a = 3, # Units of recruits/spawner
                        b = 1/1000, # Units of reciprocal spawners
                        maxage = 5,
                        rel_F = c(0, 1), # e_1 and e_2
                        p_female = 1,
                        vulPT = rep(0, maxage),
                        vulT = c(0, 0, 0.2, 0.5, 1),
                        M = c(1, 0.8, 0.6, 0.4, 0.2),
                        fec = c(0, 1500, 3000, 3200, 3500),
                        p_mature = c(0, 0.1, 0.2, 0.3, 1)) {
  
  require(salmonMSE)

  # Calculate eggs per smolt
  phi <- salmonMSE:::calc_phi(
    Mjuv = M,
    p_mature = p_mature,
    p_female = p_female,
    fec = fec,
    s_enroute = 1
  )

  # To convert Smax from units of spawners to eggs:
  # 1. Calculate unfished spawners from Ricker SRR, set R = S --> S0 = log(a)/b
  # 2. Calculate spawners per juvenile (spro)
  # 3. Calculate smolt per egg (phi)
  # 4. unfished eggs (eo) = unfished spawners * juvenile per spawner * egg per smolt
  # 5. Smax_eggs = log(a) / unfished eggs
  spro <- salmonMSE:::calc_phi(
    Mjuv = M,
    p_mature = p_mature,
    p_female = p_female,
    fec = fec,
    s_enroute = 1,
    output = "spawner"
  )
  
  Smax <- 1/b
  so <- Smax * log(a)
  eo <- so / spro * phi
  beta_eggs <- log(a)/eo
  Smax_eggs <- 1/beta_eggs

  SRRpars <- data.frame(
    kappa = a,
    Smax = Smax_eggs,
    phi = phi,
    SRrel = "Ricker"
  )

  # salmonMSE ref
  ref <- calc_MSY(
    M, fec, p_female, rel_F, vulPT, vulT, p_mature,
    s_enroute = 1, SRRpars = SRRpars
  )

  ref["Sgen"] <- calc_Sgen(
    M, fec, p_female, rel_F, vulPT, vulT, p_mature,
    s_enroute = 1, SRRpars = SRRpars, SMSY = ref["Spawners_MSY"]
  )
  ref["Catch/Return"] <- ref["KT_MSY"]/ref["Return_MSY"]

  ref_salmonMSE <- structure(
    ref[c("UPT_MSY", "UT_MSY", "Catch/Return", "Spawners_MSY", "Sgen")], 
    names = c("UMSY (preterminal)", "UMSY (terminal)", "Terminal Catch/Return", "SMSY", "Sgen")
  )
  
  # Ricker SRR calcs
  ref_Ricker <- local({
    umsy <- calc_Umsy_Ricker(log(a))
    smsy <- calc_Smsy_Ricker(log(a), b)
    sgen <- calc_Sgen_Ricker(log(a), b)
    structure(c(umsy, smsy, sgen), names = c("UMSY", "SMSY", "Sgen"))
  })
  
  output <- list(
    salmonMSE = ref_salmonMSE,
    Ricker = ref_Ricker
  )
  return(output)
}
```

### Situation 1 - single age of maturity

Situation 1 specifies the salmonMSE (age-structured) model with the
simple assumptions identified in the introduction. The salmonMSE and
Ricker (stock-recruit only) reference points are identical.

``` r
compare_ref(
  a = 3,
  b = 1/1000,
  maxage = 5,
  rel_F = c(0, 1),
  p_female = 1,
  vulPT = rep(0, 5),
  vulT = c(0, 0, 0, 0, 1),
  M = c(-log(0.01), 0, 0, 0, 0), # SAR = 0.01
  fec = c(0, 0, 0, 0, 1),
  p_mature = c(0, 0, 0, 0, 1)
)
#> Loading required package: salmonMSE
#> $salmonMSE
#>    UMSY (preterminal)       UMSY (terminal) Terminal Catch/Return 
#>             0.0000000             0.4678228             0.4678228 
#>                  SMSY                  Sgen 
#>           467.8334508           188.2447444 
#> 
#> $Ricker
#>        UMSY        SMSY        Sgen 
#>   0.4678265 467.8265256 188.2417444
```

### Situation 2 - partial age maturity, equal vulnerability

Situation 2 specifies the salmonMSE (age-structured) model with
age-specific juvenile mortality and maturity. All spawners are identical
with fecundity = 1 and fully vulnerable to the fishery.

The salmonMSE and Ricker-only reference points are identical.

``` r
compare_ref(
  a = 3,
  b = 1/1000,
  maxage = 5,
  rel_F = c(0, 1),
  p_female = 1,
  vulPT = rep(0, 5),
  vulT = rep(1, 5),
  M = c(1, 0.3, 0.2, 0.1, 0.1),
  fec = rep(1, 5),
  p_mature = c(0, 0.1, 0.2, 0.3, 1)
)
#> $salmonMSE
#>    UMSY (preterminal)       UMSY (terminal) Terminal Catch/Return 
#>             0.0000000             0.4678228             0.4678228 
#>                  SMSY                  Sgen 
#>           467.8334508           188.2447444 
#> 
#> $Ricker
#>        UMSY        SMSY        Sgen 
#>   0.4678265 467.8265256 188.2417444
```

### Situation 3 - partial age maturity, varying vulnerability by age

Situation 3 is the same as 2 except with age-varying fishery
vulnerability.

The salmonMSE and Ricker-only reference points are identical.

One nuance is how to report the exploitation rate, since it varies by
age due to differential vulnerability. UMSY reported in salmonMSE is the
maximum value experienced by an age class (typically the oldest), but
the ratio of total catch to total return may be more appropriate metric.

``` r
compare_ref(
  a = 3,
  b = 1/1000,
  maxage = 5,
  rel_F = c(0, 1),
  p_female = 1,
  vulPT = rep(0, 5),
  vulT = c(0, 0.1, 0.2, 0.4, 1),
  M = c(1, 0.3, 0.2, 0.1, 0.1),
  fec = rep(1, 5),
  p_mature = c(0, 0.1, 0.2, 0.3, 1)
)
#> $salmonMSE
#>    UMSY (preterminal)       UMSY (terminal) Terminal Catch/Return 
#>             0.0000000             0.7330604             0.4678273 
#>                  SMSY                  Sgen 
#>           467.8250871           188.2409734 
#> 
#> $Ricker
#>        UMSY        SMSY        Sgen 
#>   0.4678265 467.8265256 188.2417444
```

### Situation 4 - partial age maturity, varying fecundity and vulnerability by age

Situation 4 is the same as 3 except with age-varying fecundity.

The salmonMSE and Ricker-only reference points are **not** identical.

Preferential vulnerability by the fishery for return which are not
identical in fecundity are additional complexities not implemented in
Ricker-only reference points.

``` r
compare_ref(
  a = 3,
  b = 1/1000,
  maxage = 5,
  rel_F = c(0, 1),
  p_female = 1,
  vulPT = rep(0, 5),
  vulT = c(0, 0.1, 0.2, 0.4, 1),
  M = c(1, 0.3, 0.2, 0.1, 0.1),
  fec = c(0, 1000, 2000, 3000, 3500),
  p_mature = c(0, 0.1, 0.2, 0.3, 1)
)
#> $salmonMSE
#>    UMSY (preterminal)       UMSY (terminal) Terminal Catch/Return 
#>             0.0000000             0.6424027             0.4000761 
#>                  SMSY                  Sgen 
#>           523.7657546           260.2641560 
#> 
#> $Ricker
#>        UMSY        SMSY        Sgen 
#>   0.4678265 467.8265256 188.2417444
```

### Situation 5 - partial age maturity, with immature exploitation

Reference points can also be calculated by considering preterminal
exploitation, e.g., [Hankin and Healey
1986](https://doi.org/10.1139/f86-219), although this is not common
practice.

Preterminal mortality is not explicitly modeled in Ricker-only reference
points.

Thus, salmonMSE and Ricker-only reference points are **not** identical.

``` r
compare_ref(
  a = 3,
  b = 1/1000,
  maxage = 5,
  rel_F = c(1, 1),
  p_female = 1,
  vulPT = c(0, 0.1, 0.2, 0.3, 1),
  vulT = c(0, 0.1, 0.2, 0.4, 1),
  M = c(1, 0.3, 0.2, 0.1, 0.1),
  fec = c(0, 1000, 2000, 3000, 3500),
  p_mature = c(0, 0.1, 0.2, 0.3, 1)
)
#> $salmonMSE
#>    UMSY (preterminal)       UMSY (terminal) Terminal Catch/Return 
#>             0.4670406             0.4670406             0.2553111 
#>                  SMSY                  Sgen 
#>           535.7260998           315.6413128 
#> 
#> $Ricker
#>        UMSY        SMSY        Sgen 
#>   0.4678265 467.8265256 188.2417444
```

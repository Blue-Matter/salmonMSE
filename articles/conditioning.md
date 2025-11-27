# Conditioning model from CWT and escapement

## Introduction

Typical run reconstruction uses a historical time series of spawner and
recruit data to estimate the productivity of salmon populations. For
many small populations of interest, escapement time series may be
available but catch may not be identifiable to the local scale, for
example, if the population is a conservation unit (CU) as part of a
larger stock complex and catch composition is not identified to
individual CUs.

[Walters and Korman
(2024)](https://publications.gc.ca/site/eng/9.940685/publication.html)
demonstrated an approach for reconstruction if Coded Wire Tag (CWT) data
from an indicator hatchery are assumed to be representative to life
cycle parameters of a natural population.

The model consists of two components. First, CWT data informs natural
mortality, maturation, and exploitation rate in the marine environment.
Second, these parameters are then applied to the population of interest
and informs the size and productivity of the population that has a time
series of total escapement. The number of hatchery releases in the
system informs hatchery production, and the difference in total
escapement and hatchery production informs natural production. Both
steps are accomplished within a single model fit, which can account for
uncertainty among the various data components, and posterior
distributions of parameters are obtained by MCMC.

We utilize this approach as a conditioning model to inform stochastic
parameters for projections in salmonMSE, although use of this model is
not necessary to set up an operating model.

### Model fitting

Model fitting is performed in RTMB with
[`fit_CM()`](https://docs.salmonmse.com/reference/fit_CM.md).

The posterior can then be sampled with
[`sample_CM()`](https://docs.salmonmse.com/reference/fit_CM.md) through
rstan.

A markdown report is available with
[`report_CM()`](https://docs.salmonmse.com/reference/report_CM.md)

A subset of posterior MCMC draws of parameters to reconstruct the
historical population can be imported into the operating model with
[`CM2SOM()`](https://docs.salmonmse.com/reference/CM2SOM.md). Further
modification of the operating model and additional settings can then be
added to run the projection.

## Variable definitions

| Name                 | Definition                                         |
|:---------------------|:---------------------------------------------------|
| $C$                  | Catch                                              |
| $\text{CWT}$         | Coded wire tag                                     |
| $E$                  | Escapement                                         |
| $F$                  | Instantaneous fishing mortality                    |
| $\text{HO}$          | Hatchery origin                                    |
| $\text{HOS}$         | Hatchery origin spawner                            |
| $M$                  | Instantaneous natural mortality                    |
| $N^{\text{juv}}$     | Juvenile abundance                                 |
| $\text{NO}$          | Natural origin                                     |
| $\text{NOS}$         | Natural origin spawner                             |
| $\text{PT}$          | Preterminal fishery                                |
| $R$                  | Recruitment                                        |
| $\text{T}$           | Terminal fishery                                   |
| $m$                  | Proportion mature                                  |
| $p^{\text{female}}$  | Proportion female                                  |
| $p^{\text{spawn}}$   | Proportion spawning                                |
| $s^{\text{enroute}}$ | En-route survival (escapement to spawning grounds) |
| $v$                  | Fishery vulnerability                              |

## Life cycle parameters from CWT

Life cycle parameters are informed by CWT. Mortality rates are
parameterized in instantaneous units, which can be converted to survival
terms.

Fishing mortality $F$ is separated into a preterminal (PT) component
which acts on juvenile fish and a terminal (T) component on mature fish.
Separable effects are modeled where the fishing mortality is
year-specific and modified by age class by a vulnerability term $v$.

Survival of juvenile CWT to the next age class $a$ at the beginning of
year $y$ is calculated after exploitation, maturation $m$, and natural
mortality $M$:

$$N_{y + 1,a + 1}^{\text{juv,CWT}} = N_{y,a}^{\text{juv,CWT}}\exp\left( - v_{a}^{\text{PT}}F_{y}^{\text{PT}} \right)\left( 1 - m_{y,a} \right)\exp\left( - M_{y,a} \right)$$

An additional mortality term, associated with release mortality, can be
applied for survival in the first year of life shortly after release:

$$N_{y,a = 1}^{\text{juv,CWT}} = N_{y}^{\text{rel,CWT}}s^{\text{rel}}$$

where $s^{\text{rel}}$ is the survival after release into the wild.

The CWT return $R$ is the fraction of maturing juveniles after
preterminal exploitation

$$R_{y,a}^{\text{CWT}} = N_{y,a}^{\text{juv,CWT}}\exp\left( - v_{a}^{\text{PT}}F_{y}^{\text{PT}} \right)m_{y,a}$$

The escapement is the return that survive terminal exploitation

$$E_{y,a}^{\text{CWT}} = R_{y,a}^{\text{CWT}}\exp\left( - v_{a}^{\text{T}}F_{y}^{\text{T}} \right)$$

The catch $C$ is

$$\begin{aligned}
C_{y,a}^{\text{CWT,PT}} & {= N_{y,a}^{\text{juv,CWT}}\left( 1 - \exp\left( - v_{a}^{\text{PT}}F_{y}^{\text{PT}} \right) \right)} \\
C_{y,a}^{\text{CWT,T}} & {= R_{y,a}^{\text{CWT}}\left( 1 - \exp\left( - v_{a}^{\text{T}}F_{y}^{\text{T}} \right) \right)}
\end{aligned}$$

## Population reconstruction

The model does separate accounting of natural-origin ($\text{NO}$) and
hatchery-origin ($\text{HO}$) fish in the population of interest.

### Juvenile life stage

The age-1 HO fish is obtained from the number of releases and an
assumption about survival after release:
$N_{y,a = 1}^{\text{juv,HO}} = N_{y}^{\text{rel}}s^{\text{rel}}$.

The age-1 NO fish is equal to the smolt production, usually considered
to be the juveniles migrating out of the freshwater environment:
$N_{y,a = 1}^{\text{juv,NO}} = \text{Smolt}_{y}$

> In the first year in the marine life stage, both HO and NO fish
> experience the natural mortality rate associated with $M$. Additional
> mortality can be applied to HO fish associated with release mortality
> (not experienced by NO fish). In effect, since the CWT data inform
> hatchery survival in the first year, the release mortality assumption
> effectively allows the analyst to attribute a smaller proportion of
> the mortality experienced by hatchery fish to natural origin fish.

The abundance of juvenile fish $N$, uses the same preterminal fishing
mortality, maturity, and natural mortality paramaters estimated from the
CWT: $$\begin{aligned}
N_{y + 1,a + 1}^{\text{juv,NO}} & {= N_{y,a}^{\text{juv,NO}}\exp\left( - v_{a}^{\text{PT}}F_{y}^{\text{PT}} \right)\left( 1 - m_{y,a} \right)\exp\left( - M_{y,a} \right)} \\
N_{y + 1,a + 1}^{\text{juv,HO}} & {= N_{y,a}^{\text{juv,HO}}\exp\left( - v_{a}^{\text{PT}}F_{y}^{\text{PT}} \right)\left( 1 - m_{y,a} \right)\exp\left( - M_{y,a} \right)}
\end{aligned}$$

### Adult life stage

The return $R$, and escapement $E$ is calculated from the terminal
fishing mortality and maturity parameters estimated from the CWT:

$$\begin{aligned}
R_{y,a}^{\text{NO}} & {= N_{y,a}^{\text{juv,NO}}\exp\left( - v_{a}^{\text{PT}}F_{y}^{\text{PT}} \right)m_{y,a}} \\
R_{y,a}^{\text{HO}} & {= N_{y,a}^{\text{juv,HO}}\exp\left( - v_{a}^{\text{PT}}F_{y}^{\text{PT}} \right)m_{y,a}}
\end{aligned}$$

$$\begin{aligned}
E_{y,a}^{\text{NO}} & {= R_{y,a}^{\text{NO}}\exp\left( - v_{a}^{\text{T}}F_{y}^{\text{T}} \right)} \\
E_{y,a}^{\text{HO}} & {= R_{y,a}^{\text{HO}}\exp\left( - v_{a}^{\text{T}}F_{y}^{\text{T}} \right)}
\end{aligned}$$

The number of spawners is calculated from the escapement, the en-route
survival, and the proportion of escapement that spawn. This proportion
can be calculated from the ratio of broodtake to observed escapement
(the model is parameterized in this manner for numerical stability and
it is assumed broodtake is non-selective with respect to origin):

$$\begin{aligned}
\text{NOS}_{y,a} & {= E_{y,a}^{\text{NO}}s^{\text{enroute}}p_{y}^{\text{spawn}}} \\
\text{HOS}_{y,a} & {= E_{y,a}^{\text{HO}}s^{\text{enroute}}p_{y}^{\text{spawn}}}
\end{aligned}$$

### Egg and juvenile production

The egg production is calculated from the proportion females and
fecundity $f$ at age:

$$\text{Egg}_{y} = p^{\text{female}}\sum\limits_{a}f_{a}\left( \text{NOS}_{y,a} + \gamma \times \text{HOS}_{y,a} \right)$$

The smolt production is calculated from a Ricker stock-recruit function:

$$\text{Smolt}_{y + 1} = \alpha \times \text{Egg}_{y} \times \exp\left( - \beta \times \text{Egg}_{y} \right)\exp\left( - \delta_{y} \right)$$

where $\delta_{y}$ is an annual deviation from the density-dependent
function, expressed as an instantaneous mortality rate. Here, smolts can
be considered to be outmigrating juveniles. Thus, the egg-smolt survival
occurs over freshwater life stages.

### Ricker stock-recruit function

[Walters and Korman
(2024)](https://publications.gc.ca/site/eng/9.940685/publication.html)
show that the realized egg-smolt mortality rate can be defined as

$$M_{y}^{\text{egg}} = M_{\text{min}}^{\text{egg}} + \beta \times \text{Egg}_{y} + \delta_{y}$$

where
$M_{\text{min}}^{\text{egg}} = - \log\left( \phi^{- 1} \right) - \log(\kappa)$
is the density-independent egg mortality rate when the egg production
approaches zero, $\beta \times \text{Egg}_{y}$ is the density-dependent
mortality term, and $\delta_{y}$ is a year-specific deviation.

The Ricker $\alpha$ parameter is $\alpha = \kappa/\phi$ in units of
juvenile-per-egg, where $\kappa$ is productivity (return/spawner) and
$\phi$ is the egg-per-juvenile at replacement (corresponding to the
one-to-one return per spawner line), with
$\phi = \sum_{a}\ell_{a}m_{a}^{\text{base}}f_{a}p^{\text{female}}$.

$\ell_{a}$ is the juvenile survival at age:

$$\ell_{a} = \begin{cases}
1 & {,a = 1} \\
{\ell_{a - 1}\exp\left( - M_{a - 1}^{\text{base}} \right)\left( 1 - m_{a - 1}^{\text{base}} \right)} & {,a = 2,\ldots,A}
\end{cases}$$

The $\beta$ parameter is
$\beta = \log(\kappa)/\left( \phi \times J^{\text{rep}} \right)$ is in
units of reciprocal eggs, where $J^{\text{rep}}$ is the juvenile smolt
production at replacement. The equilibrium spawners at replacement is
$S^{\text{rep}}$ and $J^{\text{rep}} = S^{\text{rep}}/\tau$, where
$\tau = \sum_{a}\ell_{a}m_{a}^{\text{base}}p^{\text{female}}$ is the
equilibrium female spawners per juvenile.

### Realized productivity

Note that the replacement values described above correspond to initial
maturity and natural mortality values. The realized replacement, along
with the realized productivity and spawners at replacement, varies if
either maturity and natural mortality are varying in the historical time
series because these parameters change the egg production over the
lifetime of a juvenile.

To obtain the realized productivity in a particular year, we can
calculate the egg-per-juvenile
$\phi_{y} = \sum_{a}\ell_{y,a}m_{y,a}f_{a}p^{\text{female}}$ and
spawner-per-juvenile
$\tau_{y} = \sum_{a}\ell_{y,a}m_{y,a}p^{\text{female}}$ from the
historical maturity and mortality, where

$$\ell_{y,a} = \begin{cases}
1 & {,a = 1} \\
{\ell_{y,a - 1}\exp\left( - M_{y,a - 1} \right)\left( 1 - m_{y,a - 1} \right)} & {,a = 2,\ldots,A}
\end{cases}$$

The realized productivity is $\kappa_{y} = \alpha\phi_{y}$ with
corresponding
$S_{y}^{\text{rep}} = \log\left( \kappa_{y} \right) \times \tau_{y}/\left( \beta\phi_{y} \right)$
and $E_{y}^{\text{rep}} = S_{y}^{\text{rep}}\phi_{y}/\tau_{y}$
($E_{y}^{\text{rep}}$ is the egg production at replacement).

The realized productivity can be less than one (and
$S_{y}^{\text{rep}} < 0$) in an individual year. If juvenile marine
mortality is high or maturity is late, then there are too few returns to
replace the population.

### Initial abundance

As a forward-projecting model, an assumption about the equilibrium age
distribution in the first year of the model is required. This age
distribution is calculated from the initial fishing mortality rate
$F^{\text{init}}$ of preterminal and terminal fisheries.

The initial equilibrium juvenile survival at age is:

$$\ell_{a}^{\text{Finit}} = \begin{cases}
1 & {,a = 1} \\
{\ell_{a - 1}\exp\left( - v_{a - 1}^{\text{PT}}F^{\text{init,PT}} \right)\exp\left( - M_{a - 1}^{\text{base}} \right)\left( 1 - m_{a - 1}^{\text{base}} \right)} & {,a = 2,\ldots,A}
\end{cases}$$

For the hatchery-origin population, the juvenile abundance is calculated
from an initial equilibrium release assumption $N_{eq}^{\text{rel}}$:

$$N_{y = 1,a}^{\text{juv,HO}} = N_{eq}^{\text{rel}} \times s^{\text{rel}} \times \ell_{a}^{\text{Finit}}$$

For the natural-origin population, two assumptions are available.

First, the juvenile abundance is in equilibrium corresponding to the
observed escapement $E_{y = 1}$ in the first year (feasible if there is
no hatchery population at the beginning of the model):

$$N_{y = 1,a}^{\text{juv,HO}} = \frac{E_{y = 1}}{\tau^{\text{Finit}}} \times \ell_{a}^{\text{Finit}}$$

where the juvenile abundance is calculated from the ratio of initial
escapement and spawner-per-juvenile
$\tau^{\text{Finit}} = \sum_{a}\ell_{a}^{\text{Finit}}\exp\left( - v_{a}^{\text{PT}}F^{\text{init,PT}} \right)\exp\left( - v_{a}^{\text{T}}F^{\text{init,T}} \right)m_{a}^{\text{base}}p^{\text{female}}$.

Alternatively, the natural-origin initial abundance is calculated from
the Ricker stock-recruit relationship where

$$N_{y = 1,a}^{\text{juv,NO}} = \frac{\log\left( \alpha \times \phi^{\text{Finit}} \right)}{\beta \times \phi^{\text{Finit}}} \times \ell_{a}^{\text{Finit}}$$

with the corresponding egg-per-juvenile
$\phi^{\text{Finit}} = \sum_{a}\ell_{a}^{\text{Finit}}\exp\left( - v_{a}^{\text{PT}}F^{\text{init,PT}} \right)\exp\left( - v_{a}^{\text{T}}F^{\text{init,T}} \right)m_{a}^{\text{base}}f_{a}p^{\text{female}}$.

## Parameter estimation and priors

### Fishing mortality

Year-specific fishing mortality is parameterized as

$$\begin{aligned}
F_{y}^{\text{PT}} & {= \exp\left( a^{\text{PT}} \right)F_{y}^{\text{trend,PT}}\exp\left( \omega_{y}^{\text{FPT}} \right)} \\
F_{y}^{\text{T}} & {= \exp\left( a^{\text{T}} \right)F_{y}^{\text{trend,T}}\exp\left( \omega_{y}^{\text{FT}} \right)} \\
 & 
\end{aligned}$$

where $F^{\text{trend,PT}}$ is a time series of relative exploitation
provided by the analyst. The model estimates a scaling coefficient $a$
and annual deviations $\omega$ to estimate fishing mortality.

The prior for the annual deviations is

$$\begin{aligned}
\omega_{y}^{\text{FPT}} & {\sim N\left( 0,\sigma_{\text{FPT}}^{2} \right)} \\
\omega_{y}^{\text{FT}} & {\sim N\left( 0,\sigma_{\text{FT}}^{2} \right)}
\end{aligned}$$

with hyperpriors for $\sigma_{\text{FPT}}^{2}$ and
$\sigma_{\text{FT}}^{2}$:

$$\begin{aligned}
\sigma_{\text{FPT}} & {\sim \text{Gamma}(2,5)} \\
\sigma_{\text{FT}} & {\sim \text{Gamma}(2,5)}
\end{aligned}$$

### Vulnerability

Vulnerability at age are independent terms estimated in logit space:

$$\begin{aligned}
{\text{logit}\left( v_{a}^{\text{PT}} \right)} & {\sim N\left( \text{logit}\left( \mu_{a}^{vPT} \right),1.6^{2} \right)} \\
{\text{logit}\left( v_{a}^{\text{T}} \right)} & {\sim N\left( \text{logit}\left( \mu_{a}^{vT} \right),1.6^{2} \right)}
\end{aligned}$$

where $\mu$ is the prior mean specified by the analyst. Vulnerability is
fixed at zero and one for the age 1 and maximum age ($A$), respectively
($v_{1} = 0$ and $v_{A} = 1$).

A useful prior may be with $\mu = 0.5$. When transformed to normal
space, the prior density is relatively uniform between 0-1 with low
density at the bounds:

``` r
x <- seq(-5, 5, 0.1)
f_x <- dnorm(x, 0, 1.6)

y <- plogis(x)
g_y <- f_x /(y * (1 - y)) # Prior density with Jacobian transformation

par(mfcol = c(1, 2), mar = c(5, 4, 1, 1))
plot(y, g_y, type = 'l', xlab = expression(v[a]), ylab = "Prior density")
plot(y, log(g_y), type = 'l', xlab = expression(v[a]), ylab = "log prior density")
```

![](conditioning_files/figure-html/unnamed-chunk-2-1.png)

### Maturity

Maturity is estimated in logit space as deviations from base parameters
provided by the analyst. The prior density function is Gaussian with
separate standard deviation $\sigma_{a}^{m}$ by age:

$$\text{logit}\left( m_{y,a} \right) \sim N\left( \text{logit}\left( m_{a}^{\text{base}} \right),\left\lbrack \sigma_{a}^{m} \right\rbrack^{2} \right)$$

with hyperprior $\sigma_{a}^{m} \sim \text{Gamma}(2,5)$ (mode of 0.2,
mean of 0.4, and standard deviation of 0.28).

``` r
x <- seq(0, 1.5, 0.01)
f_x <- dgamma(x, 2, scale = 1/5)

#x[which.max(f_x)] # Prior mode at 0.2

par(mfcol = c(1, 2), mar = c(5, 4, 1, 1), oma = c(0, 0, 2, 0))
plot(x, f_x, type = "l", ylab = "Prior density")
plot(x, log(f_x), type = "l", ylab = "log prior density")
mtext("Gamma(2, 5) prior distribution", outer = TRUE, line = 1)
```

![](conditioning_files/figure-html/unnamed-chunk-3-1.png)

### Natural mortality

Natural mortality is parameterized as

$$M_{y,a} = \begin{cases}
{M_{a}^{\text{base}} + \sum\limits_{i}X_{i,y}\beta_{i} + M^{\text{add}} + \varepsilon_{y}} & {,a = 1} \\
{M_{a}^{\text{base}} + \sum\limits_{j}X_{j,y}\beta_{j}} & {,a = 2,\ldots,A - 1}
\end{cases}$$

From base values provided by the analyst, year-specific mortality rates
can be estimated from linear combination of environmental covariates $X$
and estimated coefficients $\beta$. Separate covariates are used for
age-1 and age-2+ fish.

For age-1, an additional scalar $M^{\text{add}}$ and annual deviations
$\varepsilon_{y}$ can be estimated from the base parameter.

Gaussion priors are used for $\varepsilon_{y}$:

$$\varepsilon_{y} \sim N\left( 0,\left\lbrack \sigma^{M} \right\rbrack^{2} \right)$$

with hyperprior $\sigma^{M} \sim \text{Gamma}(2,5)$.

### Natural production

A uniform prior is used for $\log(\kappa)$, and a lognormal prior is
used to estimate $S_{0}$.

The annual deviations in smolt production are estimated with prior
$\delta_{y} \sim N\left( 0,\sigma_{\delta}^{2} \right)$ and hyperprior
$\sigma_{\delta} \sim \text{Gamma}(2,5)$.

## Likelihoods

### Total escapement

The likelihood of the total escapement for the population of interest
uses a lognormal distribution:

$$\log\left( E_{y} \right) \sim N\left( \log\left( \sum\limits_{a}\left( {\widehat{E}}_{y,a}^{\text{NO}} + {\widehat{E}}_{y,a}^{\text{HO}} \right) \right),\left\lbrack \sigma^{E} \right\rbrack^{2} \right)$$

### pHOS

The model can also be fitted to observations of pHOS (proportion
hatchery origin spawners by brood year $t$) for the population of
interest. A logistic normal distribution is used:

$$\log\left( \frac{p\text{HOS}_{t}}{1 - p\text{HOS}_{t}} \right) \sim N\left( \log\left( \frac{{\widehat{p\text{HOS}}}_{t}}{1 - {\widehat{p\text{HOS}}}_{t}} \right),\left\lbrack \sigma^{p\text{HOS}} \right\rbrack^{2} \right)$$

where the predicted value is

$$p\text{HOS}_{t} = \frac{\sum\limits_{a}\text{HOS}_{t + a - 1} \times f_{t + a - 1}}{\sum\limits_{a}(\left\lbrack \text{NOS}_{t + a - 1} + \text{HOS}_{t + a - 1}) \times f_{t + a - 1} \right\rbrack}$$

### CWT catch and escapement at age

The likelihood of the CWT at age uses a Poisson distribution:

$$\begin{aligned}
C_{y,a}^{\text{CWT,PT}} & {\sim \text{Poisson}\left( \frac{1}{\lambda} \times {\widehat{C}}_{y,a}^{\text{CWT,PT}} \right)} \\
C_{y,a}^{\text{CWT,T}} & {\sim \text{Poisson}\left( \frac{1}{\lambda} \times {\widehat{C}}_{y,a}^{\text{CWT,T}} \right)} \\
E_{y,a}^{\text{CWT}} & {\sim \text{Poisson}\left( \frac{1}{\lambda} \times {\widehat{E}}_{y,a}^{\text{CWT}} \right)}
\end{aligned}$$

where the $\land$ symbol denotes an estimate of the corresponding
observed quantity and $\lambda$ is the expansion factor.

The expansion factor is the scalar to ensure CWT data are representative
of total recoveries (typically greater or equal to 1). For example,
$\lambda = 10$ may be appropriate if 10 percent of the fishery catch is
sampled for coded wire tags. It is assumed that the CWT data are
unexpanded numbers.

### Likelihood re-weighting

CWT likelihood components can be reweighted relative to the total
escapement with alternative expansion factors (higher expansion factors
have lower weight with increasing variance for the Poisson
distribution).

The data needed to be adjusted simultaneously to preserve the magnitude
of predicted CWT catch and escapement between alternative fits. For
example, if the true sampling rate is 0.1, then an initial fit may use
$\lambda = 10$, but the CWT can be up-weighted assuming a higher
sampling rate of 0.5 ($\lambda = 2$). The analyst should multiply the
CWT catch and escapement by 5, or equivalently, multiply the releases by
0.2 with $\lambda = 1$.

## Release strategies

The equations above assume that CWT come from a single release strategy,
for example, fed fry, smolts, or yearlings (see [James et
al. 2021](https://www.marinescience.psf.ca/wp-content/uploads/2022/03/2021PSF-HatcheryReleaseStrategies-BCReview-FinalReport-Web.pdf)
for a review).

The conditioning model can accommodate multiple release strategies. The
number of CWT releases are specified by brood year and release strategy,
and the model is fitted to CWT catch by year, age, and release strategy.

Fishing and natural mortality rates are identical by year and age but
maturity differs by release strategy.

The maturity is estimated by year and age for a particular release
strategy indexed by $r = 1$:

$$\text{logit}\left( m_{y,a,r} \right) \sim N\left( \text{logit}\left( m_{a}^{\text{base}} \right),\left\lbrack \sigma_{a}^{m} \right\rbrack^{2} \right)$$

For all other release strategies $r = 2,\ldots,n_{r}$, logit deviations
by age (constant with time) are estimated as fixed effects:

$$\text{logit}\left( m_{y,a,r} \right) = \text{logit}\left( m_{y,a,1} \right) + \delta_{a,r}$$

The user then specifies the release strategy from which the maturity
estimates then apply towards the population of interest.

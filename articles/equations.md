# Dynamics equations

salmonMSE utilizes an age-structured model in the projections. The
population is tracked by age and year but various dynamics correspond to
the salmon life stages as described below.

## Variable definitions

*Definition of variable names and the corresponding slots in either the
input (SOM) or output (SMSE) objects in salmonMSE.*

### Natural production

| Name                        | Definition                                                                                                      | Type                          | Class    | Slot               |
|:----------------------------|:----------------------------------------------------------------------------------------------------------------|:------------------------------|:---------|:-------------------|
| $\text{NOS}$                | Natural origin spawners                                                                                         | Natural production            | SMSE     | NOS                |
| $\text{Fry}^{\text{NOS}}$   | Fry production by natural origin spawners, assumed to be equal to egg production                                | Natural production            | SMSE     | Fry_NOS            |
| $\text{Smolt}^{\text{NOS}}$ | Smolt production by natural origin spawners, density-dependent                                                  | Natural production            | SMSE     | Smolt_NOS          |
| $C_{\text{egg-smolt}}$      | Carrying capacity of smolts (Beverton-Holt stock-recruit parameter)                                             | Natural production            | SOM, Bio | capacity           |
| $S_{\text{max}}$            | Spawning output that maximizes smolt production (Ricker stock-recruit parameter)                                | Natural production            | SOM, Bio | Smax               |
| $\kappa$                    | Productivity (maximum recruitment production rate), units of recruit per spawner                                | Natural production            | SOM, Bio | kappa              |
| $\phi$                      | Unfished per capita egg production rate, units of egg per smolt                                                 | Natural production            | SOM, Bio | phi                |
| $r$                         | Maturity at age, i.e., recruitment rate                                                                         | Natural production            | SOM, Bio | p_mature           |
| $\text{Fec}$                | Fecundity of spawners (eggs per female)                                                                         | Natural production            | SOM, Bio | fec                |
| $p^{\text{female}}$         | Proportion of female spawners in broodtake and spawners                                                         | Natural production            | SOM, Bio | p_female           |
| $\text{SAR}$                | Smolt-to-adult recruit survival                                                                                 | Natural production            | \-       | \-                 |
| $M$                         | Juvenile instantaneous natural mortality of juvenile (either the freshwater or marine environment by age class) | Natural production + Hatchery | SOM, Bio | Mjuv_NOS, Mjuv_HOS |
| $s_{\text{enroute}}$        | Survival of escapement to spawning grounds and hatchery                                                         | Natural production            | SOM, Bio | s_enroute          |
| $\text{NOR}$                | Natural origin return                                                                                           | Natural production            | SMSE     | Return_NOS         |
| $p_{\text{HOSeff}}$         | Proportion of effective hatchery origin spawners (vs. NOS)                                                      | Population dynamics           | SMSE     | pHOS_effective     |
| $p_{\text{HOScensus}}$      | Proportion of hatchery origin spawners (vs. NOS)                                                                | Population dynamics           | SMSE     | pHOS_census        |
| $p^{\text{WILD}}$           | Proportion of wild spawners                                                                                     | Population dynamics           | SMSE     | p_wild             |

### Habitat

| Name                                 | Definition                                                                       | Type    | Class | Slot           |
|:-------------------------------------|:---------------------------------------------------------------------------------|:--------|:------|:---------------|
| $P^{\text{inc}}$                     | Productivity for density-dependent survival: egg incubation from spawning output | Habitat | SOM   | egg_prod       |
| $C^{\text{inc}}$                     | Capacity for density-dependent survival: egg incubation from spawning output     | Habitat | SOM   | egg_capacity   |
| $P^{\text{egg-fry}}$                 | Productivity for density-dependent survival: egg to fry life stage               | Habitat | SOM   | fry_prod       |
| $C^{\text{egg-fry}}$                 | Capacity for density-dependent survival: egg to fry life stage                   | Habitat | SOM   | fry_capacity   |
| $\varepsilon_{y}^{\text{egg-fry}}$   | Deviations in density-dependent survival: egg to fry life stage                  | Habitat | SOM   | fry_sdev       |
| $P^{\text{fry-smolt}}$               | Productivity for density-dependent survival: fry to smolt life stage             | Habitat | SOM   | smolt_prod     |
| $C^{\text{fry-smolt}}$               | Capacity for density-dependent survival: fry to smolt life stage                 | Habitat | SOM   | smolt_capacity |
| $\varepsilon_{y}^{\text{fry-smolt}}$ | Deviations in density-dependent survival: fry to smolt life stage                | Habitat | SOM   | smolt_sdev     |

### Hatchery

| Name                                 | Definition                                                                                                             | Type                          | Class                           | Slot               |
|:-------------------------------------|:-----------------------------------------------------------------------------------------------------------------------|:------------------------------|:--------------------------------|:-------------------|
| $\text{HOS}$                         | Hatchery origin spawners                                                                                               | Hatchery                      | SMSE                            | HOS                |
| $\text{HOS}_{\text{eff}}$            | Effective number of HOS, spawning output discounted by $\gamma$                                                        | Hatchery                      | SMSE                            | HOSeff             |
| $\text{Fry}^{\text{HOS}}$            | Fry production by hatchery origin spawners, assumed to be equal to egg production                                      | Hatchery                      | SMSE                            | Fry_HOS            |
| $\text{Smolt}^{\text{HOS}}$          | Smolt production by hatchery origin spawners, density-dependent                                                        | Hatchery                      | SMSE                            | Smolt_HOS          |
| $\text{Fec}^{\text{brood}}$          | Fecundity of broodtake (eggs per female)                                                                               | Hatchery                      | SOM                             | fec_brood          |
| $M$                                  | Juvenile instantaneous natural mortality of juvenile (either the freshwater or marine environment by age class)        | Natural production + Hatchery | SOM, Bio                        | Mjuv_NOS, Mjuv_HOS |
| $\text{NOB}$                         | Natural origin broodtake                                                                                               | Hatchery                      | SMSE                            | NOB                |
| $\text{HOB}$                         | Hatchery origin broodtake                                                                                              | Hatchery                      | SMSE                            | HOB                |
| $\text{Stray}$                       | External strays of hatchery origin fish (considered 0% marked)                                                         | Hatchery                      | SOM                             | stray_external     |
| $\text{HOB\_stray}$                  | Broodtake from strays                                                                                                  | Hatchery                      | SMSE                            | HOB_stray          |
| $\text{Brood}^{\text{avail,import}}$ | Available imported brood (considered 100% marked HO fish)                                                              | Hatchery                      | SOM                             | brood_import       |
| $\text{Brood}^{\text{import}}$       | Realized imported fish used for brood (considered 100% marked HO)                                                      | Hatchery                      | SMSE                            | HOB_import         |
| $s_{\text{yearling}}$                | Survival of hatchery eggs to yearling life stage                                                                       | Hatchery                      | SOM                             | s_egg_smolt        |
| $s_{\text{subyearling}}$             | Survival of hatchery eggs to subyearling life stage                                                                    | Hatchery                      | SOM                             | s_egg_subyearling  |
| $p_{\text{yearling}}$                | Proportion of hatchery releases as yearling (vs. subyearling)                                                          | Hatchery                      | Internal state variable         | \-                 |
| $s_{\text{prespawn}}$                | Survival of adult broodtake in hatchery                                                                                | Hatchery                      | SOM                             | s_prespawn         |
| $n_{\text{yearling}}$                | Target number of hatchery releases as yearlings                                                                        | Hatchery                      | SOM                             | n_yearling         |
| $n_{\text{subyearling}}$             | Target number of hatchery releases as subyearlings                                                                     | Hatchery                      | SOM                             | n_subyearling      |
| $m$                                  | Mark rate of hatchery fish                                                                                             | Hatchery                      | SOM                             | m                  |
| $p_{\text{max}}^{\text{esc}}$        | Maximum proportion of total escapement (after en-route mortality) to use as broodtake                                  | Hatchery                      | SOM                             | pmax_esc           |
| $p_{\text{target}}^{\text{NOB}}$     | Target proportion of the natural origin broodtake from the escapement (after en-route mortality), i.e., NOB/NOS ratio  | Hatchery                      | SOM                             | ptarget_NOB        |
| $p_{\text{max}}^{\text{NOB}}$        | Maximum proportion of the natural origin broodtake from the escapement (after en-route mortality), i.e., NOB/NOS ratio | Hatchery                      | SOM                             | pmax_NOB           |
| $p_{\text{NOB}}$                     | Realized proportion of the total broodtake of hatchery origin (vs. natural origin)                                     | Hatchery                      | SMSE                            | pNOB               |
| $\text{HOR}$                         | Hachery origin return                                                                                                  | Hatchery                      | SMSE                            | Return_HOS         |
| $p^{\text{hatchery}}$                | Proportion of hatchery origin escapement to hatchery, available for broodtake                                          | Hatchery                      | SOM                             | phatchery          |
| $p_{\text{removal}}^{\text{HOS}}$    | Proportion of hatchery origin fish removed from spawning grounds, not available for broodtake                          | Hatchery                      | SOM                             | premove_HOS        |
| $p_{\text{removal}}^{\text{NOS}}$    | Proportion of natural origin fish removed from spawning grounds, not available for broodtake                           | Hatchery                      | SOM                             | premove_NOS        |
| $\gamma$                             | Reduced reproductive success of HOS (relative to NOS)                                                                  | Hatchery                      | SOM                             | gamma              |
| $\bar{z}$                            | Mean phenotypic value of cohort in natural and hatchery environments                                                   | Fitness                       | Internal state variable and SOM | zbar_start         |
| $\theta$                             | Optimal phenotypic value for natural and hatchery environments                                                         | Fitness                       | SOM                             | theta              |
| $\sigma^{2}$                         | Variance of phenotypic traits in population                                                                            | Fitness                       | SOM                             | phenotype_variance |
| $\omega^{2}$                         | Variance of fitness function                                                                                           | Fitness                       | SOM                             | fitness_variance   |
| $h^{2}$                              | Heritability of phenotypic traits                                                                                      | Fitness                       | SOM                             | heritability       |
| $\bar{W}$                            | Population fitness in the natural and hatchery environments                                                            | Fitness                       | SMSE                            | fitness            |
| $\ell_{i}$                           | Relative fitness loss at the life stage i (egg, fry, smolt)                                                            | Fitness                       | SOM                             | rel_loss           |
| $\text{PNI}$                         | Proportionate natural influence                                                                                        | Fitness                       | SMSE                            | PNI                |

### Harvest

| Name            | Definition                                                                                       | Type    | Class | Slot          |
|:----------------|:-------------------------------------------------------------------------------------------------|:--------|:------|:--------------|
| $m$             | Mark rate of hatchery fish (affects fishery retention of hatchery fish relative to natural fish) | Harvest | SOM   | m             |
| $u^{\text{PT}}$ | Pre-terminal fishery harvest rate                                                                | Harvest | SOM   | u_preterminal |
| $u^{\text{T}}$  | Terminal fishery harvest rate                                                                    | Harvest | SOM   | u_terminal    |
| $\delta$        | Mortality from catch and release (proportion)                                                    | Harvest | SOM   | release_mort  |
| $v$             | Relative vulnerability by age to the fishery                                                     | Harvest | SOM   | vulPT, vulT   |

## Natural production

First, we consider natural production in the absence of fitness effects
arising from hatchery production.

### Spawning output

From the spawners (NOS and HOS) of age $a$ in year $y$, the
corresponding spawning output (units of eggs) of the subsequent
generation is calculated as:

$$\begin{aligned}
\text{Egg}_{y}^{\text{NOS}} & {= \sum\limits_{a}\text{NOS}_{y,a} \times p^{\text{female}} \times \text{Fec}_{a}} \\
\text{Egg}_{y}^{\text{HOS}} & {= \sum\limits_{a}\text{HOS}_{\text{eff}y,a} \times p^{\text{female}} \times \text{Fec}_{a}}
\end{aligned}$$

where $\text{HOS}_{\text{eff}} = \gamma \times \text{HOS}$ and the
superscript denotes the parentage of the progeny.

### Smolt production - no habitat modeling

If no habitat modeling is used, then fry production is assumed to be
equal to spawning output, i.e.,
$\text{Fry}_{y + 1}^{\text{NOS}} = \text{Egg}_{y}^{\text{NOS}}$ and
$\text{Fry}_{y + 1}^{\text{HOS}} = \text{Egg}_{y}^{\text{HOS}}$.

Survival from egg to smolt life stage is density-dependent. With the
Beverton-Holt stock-recruit relationship, the age-1 smolt production is

$$\begin{aligned}
\text{Smolt}_{y + 1}^{\text{NOS}} & {= \frac{\alpha \times \text{Fry}_{y + 1}^{\text{NOS}}}{1 + \beta\left( \text{Fry}_{y + 1}^{\text{NOS}} + \text{Fry}_{y + 1}^{\text{HOS}} + n_{y + 1}^{\text{sub}} \right)}} \\
\text{Smolt}_{y + 1}^{\text{HOS}} & {= \frac{\alpha \times \text{Fry}_{y + 1}^{\text{HOS}}}{1 + \beta\left( \text{Fry}_{y + 1}^{\text{NOS}} + \text{Fry}_{y + 1}^{\text{HOS}} + n_{y + 1}^{\text{sub}} \right)}}
\end{aligned}$$

where $\alpha = \kappa/\phi$, $\beta = \alpha/C_{\text{egg-smolt}}$, the
unfished egg per smolt
$\phi = \sum_{a}\left( \prod_{i = 1}^{a - 1}\exp\left( - M_{i}^{\text{NOS}} \right)\left( 1 - r_{i} \right) \right) \times r_{a} \times p^{\text{female}} \times \text{Fec}_{a}$,
with $r_{a}$ as the maturity at age.

> Smolt production can be predicted from total adult spawners by setting
> $\text{Fec}_{a} = 1$ and $\phi = 1$.

The density-independent component of the survival equation is controlled
by $\alpha$ and the density-dependent component of survival is
controlled by $\beta$ and scaled by the total number of fry in
competition with subyearling hatchery releases (see
[Hatchery](#hatchery-production) section).

If there is knife-edge maturity, i.e., all fish mature at the terminal
age, the equation simplifies to
$\phi = \text{SAR} \times p^{\text{female}} \times \text{Fec}$, with
$\text{SAR}$ as the marine survival (between 0-1).

> With the Ricker stock-recruit relationship, smolt production is
>
> $$\begin{aligned}
> \text{Smolt}_{y + 1}^{\text{NOS}} & {= \alpha \times \text{Fry}_{y + 1}^{\text{NOS}} \times \exp\left( - \beta\left\lbrack \text{Fry}_{y + 1}^{\text{NOS}} + \text{Fry}_{y + 1}^{\text{HOS}} + n_{y + 1}^{\text{sub}} \right\rbrack \right)} \\
> \text{Smolt}_{y + 1}^{\text{HOS}} & {= \alpha \times \text{Fry}_{y + 1}^{\text{HOS}} \times \exp\left( - \beta\left\lbrack \text{Fry}_{y + 1}^{\text{NOS}} + \text{Fry}_{y + 1}^{\text{HOS}} + n_{y + 1}^{\text{sub}} \right\rbrack \right)}
> \end{aligned}$$
>
> where $\alpha = \kappa/\phi$ and $\beta = 1/S_{\text{max}}$,
> $S_{\text{max}}$ is the egg production that maximizes smolt
> production.

### Smolt production - habitat modeling

Egg to smolt production can also be modeled as a series of
density-dependent functions by life stage, following the approach of
[Jorgensen et al. 2021](https://doi.org/10.1371/journal.pone.0256792).
Three relationships are modeled.

The realized egg production (${\widetilde{\text{Egg}}}_{y}$) can be
modified from the spawning output ($\text{Egg}_{y}$) due to incubation
mortality. With a Beverton-Holt function:

$$\begin{aligned}
{\widetilde{\text{Egg}}}_{y}^{\text{NOS}} & {= \frac{P^{\text{inc}} \times \text{Egg}_{y}^{\text{NOS}}}{1 + \frac{P^{\text{inc}}}{C^{\text{inc}}}\left( \text{Egg}_{y}^{\text{NOS}} + \text{Egg}_{y}^{\text{HOS}} \right)}} \\
{\widetilde{\text{Egg}}}_{y}^{\text{HOS}} & {= \frac{P^{\text{inc}} \times \text{Egg}_{y}^{\text{HOS}}}{1 + \frac{P^{\text{inc}}}{C^{\text{inc}}}\left( \text{Egg}_{y}^{\text{NOS}} + \text{Egg}_{y}^{\text{HOS}} \right)}}
\end{aligned}$$

where productivity $P$ is the maximum survival as spawning output
approaches zero and $C$ is the asymptotic production.

> Set the capacity to infinite to model density-independence. The
> productivity parameter is then the survival to the next life stage.

Fry production is modeled as:

$$\begin{aligned}
\text{Fry}_{y + 1}^{\text{NOS}} & {= \frac{P^{\text{egg-fry}} \times {\widetilde{\text{Egg}}}_{y}^{\text{NOS}}}{1 + \frac{P^{\text{egg-fry}}}{C^{\text{egg-fry}}}\left( {\widetilde{\text{Egg}}}_{y}^{\text{NOS}} + {\widetilde{\text{Egg}}}_{y}^{\text{HOS}} \right)} \times \varepsilon_{y}^{\text{egg-fry}}} \\
\text{Fry}_{y + 1}^{\text{HOS}} & {= \frac{P^{\text{egg-fry}} \times {\widetilde{\text{Egg}}}_{y}^{\text{HOS}}}{1 + \frac{P^{\text{egg-fry}}}{C^{\text{egg-fry}}}\left( {\widetilde{\text{Egg}}}_{y}^{\text{NOS}} + {\widetilde{\text{Egg}}}_{y}^{\text{HOS}} \right)} \times \varepsilon_{y}^{\text{egg-fry}}}
\end{aligned}$$

where $\varepsilon_{y}^{\text{egg-fry}}$ is a year-specific deviation in
survival. They can be modeled as a function of a proposed time series of
environmental variables $\eta$, for example,
$\varepsilon_{y}^{\text{egg-fry}} = \prod_{j}f\left( \eta_{y,j} \right)$
or
$\varepsilon_{y}^{\text{egg-fry}} = \sum_{j}f\left( \eta_{y,j} \right)$.

Similarly, smolt production is modeled as:

$$\begin{aligned}
\text{Smolt}_{y}^{\text{NOS}} & {= \frac{P^{\text{fry-smolt}} \times \text{Fry}_{y}^{\text{NOS}}}{1 + \frac{P^{\text{fry-smolt}}}{C^{\text{fry-smolt}}}\left( \text{Fry}_{y}^{\text{NOS}} + \text{Fry}_{y}^{\text{HOS}} \right)} \times \varepsilon_{y}^{\text{fry-smolt}}} \\
\text{Smolt}_{y}^{\text{HOS}} & {= \frac{P^{\text{fry-smolt}} \times \text{Fry}_{y}^{\text{HOS}}}{1 + \frac{P^{\text{fry-smolt}}}{C^{\text{fry-smolt}}}\left( \text{Fry}_{y}^{\text{NOS}} + \text{Fry}_{y}^{\text{HOS}} \right)} \times \varepsilon_{y}^{\text{fry-smolt}}}
\end{aligned}$$

Alternative scenarios with changes in productivity or capacity
parameters can be used to evaluate changes in life stage survival from
habitat improvement or mitigation measures as part of a management
strategy, or from climate regimes (low productivity vs. high
productivity, or low capacity vs. high capacity). An increase in
capacity can arise from restoration which increases the area of suitable
habitat. An increase in productivity can arise from improvement in
habitat, e.g., sediment quality.

Approaches such as
[HARP](https://www.fisheries.noaa.gov/resource/tool-app/habitat-assessment-and-restoration-planning-harp-model)
and
[CEMPRA](https://www.essa.com/explore-essa/projects/cumulative-effects-model-for-priority-of-recovery-actions-cempra/)
can inform productivity and capacity parameters across these life stages
as quantitative relationships between habitat variables.

> For all life stages, a hockey-stick formulation is also possible. For
> example:
>
> $${\widetilde{\text{Egg}}}_{y}^{\text{NOS}} = \begin{cases}
> {P^{\text{inc}} \times \text{Egg}_{y}^{\text{NOS}}} & {,\text{Egg}_{y}^{\text{NOS}} \leq C_{y}^{\text{inc*}}/P^{\text{inc}}} \\
> C_{y}^{\text{inc*}} & {,\text{otherwise}} \\
>  & 
> \end{cases}$$
>
> where
> $C_{y}^{\text{inc*}} = C^{\text{inc}} \times \text{Egg}_{y}^{\text{NOS}}/\left( \text{Egg}_{y}^{\text{NOS}} + \text{Egg}_{y}^{\text{HOS}} \right)$
> is the capacity apportioned to natural spawners based on relative
> abundance.

## Hatchery production

Hatchery production is controlled by several sets of variables specified
by the analyst, roughly following the
[AHA](https://www.streamnet.org/home/data-maps/hatchery-reform/hsrg-tools/)
approach.

The first consideration is to specify the target number of annual
releases of sub-yearlings $n_{\text{target}}^{\text{subyearling}}$ and
yearlings $n_{\text{target}}^{\text{yearling}}$.

Yearlings are intended to represent hatchery releases that immediately
leave freshwater environment, while subyearlings are subject to
density-dependent survival in competition with natural production of
fry, e.g., they reside in freshwater environment for a period of time
before leaving.

Going backwards, the corresponding number of eggs needed to reach the
target number depends on the egg survival to those life stages in the
hatchery. The corresponding number of broodtake is calculated from
target egg production based on the brood fecundity and hatchery survival
of broodtake, which is non-selective with respect to age.

An additional consideration is the composition (natural vs. hatchery
origin) of in-river broodtake. To minimize genetic drift of the
population due to hatchery production, it is desirable to maintain a
high proportion of natural origin broodtake. This is controlled by
$p_{\text{target}}^{\text{NOB}}$, the desired proportion of natural
broodtake relative to all broodtake (any specified amount of available
imported brood is considered hatchery-origin for this purpose), but can
be exceeded if there is insufficient escapement of natural origin fish.

The ability to meet this target depends on the mark rate of hatchery
origin fish. Thus, $p_{\text{target}}^{\text{NOB}}$ represents ratio of
unmarked fish in the projection (imported brood is considered marked for
this calculation, strays are considered unmarked), and the realized
$p^{\text{NOB}}$ is reduced by the mark rate.

Another consideration for broodtake dynamics is to maintain high
spawning of natural origin fish. This is controlled by
$p_{\text{max}}^{\text{NOB}}$, the maximum allowable proportion of the
natural origin escapement to be used as broodtake. This value is never
exceeded.

> To set up a segregated hatchery program, set
> $p_{\text{max}}^{\text{NOB}} = 0$. Otherwise, these equations set up
> an integerated hatchery.

The following equations then generate the annual broodtake and hatchery
production from the state variables given these constraints.

### Broodtake

The annual target egg production for the hatchery is calculated from the
target releases as

$$\text{Egg}_{\text{target,broodtake}} = \frac{n_{\text{target}}^{\text{yearling}}}{s^{\text{yearling}}} + \frac{n_{\text{target}}^{\text{subyearling}}}{s^{\text{subyearling}}}$$

where $s$ is the corresponding survival term from the egg life stage.

The broodtake is back-calculated from the target egg production. The
composition of natural and hatchery origin broodtake (NOB and HOB,
respectively) is dependent on the mark rate $m$ and the target
proportion of NOB $p_{\text{target}}^{\text{NOB}}$. When the mark rate
is 1, then the realized pNOB should be equal to
$p_{\text{target}}^{\text{NOB}}$ provided there is sufficient
escapement. If the mark rate is less than one, then
$p_{\text{target}}^{\text{NOB}}$ reflects the proportion of unmarked
fish in the broodtake, some which are hatchery origin. Thus, the
realized pNOB is less than $p_{\text{target}}^{\text{NOB}}$. If the mark
rate is zero, then broodtake is non-selective with pNOB equal to the
proportion of natural origin escapement.

From the escapement in year $y$, some proportion $p^{\text{broodtake}}$
is used as broodtake:

$$\begin{aligned}
\text{NOB}_{y,a} & {= p_{y}^{\text{broodtake,unmarked}} \times \text{NOR}_{y,a}^{\text{escapement}} \times s_{\text{enroute}} \times p_{\text{max}}^{\text{esc}}} \\
\text{HOB}_{y,a}^{\text{unmarked}} & {= p_{y}^{\text{broodtake,unmarked}} \times (1 - m) \times p^{\text{hatchery}} \times \text{HOR}_{y,a}^{\text{escapement}} \times s_{\text{enroute}} \times p_{\text{max}}^{\text{esc}}} \\
\text{HOB}_{y,a}^{\text{marked}} & {= p_{y}^{\text{broodtake,marked}} \times m \times p^{\text{hatchery}} \times \text{HOR}_{y,a}^{\text{escapement}} \times s_{\text{enroute}} \times p_{\text{max}}^{\text{esc}}}
\end{aligned}$$

The proportion of the available hatchery fish for broodtake is also
reduced by $p^{\text{hatchery}}$, which can include fish swimming back
to the hatchery or removed from spawning grounds.

Additionally, some proportion of imported fish and strays may be used as
brood:

$$\begin{aligned}
\text{Brood}_{y,a}^{\text{import}} & {= p_{y}^{\text{broodtake,marked}}\sum\limits_{a}\text{Brood}_{a}^{\text{avail,import}}} \\
\text{HOB}_{y,a}^{\text{stray}} & {= p_{y}^{\text{broodtake,unmarked}} \times \text{Stray}_{y,a} \times s_{\text{enroute}}} \\
 & 
\end{aligned}$$

The availability of both natural and hatchery origin fish depends on the
escapement reduced by en-route mortality and can be capped by some
proportion denoted by the $p_{\text{max}}^{\text{esc}}$ parameter.

> To exclusively use imported brood, set
> $p_{\text{max}}^{\text{esc}} = 0$.

The realized hatchery egg production is

$$\begin{aligned}
\text{Egg}_{\text{y}}^{\text{NOB}} & {= \sum\limits_{a}\text{NOB}_{y,a} \times s^{\text{prespawn}} \times p^{\text{female}} \times \text{Fec}_{a}^{\text{brood}}} \\
\text{Egg}_{\text{y}}^{\text{HOB}} & {= \sum\limits_{a}\left( \text{HOB}_{y,a}^{\text{marked}} + \text{HOB}_{y,a}^{\text{unmarked}} \right) \times s^{\text{prespawn}} \times p^{\text{female}} \times \text{Fec}_{a}^{\text{brood}}} \\
\text{Egg}_{\text{y}}^{\text{import}} & {= \sum\limits_{a}\text{Brood}_{y,a}^{\text{import}} \times s^{\text{prespawn}} \times p^{\text{female}} \times \text{Fec}_{a}^{\text{brood}}} \\
\text{Egg}_{\text{y}}^{\text{stray}} & {= \sum\limits_{a}\text{HOB}_{y,a}^{\text{stray}} \times s^{\text{prespawn}} \times p^{\text{female}} \times \text{Fec}_{a}^{\text{brood}}} \\
 & 
\end{aligned}$$

where hatchery egg production is subject to a survival term
$s^{\text{prespawn}}$.

The proportion $p_{y}^{\text{broodtake}}$ is solved annually to satisfy
the following conditions:

$\frac{\sum_{a}\left( \text{NOB}_{y,a} + \text{HOB}_{y,a}^{\text{unmarked}} + \text{HOB}_{y,a}^{\text{stray}} \right)}{\sum_{a}\left( \text{NOB}_{y,a} + \text{HOB}_{y,a}^{\text{unmarked}} + \text{HOB}_{y,a}^{\text{marked}} + \text{HOB}_{y,a}^{\text{stray}} + \text{Brood}_{y,a}^{\text{import}} \right)} = p_{\text{target}}^{\text{NOB}}$

$0 < p_{y}^{\text{broodtake,marked}} \leq 1$

$0 < p_{y}^{\text{broodtake,unmarked}} \leq p_{\text{max}}^{\text{NOB}}$

$\text{Egg}_{\text{y}}^{\text{NOB}} + \text{Egg}_{\text{y}}^{\text{HOB}} + \text{Egg}_{\text{y}}^{\text{import}} + \text{Egg}_{y}^{\text{stray}} = \text{Egg}_{\text{broodtake}}$

The target ratio $p_{\text{target}}^{\text{NOB}}$ reflects the objective
to maintain a high proportion of natural origin fish in the broodtake,
where its implementation is dependent on the mark rate. The maximum
removal rate of natural origin fish $p_{\text{max}}^{\text{NOB}}$ or
escapement $p_{\text{max}}^{\text{esc}}$ ensures that there is high
abundance of natural origin spawners.

The total egg production in a given year can fail to reach the target if
there is insufficient unmarked escapement. In this case, the unmarked
take is set to the maximum removal rate
($p_{y}^{\text{broodtake,unmarked}} = p_{\text{max}}^{\text{NOB}}$), and
the remaining deficit in egg production is met using HOB (including
strays and imports).

### Smolt releases

After the total hatchery egg production is calculated, the production of
yearlings and subyearlings is calculated to ensure the annual ratio is
equal to the target ratio. To do so, the parameter
$p_{y}^{\text{egg,yearling}}$ is solved subject to the following
conditions:

$\text{Egg}_{\text{brood,y}} = \text{Egg}_{\text{y}}^{\text{NOB}} + \text{Egg}_{\text{y}}^{\text{HOB}} + \text{Egg}_{y}^{\text{import}} + \text{Egg}_{y}^{\text{stray}}$

$n_{y + 1}^{\text{yearling}} = p_{y}^{\text{egg,yearling}} \times \text{Egg}_{\text{brood,y}} \times s^{\text{yearling}}$

$n_{y + 1}^{\text{subyearling}} = \left( 1 - p_{y}^{\text{egg,yearling}} \right) \times \text{Egg}_{\text{brood,y}} \times s^{\text{subyearling}}$

$\frac{n_{y}^{\text{yearling}}}{n_{y}^{\text{subyearling}} + n_{y}^{\text{yearling}}} = \frac{n_{\text{target}}^{\text{yearling}}}{n_{\text{target}}^{\text{subyearling}} + n_{\text{target}}^{\text{yearling}}}$

From the total broodtake, the smolt releases is calculated as

$$\text{Smolt}_{y + 1}^{\text{Rel}} = n_{y + 1}^{\text{yearling}} + \frac{\alpha \times n_{y + 1}^{\text{subyearling}}}{1 + \beta\left( \text{Fry}_{y + 1}^{\text{NOS}} + \text{Fry}_{y + 1}^{\text{HOS}} + n_{y + 1}^{\text{subyearling}} \right)}$$

or

$$\text{Smolt}_{y + 1}^{\text{Rel}} = n_{y + 1}^{\text{yearling}} + \alpha \times n_{y + 1}^{\text{subyearling}} \times \exp\left( - \beta\left( \text{Fry}_{y + 1}^{\text{NOS}} + \text{Fry}_{y + 1}^{\text{HOS}} + n_{y + 1}^{\text{subyearling}} \right) \right)$$

## Pre-terminal fishery

Let $N_{y,a}^{\text{juv}}$ be the juvenile abundance in the population
and
$N_{y,a = 1}^{\text{juv,NOS}} = \text{Smolt}_{y}^{\text{NOS}} + \text{Smolt}_{y}^{\text{HOS}}$
and $N_{y,a = 1}^{\text{juv,HOS}} = \text{Smolt}^{\text{Rel}}$. The
superscript for the smolt variable corresponds to the parentage while
the superscript for $N$ denotes the origin of the current cohort.

Harvest $u^{\text{PT}}$ in the pre-terminal ($\text{PT}$) fishery,
assuming no mark-selective fishing, is modeled as a seasonal process and
occurs in the first half of the year.

The kept catch $K$ is

$$\begin{aligned}
K_{y,a}^{\text{NOS,PT}} & {= \left( 1 - \exp\left( - v_{a}^{\text{PT}}F_{y}^{\text{PT}} \right) \right)N_{y,a}^{\text{juv,NOS}}} \\
K_{y,a}^{\text{HOS,PT}} & {= \left( 1 - \exp\left( - v_{a}^{\text{PT}}F_{y}^{\text{PT}} \right) \right)N_{y,a}^{\text{juv,HOS}}} \\
 & 
\end{aligned}$$

If management is by target harvest rate, then the specified harvest rate
$u^{\text{PT}}$ in the pre-terminal ($\text{PT}$) fishery is converted
to the apical instantaneous fishing mortality rates as
$F^{\text{PT}} = - \log\left( 1 - u^{\text{PT}} \right)$.

If management is by target catch $K^{\text{PT,target}}$, then the
fishing mortality is solved such that
$K^{\text{PT,target}} = \sum_{a}K_{y,a}^{\text{NOS,PT}} + \sum_{a}K_{y,a}^{\text{HOS,PT}}$.

Alternative equations are used if [mark-selective
fishing](#mark-selective-fishing) is implemented.

## Recruitment and maturity

The recruitment is calculated from the survival of juvenile fish after
pre-terminal harvest and maturation:

$$\begin{aligned}
\text{NOR}_{y,a} & {= N_{y,a}^{\text{juv,NOS}}\exp\left( - v_{a}F_{y}^{\text{PT}} \right)r_{y,a}} \\
\text{HOR}_{y,a} & {= N_{y,a}^{\text{juv,HOS}}\exp\left( - v_{a}F_{y}^{\text{PT}} \right)r_{y,a}}
\end{aligned}$$

The juvenile abundance in the following year consists of fish that did
not mature and subsequently survived natural mortality $M$:

$$\begin{aligned}
N_{y + 1,a + 1}^{\text{juv,NOS}} & {= N_{y,a}^{\text{juv,NOS}}\exp\left( - \left\lbrack v_{a}F_{y}^{\text{PT}} + M_{y,a}^{\text{NOS}} \right\rbrack \right)\left( 1 - r_{y,a} \right)} \\
N_{y + 1,a + 1}^{\text{juv,HOS}} & {= N_{y,a}^{\text{juv,HOS}}\exp\left( - \left\lbrack v_{a}F_{y}^{\text{PT}} + M_{y,a}^{\text{HOS}} \right\rbrack \right)\left( 1 - r_{y,a} \right)}
\end{aligned}$$

Natural mortality is specified by age class. Accordingly, this mortality
corresponds to either the freshwater or marine survival depending on age
class.

### Terminal fishery

Assuming no mark-selective fishing, the retained catch of the terminal
($\text{T}$) fishery is calculated similarly as with the pre-terminal
fishery:

$$\begin{aligned}
K_{y,a}^{\text{NOS,T}} & {= \left( 1 - \exp\left( - v_{a}^{\text{T}}F_{y}^{\text{T}} \right) \right)\text{NOR}_{y,a}} \\
K_{y,a}^{\text{HOS,T}} & {= \left( 1 - \exp\left( - v_{a}^{\text{T}}F_{y}^{\text{T}} \right) \right)\text{HOR}_{y,a}} \\
 & 
\end{aligned}$$

If management is by target harvest rate, then the specified harvest rate
$u^{\text{T}}$ in the terminal ($\text{PT}$) fishery is converted to the
apical instantaneous fishing mortality rates as
$F^{\text{T}} = - \log\left( 1 - u^{\text{T}} \right)$.

If management is by target catch $K^{\text{T,target}}$, then the fishing
mortality is solved such that
$K^{\text{T,target}} = \sum_{a}K_{y,a}^{\text{NOS,T}} + \sum_{a}K_{y,a}^{\text{HOS,T}}$.

Alternative equations are used if [mark-selective
fishing](#mark-selective-fishing) is implemented.

## Escapement and spawners

The escapement consists of the survivors of the terminal fishery:

$$\begin{aligned}
\text{NOR}_{y,a}^{\text{escapement}} & {= \text{NOR}_{y,a}\exp\left( - v_{a}F_{y}^{\text{T}} \right)} \\
\text{HOR}_{y,a}^{\text{escapement}} & {= \text{HOR}_{y,a}\exp\left( - v_{a}F_{y}^{\text{T}} \right)}
\end{aligned}$$

The abundance of natural origin spawners consists of the escapement that
survive migration to the spawning ground ($s_{\text{enroute}}$) and are
not removed for brood:

$$\text{NOS}_{y,a} = \text{NOR}_{y,a}^{\text{escapement}} \times s_{\text{enroute}} - \text{NOB}_{y,a}$$

The hatchery origin spawners is the escapement of local origin that
survive migration, do not return to the hatchery (either by swim-in
facilities or in-river collection), and are not removed from the
spawning ground (through proportion $p_{\text{removal}}^{\text{HOS}}$
and discounted by the mark rate, these animals are not available for
brood). Strays not used for brood are also included as hatchery
spawners.

$$\begin{aligned}
\text{HOS}_{y,a} & {= \text{HOS}_{y,a}^{\text{local}} + \text{HOS}_{y,a}^{\text{stray}}} \\
 & {= \left( 1 - p^{\text{hatchery}} \right)\left( 1 - p_{\text{removal}}^{\text{HOS}} \times m \right)\text{HOR}_{y,a}^{\text{escapement}} \times s_{\text{enroute}} + \left( \text{Stray}_{y,a} - \text{HOB}_{y,a}^{\text{stray}} \right)}
\end{aligned}$$

## Fitness effects on survival

Reproductive success of first generation hatchery fish has been observed
to be lower than their natural counterparts, and is accounted for in the
$\gamma$ parameter (see review in [Withler et
al. 2018](https://www.dfo-mpo.gc.ca/csas-sccs/Publications/ResDocs-DocRech/2018/2018_019-eng.html)).

Through genetic and epigenetic factors, survival of hatchery juveniles
in the hatchery environment selects for fish with a phenotype best
adapted for that environment, and likewise for juveniles spawned in the
natural environment. Since these traits are heritable, the fitness of
the natural population can shift away from the optimum for the natural
environment towards that of the hatchery environment on an evolutionary
time scale, i.e., over a number of generations, when hatchery fish are
allowed to spawn.

As described in [Ford
2002](https://doi.org/10.1046/j.1523-1739.2002.00257.x) and derived in
[Lande 1976](https://doi.org/10.1111/j.1558-5646.1976.tb00911.x), the
fitness loss function $W$ for an individual with phenotypic trait value
$z$ in a given environment is

$$W(z) = \exp\left( \frac{- (z - \theta)^{2}}{2\omega^{2}} \right)$$

where $\theta$ is the optimum for that environment and $\omega^{2}$ is
the fitness variance.

If the phenotypic trait value $z$ in the population is a random normal
variable with mean $\bar{z}$ and variance $\sigma^{2}$, then the mean
fitness of the population in generation $g$ is
$\bar{W}(z) = \int W(z)f(z)dz$, where $f(z)$ is the Gaussian probability
density function. The solution is proportional to

$$\bar{W}(z) \propto \exp\left( \frac{- \left( \bar{z} - \theta \right)^{2}}{2\left( \omega^{2} + \sigma^{2} \right)} \right)$$

The mean phenotype $\bar{z}$ is calculated iteratively, where the change
$\Delta\bar{z}$ from generation $g - 1$ to $g$ is

$$\begin{aligned}
{\Delta\bar{z}} & {= {\bar{z}}_{g} - {\bar{z}}_{g - 1} = \left( {\bar{z}}_{g - 1}^{\prime} - {\bar{z}}_{g - 1} \right)h^{2}} \\
{\bar{z}}_{g} & {= {\bar{z}}_{g - 1} + \left( {\bar{z}}_{g - 1}^{\prime} - {\bar{z}}_{g - 1} \right)h^{2}} \\
 & 
\end{aligned}$$

where $h^{2}$ is the heritability of $z$ and
${\bar{z}}_{g - 1}^{\prime}$ is the trait value after applying the
fitness function, defined as:

$$\begin{aligned}
{\bar{z}}_{g - 1}^{\prime} & {= \frac{1}{{\bar{W}}_{g - 1}}\int W_{g - 1}(z) \times zf(z)dz} \\
 & {= \frac{{\bar{z}}_{g - 1}\omega^{2} + \theta\sigma^{2}}{\omega^{2} + \sigma^{2}}}
\end{aligned}$$

Let ${\bar{z}}_{g - 1}^{\prime}(\theta)$ be a function that returns the
mean trait value after selection in an environment with optimum value
$\theta$. With a hatchery program, the mean trait value of the progeny
in the natural environment is a weighted average of the mean trait value
in natural and hatchery origin spawners, with selection in the natural
environment, i.e., with optimum trait value $\theta^{\text{natural}}$:

$$\begin{aligned}
{{\bar{z}}_{g}^{\text{natural}} =} & {\left( 1 - p_{g - 1}^{\text{HOSeff}} \right) \times \left( {\bar{z}}_{g - 1}^{\text{natural}} + \left\lbrack {\bar{z}}_{g - 1}^{\prime\text{natural}}\left( \theta^{\text{natural}} \right) - {\bar{z}}_{g - 1}^{\text{natural}} \right\rbrack h^{2} \right) +} \\
 & {p_{g - 1}^{\text{HOSeff}} \times \left( {\bar{z}}_{g - 1}^{\text{hatchery}} + \left\lbrack {\bar{z}}_{g - 1}^{\prime\text{hatchery}}\left( \theta^{\text{natural}} \right) - {\bar{z}}_{g - 1}^{\text{hatchery}} \right\rbrack h^{2} \right)}
\end{aligned}$$

where
$p^{\text{HOSeff}} = \gamma \times \text{HOS}/\left( \text{NOS} + \gamma \times \text{HOS} \right)$.

Similarly, the mean trait value in the hatchery environment
${\bar{z}}_{g}^{\text{hatchery}}$ is a weighted average of the mean
trait value of the natural and hatchery broodtake, with selection in the
hatchery environment, i.e., with optimum trait value
$\theta^{\text{hatchery}}$:

$$\begin{aligned}
{{\bar{z}}_{g}^{\text{hatchery}} =} & {p_{g - 1}^{\text{NOB}} \times \left( {\bar{z}}_{g - 1}^{\text{natural}} + \left\lbrack {\bar{z}}_{g - 1}^{\prime\text{natural}}\left( \theta^{\text{hatchery}} \right) - {\bar{z}}_{g - 1}^{\text{natural}} \right\rbrack h^{2} \right) +} \\
 & {\left( 1 - p_{g - 1}^{\text{NOB}} \right) \times \left( {\bar{z}}_{g - 1}^{\text{hatchery}} + \left\lbrack {\bar{z}}_{g - 1}^{\prime\text{hatchery}}\left( \theta^{\text{hatchery}} \right) - {\bar{z}}_{g - 1}^{\text{hatchery}} \right\rbrack h^{2} \right)}
\end{aligned}$$

where
$p^{\text{NOB}} = \text{NOB}/\left( \text{NOB} + \text{HOB} \right)$.

The fitness variance $\omega^{2}$ and phenotype variance $\sigma^{2}$
are identical in the two environments.

The mean fitness of generation $g$ in the natural environment is then:

$${\bar{W}}_{g}^{\text{natural}} = \exp\left( \frac{- \left( {\bar{z}}_{g}^{\text{natural}} - \theta^{\text{natural}} \right)^{2}}{2\left( \omega^{2} + \sigma^{2} \right)} \right)$$

### Mixed brood-year return

If a mixed-brood year return in year $y$ across several ages $a$
produces the smolt cohort in year $y + 1$, then the mean trait value in
the progeny is calculated from a weighted average by brood year and age
class fecundity:

$$\begin{aligned}
{{\bar{z}}_{y + 1}^{\text{natural}} =} & {\sum\limits_{a}p_{y,a}^{\text{NOS}} \times \left( {\bar{z}}_{y - a + 1}^{\text{natural}} + \left\lbrack {\bar{z}}_{y - a + 1}^{\prime\text{natural}}\left( \theta^{\text{natural}} \right) - {\bar{z}}_{y - a + 1}^{\text{natural}} \right\rbrack h^{2} \right) +} \\
 & {\sum\limits_{a}p_{y,a}^{\text{HOSeff}} \times \left( {\bar{z}}_{y - a + 1}^{\text{hatchery}} + \left\lbrack {\bar{z}}_{y - a + 1}^{\prime\text{hatchery}}\left( \theta^{\text{natural}} \right) - {\bar{z}}_{y - a + 1}^{\text{hatchery}} \right\rbrack h^{2} \right)}
\end{aligned}$$

$$\begin{aligned}
{{\bar{z}}_{y + 1}^{\text{hatchery}} =} & {\sum\limits_{a}p_{y,a}^{\text{NOB}} \times \left( {\bar{z}}_{y - a + 1}^{\text{natural}} + \left\lbrack {\bar{z}}_{y - a + 1}^{\prime\text{natural}}\left( \theta^{\text{hatchery}} \right) - {\bar{z}}_{y - a + 1}^{\text{natural}} \right\rbrack h^{2} \right) +} \\
 & {\sum\limits_{a}p_{y,a}^{\text{HOB}} \times \left( {\bar{z}}_{y - a + 1}^{\text{hatchery}} + \left\lbrack {\bar{z}}_{y - a + 1}^{\prime\text{hatchery}}\left( \theta^{\text{hatchery}} \right) - {\bar{z}}_{y - a + 1}^{\text{hatchery}} \right\rbrack h^{2} \right)}
\end{aligned}$$

where

$p_{y,a}^{\text{NOS}} = \frac{\text{Fec}_{a} \times \text{NOS}_{y,a}}{\sum_{a}\text{Fec}_{a}\left( \text{NOS}_{y,a} + \gamma \times \text{HOS}_{y,a} \right)}$

$p_{y,a}^{\text{HOSeff}} = \frac{\text{Fec}_{a} \times \gamma \times \text{HOS}_{y,a}}{\sum_{a}\text{Fec}_{a}\left( \text{NOS}_{y,a} + \gamma \times \text{HOS}_{y,a} \right)}$

$p_{y,a}^{\text{NOB}} = \frac{\text{Fec}_{a}^{\text{brood}} \times \text{NOB}_{y,a}}{\sum_{a}\text{Fec}_{a}^{\text{brood}}\left( \text{NOB}_{y,a} + \text{HOB}_{y,a} + \text{Brood}_{y,a}^{\text{import}} \right)}$

$p_{y,a}^{\text{HOB}} = \frac{\text{Fec}_{a}^{\text{brood}} \times \left( \text{HOB}_{y,a} + \text{Brood}_{y,a}^{\text{import}} \right)}{\sum_{a}\text{Fec}_{a}^{\text{brood}}\left( \text{NOB}_{y,a} + \text{HOB}_{y,a} + \text{Brood}_{y,a}^{\text{import}} \right)}$

Effective proportions, i.e., weighting by age-class fecundity, accounts
for older age classes that are more fecund and more likely to contribute
to the production of next cohort.

### Fitness loss

Fitness can reduce survival in the egg, fry, and immature life stages.

If no habitat model is used, then the egg-fry survival is reduced by the
fitness loss function:

$$\begin{aligned}
\text{Fry}_{y + 1}^{\text{NOS}} & {= \text{Egg}_{y}^{\text{NOS}} \times \left( W_{y}^{\text{nat.}} \right)^{\ell_{\text{egg}}}} \\
\text{Fry}_{y + 1}^{\text{HOS}} & {= \text{Egg}_{y}^{\text{HOS}} \times \left( W_{y}^{\text{nat.}} \right)^{\ell_{\text{egg}}}}
\end{aligned}$$

and the smolt production function is adjusted by loss in productivity
and capacity, with $\alpha$ and $\beta$ adjusted accordingly as:

$$\begin{aligned}
\text{Smolt}_{y + 1}^{\text{NOS}} & {= \frac{\alpha\prime_{y + 1} \times \text{Fry}_{y + 1}^{\text{NOS}}}{1 + \beta\prime_{y + 1}\left( \text{Fry}_{y + 1}^{\text{NOS}} + \text{Fry}_{y + 1}^{\text{HOS}} + n_{y + 1}^{\text{sub}} \right)}} \\
\text{Smolt}_{y + 1}^{\text{HOS}} & {= \frac{\alpha\prime_{y + 1} \times \text{Fry}_{y + 1}^{\text{HOS}}}{1 + \beta\prime_{y + 1}\left( \text{Fry}_{y + 1}^{\text{NOS}} + \text{Fry}_{y + 1}^{\text{HOS}} + n_{y + 1}^{\text{sub}} \right)}}
\end{aligned}$$

with
$\alpha\prime_{y + 1} = \left( W_{y}^{\text{nat.}} \right)^{\ell_{\text{fry}}} \times \kappa/\phi$
and
$\beta\prime_{y + 1} = \alpha/\left( C_{\text{egg-smolt}} \times \left( W_{y}^{\text{nat.}} \right)^{\ell_{\text{fry}}} \right)$.

> With the Ricker density-dependent survival, the beta parameter is
> adjusted with
> $\beta_{y}^{*} = 1/\left\lbrack S_{\text{max}} \times \left( W_{y}^{\text{nat.}} \right)^{\ell_{\text{fry}}} \right\rbrack$.

In the marine life stage, the increase in natural mortality is:

$$M_{y,a}^{\text{NOS}} = - \log\left( \exp\left( - M_{y,a}^{\text{base,NOS}} \right) \times \left( W_{y - a}^{\text{nat.}} \right)^{\ell_{\text{juv}}} \right)$$

In the marine environment, age-specific natural survival is reduced
proportional to the fitness loss term and modeled as a cohort effect.

Parameter $\ell_{i}$ is the proportion of the fitness loss apportioned
to life stage $i$ (either egg, fry, or juvenile-marine), with
$\sum_{i}\ell_{i} = 1$.

If habitat variables are modeled, then the egg and fry fitness losses
adjust the productivity and capacity of the corresponding life stage:

$$\begin{aligned}
P_{y}^{\text{egg-fry}} & {= P^{\text{egg-fry}} \times \left( W_{y}^{\text{nat.}} \right)^{\ell_{\text{egg}}}} \\
P_{y}^{\text{fry-smolt}} & {= P^{\text{fry-smolt}} \times \left( W_{y}^{\text{nat.}} \right)^{\ell_{\text{fry}}}}
\end{aligned}$$

### PNI

PNI (proportionate natural influence) is an approximation of the rate of
gene flow from the hatchery to the natural environment, calculated for
the progeny in year $y + 1$ from the parental composition of year $y$:

$$\text{PNI}_{y + 1} = \frac{\sum\limits_{a}p_{y,a}^{\text{NOB}}}{\sum\limits_{a}p_{y,a}^{\text{NOB}} + \sum\limits_{a}p_{y,a}^{\text{HOSeff}}}$$

Generally, a combination of minimizing hatchery releases, increasing
natural broodtake, and reducing the number hatchery origin spawners
maintains high PNI, i.e., low rate of gene flow from the hatchery to
natural environment.

> If there is no natural origin broodtake, i.e., all brood is imported,
> then PNI is calculated with equation 6 of [Withler et
> al. 2018](https://www.dfo-mpo.gc.ca/csas-sccs/Publications/ResDocs-DocRech/2018/2018_019-eng.html):
>
> $$\text{PNI}_{y + 1} = \frac{h^{2}}{h^{2} + \left( 1 - h^{2} + \omega^{2} \right)\sum\limits_{a}p_{y,a}^{\text{HOSeff}}}$$

### Wild salmon

With single brood-year returns, the proportion of wild salmon, natural
origin spawners whose parents were also natural spawners, can be
calculated as

$$p_{g}^{\text{WILD}} = \left( 1 - p_{g}^{\text{HOScensus}} \right) \times \frac{\left( 1 - p_{g - 1}^{\text{HOScensus}} \right)^{2}}{\left( 1 - p_{g - 1}^{\text{HOScensus}} \right)^{2} + 2\gamma \times p_{g - 1}^{\text{HOScensus}}\left( 1 - p_{g - 1}^{\text{HOScensus}} \right) + \gamma^{2}\left( p_{g - 1}^{\text{HOScensus}} \right)^{2}}$$

where
$p^{\text{HOScensus}} = \text{HOS}/\left( \text{HOS} + \text{NOS} \right)$.

The first term is the proportion of natural spawners in the current
generation $g$.

The ratio comprising the second term discounts the proportion of the
current generation to include natural spawners whose parents were both
natural spawners. Assuming non-assortative mating, the three terms in
the denominator gives the composition of generation $g$ whose parents
who are both natural origin, mixed origin (one parent in natural origin
and the other is hatchery origin), and both hatchery origin.

To generalize for mixed-brood year return, we calculate the probability
weighted across brood-years and age class fecundity:

$$p_{y}^{\text{WILD}} = \sum\limits_{a}\frac{\text{NOS}_{y,a}}{\sum\limits_{a\prime}\left( \text{NOS}_{y,a\prime} + \text{HOS}_{y,a\prime} \right)} \times \frac{\left( \sum\limits_{a\prime}p_{y - a,a\prime}^{\text{NOScensus}} \right)^{2}}{\left( \sum\limits_{a\prime}p_{y - a,a\prime}^{\text{NOScensus}} \right)^{2} + 2\gamma \times \left( \sum\limits_{a\prime}p_{y - a,a\prime}^{\text{NOScensus}} \right)\left( \sum\limits_{a\prime}p_{y - a,a\prime}^{\text{HOScensus}} \right) + \gamma^{2}\left( \sum\limits_{a\prime}p_{y - a,a\prime}^{\text{HOScensus}} \right)^{2}}$$

where

$p_{y,a}^{\text{NOScensus}} = \frac{\text{Fec}_{a} \times \text{NOS}_{y,a}}{\sum_{a}{\text{Fec}_{a}(\text{NOS}_{y,a}} + \text{HOS}_{y,a})}$

$p_{y,a}^{\text{HOScensus}} = \frac{\text{Fec}_{a} \times \text{HOS}_{y,a}}{\sum_{a}{\text{Fec}_{a}(\text{NOS}_{y,a}} + \text{HOS}_{y,a})}$

The probability of finding a wild salmon in year $y$ is the sum of
probabilities of finding a wild salmon over all ages. For each age $a$,
the first ratio is the probability of finding a natural spawner in year
$y$. The second ratio is the probability of mating success from two
parental natural spawners in year $y - a$ using a Punnett square,
assuming non-assortative mating across age and origin. The summation
across dummy age variable $a\prime$ calculates the total proportion of
spawners in a given year.

Effective proportions, i.e., weighting by age-class fecundity, in the
parental composition accounts for older age classes that are more fecund
and more likely to contribute to the production of offspring.

## Mark-selective fishing

If the mark rate $m$ of hatchery fish is greater than zero, then
mark-selective fishing can be implemented for both the pre-terminal and
terminal fisheries with no retention on natural-origin fish. The mark
rate is a proxy for retention and the harvest rate $u^{\text{harvest}}$
corresponds to the apical fishing mortality of the kept catch. The
exploitation rate $u^{\text{exploit}}$ is calculated from kept catch of
hatchery-origin fish and dead releases of natural-origin fish.
Exploitation rates differ between hatchery and natural origin fish
because there is no retention of the latter.

Fishing mortality on hatchery origin fish is partitioned into two parts,
for kept and release catch:

$$\begin{aligned}
F^{\text{kept,HO}} & {= mE} \\
F^{\text{rel,HO}} & {= (1 - m)\delta E}
\end{aligned}$$

where $\delta$ is the proportion of released fish that die, i.e.,
release mortality.

$E$ is an index of fishing effort, also referred to as the encounter
rate by the fishery, that links together $F^{\text{kept,HO}}$ and
$F^{\text{rel,NO}}$.

Natural-origin fish experience exclusively fishing mortality due to
deaths from release:

$$F^{\text{rel,NO}} = \delta E$$

Intuitively, fishing effort can increase in a mark-selective fishery
compared to a non-selective fishery. For example, if the mark rate is 20
percent, then the fishing effort could be 500 percent higher than in a
non-selective fishery in order to attain the kept quota or bag limit.
Additional catch and release mortality then occurs for un-marked fish,
according to $\delta$.

The kept catch $K$ are

$$\begin{aligned}
K_{y,a}^{\text{HOS,PT}} & {= \frac{F_{y}^{\text{kept,HO,PT}}}{F_{y}^{\text{kept,HO,PT}} + F_{y}^{\text{rel,HO,PT}}}\left( 1 - \exp\left( - v_{a}^{\text{PT}}\left\lbrack F^{\text{kept,HO,PT}} + F^{\text{rel,HO,PT}} \right\rbrack \right) \right)N_{y,a}^{\text{juv,HOS}}} \\
K_{y,a}^{\text{HOS,T}} & {= \frac{F_{y}^{\text{kept,HO,T}}}{F_{y}^{\text{kept,HO,T}} + F_{y}^{\text{rel,HO,T}}}\left( 1 - \exp\left( - v_{a}^{\text{T}}\left\lbrack F^{\text{kept,HO,T}} + F^{\text{rel,HO,T}} \right\rbrack \right) \right)\text{HOR}_{y,a}}
\end{aligned}$$

If management is by target harvest rate for the preterminal and terminal
fisheries, the corresponding $E$ is solved to satisfy:

$$\begin{aligned}
u^{\text{PT}} & {= 1 - \exp\left( - F^{\text{kept,HO,PT}} \right)} \\
u^{\text{T}} & {= 1 - \exp\left( - F^{\text{kept,HO,T}} \right)}
\end{aligned}$$

If management is by target catch, the corresponding $E$ is solved to
satisfy:

$$\begin{aligned}
K^{\text{PT,target}} & {= \sum\limits_{a}K_{y,a}^{\text{HOS,PT}}} \\
K^{\text{T,target}} & {= \sum\limits_{a}K_{y,a}^{\text{HOS,T}}}
\end{aligned}$$

The exploitation rate for natural origin fish is calculated from dead
discards. The exploitation rate for hatchery origin fish is calculated
from kept catch and dead discards:

$$\begin{aligned}
u_{y,a}^{\text{exploit,NOS,PT}} & {= 1 - \exp\left( - v_{a}F_{y}^{\text{rel.,NO,PT}} \right)} \\
u_{y,a}^{\text{exploit,HOS,PT}} & {= 1 - \exp\left( - v_{a}\left\lbrack F_{y}^{\text{kept,HO,PT}} + F_{y}^{\text{rel.,HO,PT}} \right\rbrack \right)}
\end{aligned}$$

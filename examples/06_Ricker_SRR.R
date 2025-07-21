
# Compare Ricker reference points from salmonMSE::calc_ref(), using numerical optimization
# with those using the Lambert_W function (Scheuerell 2016)

#### Define functions ----
# install.packages("gsl")
# remotes::install_github("Blue-Matter/salmonMSE")
source("https://raw.githubusercontent.com/Pacific-salmon-assess/samEst/refs/heads/main/R/RP_functions.R")


compare_ref <- function(alpha = 3, # Units of recruits/spawner
                        Smax = 1000, # Units of spawners
                        maxage = 5,
                        rel_F = c(0, 1),
                        p_female = 1,
                        vulPT = rep(0, maxage),
                        vulT = c(0, 0, 0.2, 0.5, 1),
                        Mjuv = c(1, 0.8, 0.6, 0.4, 0.2),
                        fec = c(0, 1500, 3000, 3200, 3500),
                        s_enroute = 1,
                        p_mature = c(0, 0.1, 0.2, 0.3, 1)) {

  # Calculate eggs per smolt
  phi <- salmonMSE:::calc_phi(
    Mjuv = Mjuv,
    p_mature = p_mature,
    p_female = p_female,
    fec = matrix(fec, maxage, 1),
    s_enroute = s_enroute,
    n_g = 1,
    p_LHG = 1
  )

  # To convert Smax from units of spawners to eggs:
  # 1. Calculate unfished spawners (from Ricker SRR, set R = S --> S0 = exp(1/Smax)/alpha)
  # 2. Calculate spawners per smolt (spro)
  # 3. Calculate egg per smolt (phi)
  # 4. unfished eggs (eo) = unfished spawners * smolt per spawner * egg per smolt
  # 5. Smax_eggs = log(alpha) / unfished eggs
  spro <- salmonMSE:::calc_phi(
    Mjuv = Mjuv,
    p_mature = p_mature,
    p_female = p_female,
    fec = matrix(fec, maxage, 1),
    s_enroute = s_enroute,
    n_g = 1,
    p_LHG = 1,
    output = "spawner"
  )

  so <- Smax * log(alpha)
  eo <- so / spro * phi
  beta_eggs <- log(alpha)/eo
  Smax_eggs <- 1/beta_eggs

  SRRpars <- data.frame(
    kappa = alpha,
    Smax = Smax_eggs,
    phi = phi,
    SRrel = "Ricker"
  )

  ref <- salmonMSE:::calc_MSY(
    Mjuv, fec, p_female, rel_F, vulPT, vulT, p_mature,
    s_enroute, SRRpars = SRRpars
  )

  ref["Sgen"] <- salmonMSE:::calc_Sgen(
    Mjuv, fec, p_female, rel_F, vulPT, vulT, p_mature,
    s_enroute, SRRpars = SRRpars, SMSY = ref["Spawners_MSY"]
  )

  ref_salmonMSE <- ref[c("UPT_MSY", "UT_MSY", "Spawners_MSY", "Sgen")] |>
    structure(names = c("UMSY (preterminal)", "UMSY (terminal)", "SMSY", "Sgen"))

  umsy <- umsyCalc(log(alpha))
  smsy <- smsyCalc(log(alpha), 1/Smax)
  sgen <- sgenCalcDirect(log(alpha), 1/Smax)

  ref_lambert <- structure(c(umsy, smsy, sgen), names = c("UMSY", "SMSY", "Sgen"))

  list(salmonMSE = ref_salmonMSE, lambert = ref_lambert)

}


#### Simple life cycle ----
maxage <- 5
compare_ref(
  alpha = 3,
  Smax = 1000,
  maxage = 5,
  rel_F = c(0, 1),
  p_female = 1,
  vulPT = rep(0, maxage),
  vulT = c(0, 0, 0, 0, 1),
  Mjuv = c(-log(0.01), 0, 0, 0, 0), # SAR = 0.01
  fec = c(0, 0, 0, 0, 1),
  s_enroute = 1,
  p_mature = c(0, 0, 0, 0, 1)
)

#### Age-specific M, selectivity, maturity ----
maxage <- 5
compare_ref(
  alpha = 3,
  Smax = 1000,
  maxage = 5,
  rel_F = c(0, 1),
  p_female = 1,
  vulPT = rep(0, maxage),
  vulT = c(0, 0.1, 0.2, 0.4, 1),
  Mjuv = c(1, 0.1, 0.1, 0.1, 0.1), # SAR = 0.01
  fec = c(0, 1000, 2000, 3000, 3500),
  s_enroute = 1,
  p_mature = c(0, 0.1, 0.2, 0.3, 1)
)


#### Age-specific M, selectivity, maturity with all exploitation in pre-terminal fishery instead of terminal ----
maxage <- 5
compare_ref(
  alpha = 3,
  Smax = 1000,
  maxage = 5,
  rel_F = c(1, 0),
  p_female = 1,
  vulPT = c(0, 0.1, 0.2, 0.4, 1),
  vulT = rep(0, maxage),
  Mjuv = c(1, 0.1, 0.1, 0.1, 0.1), # SAR = 0.01
  fec = c(0, 1000, 2000, 3000, 3500),
  s_enroute = 1,
  p_mature = c(0, 0.1, 0.2, 0.3, 1)
)



# Compare Ricker reference points from salmonMSE::calc_ref(), using numerical optimization
# with those using the Lambert_W function (Scheuerell 2016)

#### Define functions ----
# install.packages("gsl")
# remotes::install_github("Blue-Matter/salmonMSE", ref = "dev")
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
                        p_mature = c(0, 0.1, 0.2, 0.3, 1),
                        maximize = c("MER", "MSY")) {

  maximize <- match.arg(maximize)

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
    fec = fec,
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

  # MER
  ref <- salmonMSE:::calc_MSY(
    Mjuv, fec, p_female, rel_F, vulPT, vulT, p_mature,
    s_enroute, SRRpars = SRRpars, maximize = "MER"
  )

  ref["Sgen"] <- salmonMSE:::calc_Sgen(
    Mjuv, fec, p_female, rel_F, vulPT, vulT, p_mature,
    s_enroute, SRRpars = SRRpars, SMSY = ref["Spawners_MSY"]
  )

  ref["Catch/Return"] <- ref["KT_MSY"]/ref["Return_MSY"]

  ref_salmonMSE <- ref[c("UPT_MSY", "UT_MSY", "Catch/Return", "Spawners_MSY", "Sgen")] |>
    structure(names = c("UMSY (preterminal)", "UMSY (terminal)", "Catch/Return", "SMSY", "Sgen"))

  # MSY
  ref_MSY <- salmonMSE:::calc_MSY(
    Mjuv, fec, p_female, rel_F, vulPT, vulT, p_mature,
    s_enroute, SRRpars = SRRpars, maximize = "MSY"
  )

  ref_MSY["Sgen"] <- salmonMSE:::calc_Sgen(
    Mjuv, fec, p_female, rel_F, vulPT, vulT, p_mature,
    s_enroute, SRRpars = SRRpars, SMSY = ref_MSY["Spawners_MSY"]
  )

  ref_MSY["Catch/Return"] <- ref_MSY["KT_MSY"]/ref_MSY["Return_MSY"]

  MSY_salmonMSE <- ref_MSY[c("UPT_MSY", "UT_MSY", "Catch/Return", "Spawners_MSY", "Sgen")] |>
    structure(names = c("UMSY (preterminal)", "UMSY (terminal)", "Catch/Return", "SMSY", "Sgen"))

  # Lambert naive calculations
  ref_lambert <- local({
    umsy <- umsyCalc(log(alpha))
    smsy <- smsyCalc(log(alpha), 1/Smax)
    sgen <- sgenCalcDirect(log(alpha), 1/Smax)
    structure(c(umsy, smsy, sgen), names = c("UMSY", "SMSY", "Sgen"))
  })

  # Lambert adjusted calculations
  juv_surv <- salmonMSE:::calc_survival(Mjuv, p_mature)
  phi_return <- sum(juv_surv * p_mature) # Recruit/smolt

  alpha2 <- alpha/phi # smolt/egg
  alpha_prime <- alpha2 * phi_return # recruit/egg
  #alpha_prime <- alpha2 * phi_return * eo/so# This just returns the same numbers (alpha = alpha_prime)
  if (alpha_prime < 1) warning("alpha_prime < 1")

  ref_lambert2 <- local({
    umsy <- umsyCalc(log(alpha_prime))
    smsy <- smsyCalc(log(alpha_prime), 1/Smax)
    sgen <- sgenCalcDirect(log(alpha_prime), 1/Smax)
    structure(c(umsy, smsy, sgen), names = c("UMSY", "SMSY", "Sgen"))
  })

  list(`salmonMSE MER` = ref_salmonMSE,
       `salmonMSE MSY` = MSY_salmonMSE,
       `lambert adjusted` = ref_lambert2,
       `lambert naive` = ref_lambert)

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
  Mjuv = c(1, 0.1, 0.1, 0.1, 0.1),
  fec = c(0, 1000, 2000, 3000, 3500),
  s_enroute = 1,
  p_mature = c(0, 0.1, 0.2, 0.3, 1)
)

#### Age-specific M, maturity but full vulnerability ----
maxage <- 5
compare_ref(
  alpha = 3,
  Smax = 1000,
  maxage = 5,
  rel_F = c(0, 1),
  p_female = 1,
  vulPT = rep(0, maxage),
  vulT = rep(1, maxage),
  Mjuv = c(1, 0.1, 0.1, 0.1, 0.1),
  fec = c(0, 1000, 2000, 3000, 3500),
  s_enroute = 1,
  p_mature = c(0, 0.1, 0.2, 0.3, 1)
)




#### Age-specific M, selectivity, maturity with all exploitation in pre-terminal fishery instead of terminal (for actual yield curve) ----
#### MER requires an assumption about terminal fishery vulnerability
maxage <- 5
compare_ref(
  alpha = 3,
  Smax = 1000,
  maxage = 5,
  rel_F = c(1, 0),
  p_female = 1,
  vulPT = c(0, 0.1, 0.2, 0.4, 1),
  vulT = c(0, 0.1, 0.2, 0.4, 1), # Needs some assumption for MER calculation
  Mjuv = c(1, 0.1, 0.1, 0.1, 0.1), # SAR = 0.01
  fec = c(0, 1000, 2000, 3000, 3500),
  s_enroute = 1,
  p_mature = c(0, 0.1, 0.2, 0.3, 1)
)


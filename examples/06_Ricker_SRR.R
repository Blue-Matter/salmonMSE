
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
                        maximize = c("MER", "MSY"),
                        Sgen_nyears = maxage,
                        Sgen_output = FALSE) {

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
    s_enroute, SRRpars = SRRpars, SMSY = ref["Spawners_MSY"],
    nyears = Sgen_nyears
  )

  ref["Catch/Return"] <- ref["KT_MSY"]/ref["Return_MSY"]

  ref_salmonMSE <- ref[c("UPT_MSY", "UT_MSY", "Catch/Return", "Spawners_MSY", "Sgen")] |>
    structure(names = c("UMSY (preterminal)", "UMSY (terminal)", "Catch/Return", "SMSY", "Sgen"))

  # MSY
  ref_MSY <- salmonMSE:::calc_MSY(
    Mjuv, fec, p_female, rel_F, vulPT, vulT, p_mature,
    s_enroute, SRRpars = SRRpars, maximize = "MSY"
  )
  Sgen <- salmonMSE:::calc_Sgen(
    Mjuv, fec, p_female, rel_F, vulPT, vulT, p_mature,
    s_enroute, SRRpars = SRRpars, SMSY = ref_MSY["Spawners_MSY"],
    nyears = Sgen_nyears
  )

  ref_MSY["Sgen"] <- Sgen

  ref_MSY["Catch/Return"] <- ref_MSY["KT_MSY"]/ref_MSY["Return_MSY"]

  MSY_salmonMSE <- ref_MSY[c("UPT_MSY", "UT_MSY", "Catch/Return", "Spawners_MSY", "Sgen")] |>
    structure(names = c("UMSY (preterminal)", "UMSY (terminal)", "Catch/Return", "SMSY", "Sgen"))

  # Fsearch
  #FM <- seq(0, 3, 0.005)
  #Fsearch <- sapply(FM, function(i) {
  #  salmonMSE:::.calc_eq(
  #    .F = i,
  #    matrix(Mjuv, maxage, 1),
  #    fec, p_female = 1, rel_F = c(0, 1),
  #    vulPT = vulPT,
  #    vulT = vulT,
  #    p_mature = matrix(p_mature, maxage, 1),
  #    s_enroute = s_enroute,
  #    n_g = 1, p_LHG = 1,
  #    SRRpars = SRRpars, opt = FALSE
  #  ) %>% unlist()
  #})
  #i <- Fsearch["KT", ] >= 0
  #Fsearch <- Fsearch[, i]
  #FM <- FM[i]
  #maxF <- length(FM)
  #Fsearch["Return", maxF]/Fsearch["Spawners", maxF] # Equals to alpha

  # Lambert naive calculations
  ref_lambert <- local({
    umsy <- umsyCalc(log(alpha))
    smsy <- smsyCalc(log(alpha), 1/Smax)
    sgen <- sgenCalcDirect(log(alpha), 1/Smax)
    structure(c(umsy, smsy, sgen), names = c("UMSY", "SMSY", "Sgen"))
  })

  # Lambert adjusted calculations
  #juv_surv <- salmonMSE:::calc_survival(Mjuv, p_mature)
  #phi_return <- sum(juv_surv * p_mature) # Recruit/smolt
  #alpha2 <- alpha/phi # smolt/egg
  #alpha_prime <- alpha2 * phi_return # recruit/egg
  #alpha_prime <- alpha2 * phi_return * eo/so# This just returns the same numbers (alpha = alpha_prime)
  #if (alpha_prime < 1) warning("alpha_prime < 1")
  #ref_lambert2 <- local({
  #  umsy <- umsyCalc(log(alpha_prime))
  #  smsy <- smsyCalc(log(alpha_prime), 1/Smax)
  #  sgen <- sgenCalcDirect(log(alpha_prime), 1/Smax)
  #  structure(c(umsy, smsy, sgen), names = c("UMSY", "SMSY", "Sgen"))
  #})

  output <- list(
    `salmonMSE MER` = ref_salmonMSE,
    `salmonMSE MSY` = MSY_salmonMSE,
    #`lambert adjusted` = ref_lambert2,
    `lambert naive` = ref_lambert
  )
  if (Sgen_output) {
    attr(output, "Sgen") <- attr(Sgen, "Sgen")
    attr(output, "Sgen_proj") <- attr(Sgen, "proj")
  }

  output

}


#### Simple life cycle ----
maxage <- 5
output <- compare_ref(
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
  p_mature = c(0, 0, 0, 0, 1),
  Sgen_output = TRUE,
  Sgen_nyears = 1 # Try Sgen with 1 year projection, this is closer to the lambert calculation
)
output[1:3]
attr(output, "Sgen")
attr(output, "Sgen_proj")

#### Age-specific M and maturity ----
# There's no difference in spawners among ages
# No difference between Lambert and numerical MSY, difference in Sgen
maxage <- 5
output <- compare_ref(
  alpha = 3,
  Smax = 1000,
  maxage = 5,
  rel_F = c(0, 1),
  p_female = 1,
  vulPT = rep(0, maxage),
  vulT = rep(1, maxage),
  Mjuv = c(1, 0.3, 0.2, 0.1, 0.1),
  fec = rep(1, maxage),
  s_enroute = 1,
  p_mature = c(0, 0.1, 0.2, 0.3, 1),
  Sgen_output = TRUE,
  Sgen_nyears = 1 # Try Sgen with 1 year projection, this is closer to the lambert calculation
)
output[1:3]
attr(output, "Sgen")
attr(output, "Sgen_proj")

#### Age-specific M and maturity (and fecundity) ----
# Older spawners are more fecund
# No difference between Lambert and numerical MSY, difference in Sgen
maxage <- 5
compare_ref(
  alpha = 3,
  Smax = 1000,
  maxage = 5,
  rel_F = c(0, 1),
  p_female = 1,
  vulPT = rep(0, maxage),
  vulT = rep(1, maxage),
  Mjuv = c(1, 0.3, 0.2, 0.1, 0.1),
  fec = c(0, 1000, 2000, 3000, 3500),
  s_enroute = 1,
  p_mature = c(0, 0.1, 0.2, 0.3, 1),
  Sgen_nyears = 1
)

#### Age-specific M and maturity (and fecundity) and partial vulnerability of return ----
# Older spawners are more fecund
# Difference between Lambert and numerical MSY, difference in Sgen
maxage <- 5
compare_ref(
  alpha = 3,
  Smax = 1000,
  maxage = 5,
  rel_F = c(0, 1),
  p_female = 1,
  vulPT = rep(0, maxage),
  vulT = c(0, 0.1, 0.2, 0.4, 1),
  Mjuv = c(1, 0.3, 0.2, 0.1, 0.1),
  fec = c(0, 1000, 2000, 3000, 3500),
  s_enroute = 1,
  p_mature = c(0, 0.1, 0.2, 0.3, 1),
  Sgen_nyears = 1
)

# Set fecundity = 1 for all ages
maxage <- 5
compare_ref(
  alpha = 3,
  Smax = 1000,
  maxage = 5,
  rel_F = c(0, 1),
  p_female = 1,
  vulPT = rep(0, maxage),
  vulT = c(0, 0.1, 0.2, 0.4, 1),
  Mjuv = c(1, 0.3, 0.2, 0.1, 0.1),
  fec = rep(1, maxage),
  #fec = c(0, 1000, 2000, 3000, 3500),
  s_enroute = 1,
  p_mature = c(0, 0.1, 0.2, 0.3, 1),
  Sgen_nyears = 1
)

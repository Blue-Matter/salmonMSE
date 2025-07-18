
# Compare Ricker reference points from salmonMSE::calc_ref(), using numerical optimization
# with those using the Lambert_W function (Scheuerell 2016)

#### Define functions ----
# install.packages("gsl")
# remotes::install_github("Blue-Matter/salmonMSE")
source("https://raw.githubusercontent.com/Pacific-salmon-assess/samEst/refs/heads/main/R/RP_functions.R")


compare_ref <- function(alpha = 3,
                        Smax = 1000,
                        maxage = 5,
                        rel_F = c(0, 1),
                        p_female = 1,
                        vulPT = rep(0, maxage),
                        vulT = c(0, 0, 0.2, 0.5, 1),
                        Mjuv = c(1, 0.8, 0.6, 0.4, 0.2),
                        fec = c(0, 1500, 3000, 3200, 3500),
                        s_enroute = 1,
                        p_mature = c(0, 0.1, 0.2, 0.3, 1)) {

  phi <- salmonMSE:::calc_phi(
    Mjuv = Mjuv,
    p_mature = p_mature,
    p_female = p_female,
    fec = matrix(fec, maxage, 1),
    s_enroute = s_enroute,
    n_g = 1,
    p_LHG = 1
  )

  SRRpars <- data.frame(
    kappa = alpha,
    Smax = Smax,
    phi = phi,
    SRrel = "Ricker"
  )

  ref <- salmonMSE:::calc_MSY(
    matrix(Mjuv, maxage, 1), fec, p_female, rel_F, vulPT, vulT, matrix(p_mature, maxage, 1),
    s_enroute, n_g = 1, p_LHG = 1,
    SRRpars
  )

  ref["Sgen"] <- salmonMSE:::calc_Sgen(
    matrix(Mjuv, maxage, 1), fec, p_female, rel_F, vulPT, vulT, matrix(p_mature, maxage, 1),
    s_enroute, n_g = 1, p_LHG = 1,
    SRRpars, SMSY = ref["Spawners_MSY"]
  )

  ref_salmonMSE <- ref[c("UT_MSY", "Spawners_MSY", "Sgen")] |> structure(names = c("UMSY", "SMSY", "Sgen"))

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



#' Calculate abundance from density-dependent mortality
#'
#' Calculates the abundance of survivors after applying either a Beverton-Holt or Ricker stock-recruit relationship.
#'
#' @param N1 Numeric, the initial abundance that scales the density-independent survival term
#' @param N2 Numeric, the initial abundance that scales the density-dependent survival term
#' @param p Numeric, the productivity parameter that sets the maximum survival as the initial abundance approaches zero
#' @param capacity Numeric, the capacity parameter that set the maximum survivors
#' @param type Character, the functional form of the stock-recruit relationship
#' @returns
#' - `calc_SRR` calculates the abundance of survivors
#' - `calc_SRRpars` calculates the alpha and beta terms when the productivity parameter is in terms of abundance but N1 and N2 is in terms of
#'
#' @details
#' The Beverton-Holt stock recruit relationship is of the following form:
#' \deqn{\textrm{Smolt} = \dfrac{\alpha N_1}{1 + \beta N_2}}
#' where \eqn{\alpha = P}, \eqn{\beta = P/C}.
#'
#' The Ricker stock recruit relationship is of the following form:
#' \deqn{\textrm{Smolt} = \alpha N_1 \exp(-\beta N_2)}
#' where \eqn{\alpha = P}, \eqn{\beta = P/(Ce)}, \eqn{e} is Euler's number.
#'
#' Productivity \eqn{P} is in terms of abundance per unit of \eqn{N_1} and \eqn{N_2}.
#'
#' @seealso [calc_SRRpars()]
#' @export
calc_SRR <- function(N1, N2 = N1, p, capacity, type = c("BH", "Ricker")) {
  type <- match.arg(type)

  pars <- calc_SRRpars(p, capacity)
  alpha <- pars[1]
  beta <- pars[2]
  if (type == "BH") {
    out <- alpha * N1 / (1 + beta * N2)
  } else {
    out <- alpha * N1 * exp(-beta * N2)
  }
  return(out)
}


#' Convert density-dependent survival parameters
#'
#' Converts from capacity/productivity parameters to alpha/beta stock-recruit parameters where productivity is in terms
#' of smolts per spawner and alpha is terms of smolts per egg.
#'
#' @param p Numeric, the productivity parameter that sets the maximum survival as the initial abundance approaches zero
#' @param capacity Numeric, the capacity parameter that set the maximum survivors
#' @param f Fecundity, the spawning output per mature female
#' @param p_female The proportion of females per spawner
#' @param type Character, the functional form of the stock-recruit relationship
#' @returns A vector for alpha and beta value, respectively
#' @details
#' \deqn{\alpha = \dfrac{P}{f \times p_{female}}}
#'
#' For the Beverton-Holt stock recruit relationship:
#' \deqn{\beta = \dfrac{\alpha}{C}}
#'
#' For the Ricker stock recruit relationship:
#' \deqn{\beta = \dfrac{\alpha}{Ce}}, \eqn{e} is Euler's number.
#' @seealso [calc_SRR()]
#' @export
calc_SRRpars <- function(p, capacity, f = 1, p_female = 1, type = c("BH", "Ricker")) {
  type <- match.arg(type)
  alpha <- p/f/p_female
  beta <- switch(
    type,
    "BH" = alpha/capacity,
    "Ricker" = exp(-1) * alpha/capacity
  )
  c(alpha, beta)
}

make_SRR <- function(Bio, Habitat) {

  a <- Bio@kappa/Bio@phi
  if (Bio@SRrel == "BH") {
    b <- a/Bio@capacity_smolt
    SRRpars <- data.frame(
      a = a, b = b, phi = Bio@phi, SPRcrash = 1/a/Bio@phi, SRrel = Bio@SRrel, kappa = Bio@kappa,
      capacity_smolt = Bio@capacity_smolt, kappa_improve = Habitat@kappa_improve, capacity_smolt_improve = Habitat@capacity_smolt_improve
    )
  } else {
    b <- 1/Bio@Smax
    SRRpars <- data.frame(
      a = a, b = b, phi = Bio@phi, SPRcrash = 1/a/Bio@phi, SRrel = Bio@SRrel, kappa = Bio@kappa,
      Smax = Bio@Smax, kappa_improve = Habitat@kappa_improve, capacity_smolt_improve = Habitat@capacity_smolt_improve
    )
  }
  SRRfun <- function(SB, SRRpars) {
    if (SRRpars$SRrel == "BH") {
      SRRpars$a * SB / (1 + SRRpars$b * SB)
    } else {
      SRRpars$a * SB * exp(-SRRpars$b * SB)
    }
  }
  relRfun <- function(SSBpR, SRRpars) {
    if (SRRpars$SRrel == "BH") {
      (SRRpars$a * SSBpR - 1)/SRRpars$b/SSBpR
    } else {
      log(SRRpars$a * SSBpR)/SRRpars$b/SSBpR
    }
  }
  SPRcrashfun <- function(SSBpR0, SRRpars) SRRpars$SPRcrash
  list(SRRfun = SRRfun, SRRpars = SRRpars, relRfun = relRfun, SPRcrashfun = SPRcrashfun)
}

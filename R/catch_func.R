
#' Catch function
#'
#' Internal function that does catch calculations for the preterminal and terminal marine fisheries.
#'
#' @param NO Array `[ns, nage, n_g]` of natural-origin fish
#' @param HO Array `[ns, nage, n_r]` of hatchery-origin fish
#' @param type Character. Whether to calculate removals from harvest rate `"u"`or catch `"catch"`
#' @param U Numeric or function. Harvest rate of fishery
#' @param K Numeric or function. Total catch of the fishery
#' @param V Matrix `[ns, nage]` Relative vulnerability by age class to fishery
#' @param MSF Logical, whether fishing is mark-selective
#' @param m Numeric vector length `[ns]`, mark rate of hatchery-origin fish. Only used if `MSF = TRUE`.
#' @param release_mort Numeric vector length `[ns]`, proportion of released unmarked fish that die. Only used if `MSF = TRUE`.
#' @param p_mature_NO Optional, array `[ns, nage, n_g]` of maturity for natural-origin fish. Only use to calculate adult equivalents escapement.
#' @param p_mature_HO Optional, array `[ns, nage, n_r]` of maturity for hatchery-origin fish. Only use to calculate adult equivalents escapement.
#' @param AEQ_NO Optional, array `[ns, nage, n_g]` of adult equivalents for natural-origin fish.
#' @param AEQ_HO Optional, array `[ns, nage, n_r]` of adult equivalents for hatchery-origin fish.
#' @returns Named list:
#' - `K_NO` kept natural-origin catch, same dimension as `NO`
#' - `K_HO` kept hatchery-origin catch, same dimension as `HO`
#' - `D_NO` discarded (live + dead) natural-origin catch, same dimension as `NO`
#' - `D_HO` discarded (live + dead) hatchery-origin catch, same dimension as `HO`
#' - `DD_NO` dead discarded natural-origin catch, same dimension as `NO`
#' - `DD_HO` dead discarded hatchery-origin catch, same dimension as `HO`
#' - `Nsurv_NO` natural-origin survivors, same dimension as `NO`
#' - `Nsurv_HO` hatchery-origin survivors, same dimension as `HO`
#' - `U_NO` natural-origin harvest rate (ratio of kept catch and abundance), vector `ns`
#' - `U_HO` hatchery-origin harvest rate (ratio of dead catch and abundance), vector `ns`
#' - `Ex_NO` natural-origin exploitation rate (ratio of dead catch and abundance), vector `ns`
#' - `Ex_HO` hatchery-origin exploitation rate (ratio of dead catch and abundance), vector `ns`
#' @keywords internal
catch_func <- function(NO, HO, type = c("u", "catch"), U, K, V, MSF = FALSE, m = 1, release_mort = 0,
                       p_mature_NO = array(1, dim(NO)), p_mature_HO = array(1, dim(HO)),
                       AEQ_NO = array(1, dim(NO)), AEQ_HO = array(1, dim(HO))) {
  type <- match.arg(type)
  u_solve <- K_solve <- NA_real_

  if (type == "u") {
    if (is.numeric(U)) {
      u_solve <- U
    } else if (is.function(U)) {
      NO_obs <- apply(NO, 2, sum)
      HO_obs <- apply(HO, 2, sum)

      u_solve <- U(
        NO_obs,
        HO_obs,
        m
      )
    }
  } else {
    if (is.numeric(K)) {
      K_solve <- K
    } else if (is.function(K)) {
      NO_obs <- apply(NO, 2, sum)
      HO_obs <- apply(HO, 2, sum)

      K_solve <- K(
        NO_obs,
        HO_obs,
        m
      )
    }
  }

  # Solve for fishing effort or encounter rate that achieves the harvest rate (u_solver) or catch target (K_solve)
  ns <- dim(NO)[1]
  do_solver <- (type == "u" && u_solve > 0) || (type == "K" && K_solve > 0)
  if (do_solver && sum(NO, HO)) {

    if (MSF) {
      N_solve <- HO
      p_mature_solve <- p_mature_HO
      ret_solve <- m
      p_mature_solve <- p_mature_HO
      AEQ_solve <- AEQ_HO
    } else {
      N_solve <- abind::abind(NO, HO, along = 3)
      p_mature_solve <- abind::abind(p_mature_NO, p_mature_HO, along = 3)
      ret_solve <- rep(1, ns)
      AEQ_solve <- abind::abind(AEQ_NO, AEQ_HO, along = 3)
    }
    V_solve <- array(V, dim(N_solve))

    Emax <- 20
    opt <- try(
      uniroot(
        Effort_solver, interval = c(.Machine$double.eps, Emax),
        N = N_solve, vul = V_solve, ret = ret_solve, release_mort = release_mort,
        type = type, u = u_solve, K = K_solve,
        p_mature = p_mature_solve, AEQ = AEQ_solve
      ),
      silent = TRUE
    )
    if (is.character(opt)) {
      stop("Cannot calculate fishery effort")
      #Effort <- Emax
    } else {
      Effort <- opt$root
    }

  } else {
    Effort <- 0
  }

  ns <- dim(NO)[1]
  if (MSF) {
    ret_NO <- rep(0, ns)
    ret_HO <- m
  } else {
    ret_NO <- ret_HO <- rep(1, ns)
  }

  # K = kept catch, D = discards, DD = dead discards, U = harvest rate, Ex = exploitation rate
  K_NO <- D_NO <- DD_NO <- DK_NO <- Nsurv_NO <- array(0, dim(NO))
  K_HO <- D_HO <- DD_HO <- DK_HO <- Nsurv_HO <- array(0, dim(HO))

  # Fishing mortality
  Fret_NO <- array(V * ret_NO * Effort, dim(NO))
  Fret_HO <- array(V * ret_HO * Effort, dim(HO))

  Frel_NO <- array(V * (1 - ret_NO) * Effort, dim(NO))
  Frel_HO <- array(V * (1 - ret_HO) * Effort, dim(HO))

  # Kept catch
  K_NO[] <- (1 - exp(-Fret_NO)) * NO
  K_HO[] <- (1 - exp(-Fret_HO)) * HO

  # Discards
  D_NO[] <- (1 - exp(-Frel_NO)) * NO
  D_HO[] <- (1 - exp(-Frel_HO)) * HO

  # Dead discards
  DD_NO[] <- D_NO * release_mort
  DD_HO[] <- D_HO * release_mort

  # All dead catch
  DK_NO[] <- K_NO + DD_NO
  DK_HO[] <- K_HO + DD_HO

  # Survivors
  Nsurv_NO <- NO - K_NO - DD_NO
  Nsurv_HO <- HO - K_HO - DD_HO

  # Harvest rate and exploitation rates
  U_NO <- apply(K_NO * AEQ_NO, 1, sum)/apply(DK_NO * AEQ_NO + Nsurv_NO * p_mature_NO, 1, sum)
  U_HO <- apply(K_HO * AEQ_HO, 1, sum)/apply(DK_HO * AEQ_HO + Nsurv_HO * p_mature_HO, 1, sum)

  Ex_NO <- apply(DK_NO * AEQ_NO, 1, sum)/apply(DK_NO * AEQ_NO + Nsurv_NO * p_mature_NO, 1, sum)
  Ex_HO <- apply(DK_HO * AEQ_HO, 1, sum)/apply(DK_HO * AEQ_HO + Nsurv_HO * p_mature_HO, 1, sum)

  # Divide by zero
  U_NO[is.na(U_NO)] <- U_HO[is.na(U_HO)] <- Ex_NO[is.na(Ex_NO)] <- Ex_HO[is.na(Ex_HO)] <- 0

  output <- list(
    K_NO = K_NO,
    K_HO = K_HO,
    D_NO = D_NO,
    D_HO = D_HO,
    DD_NO = DD_NO,
    DD_HO = DD_HO,
    Nsurv_NO = Nsurv_NO,
    Nsurv_HO = Nsurv_HO,
    U_NO = U_NO,
    U_HO = U_HO,
    Ex_NO = Ex_NO,
    Ex_HO = Ex_HO
  )

  return(output)
}

#' Calculate F from harvest rate
#'
#' Solves for apical instantaneous fishing mortality rate (F), proportional to fishing effort, from harvest rate (total retained catch over total abundance).
#' The apical F can be greater than the realized F, if retention < 1. Used by [AHA()]
#'
#' @param u Harvest rate, between 0-1
#' @param K Catch, between 0-Inf
#' @param type Character, either `"catch"`, or `"u"`, whether to solve for catch or harvest rate, respectively
#' @param M Instantaneous natural mortality rate
#' @param N Abundance
#' @param vul Vulnerability
#' @param ret Retention rate
#' @param release_mort Release mortality as a proportion, between 0-1. Only relevant if `ret < 1`.
#' @param Fmax Maximum allowable value of F
#' @return Numeric for the apical F
#'
#' @keywords internal
get_F <- function(u = 0, K = 0, type = c("u", "catch"), M, N = 1, vul = 1, ret = 1, release_mort = 0,
                  Fmax = 20) {
  type <- unique(type)

  type <- match.arg(type)
  Fout <- 0

  solve_u <- type == "u" && u > 0
  solve_K <- type == "catch" && K > 0

  if (solve_u || solve_K) {
    .F <- try(
      uniroot(F_solver, interval = c(0, Fmax), M = M, N = N, vul = vul, ret = ret, release_mort = release_mort, u = u, K = K, type),
      silent = TRUE
    )

    if (is.character(.F)) {
      Fout <- Fmax
    } else {
      Fout <- .F$root
    }
  }

  return(Fout)
}


#' Solve for fishing effort
#'
#' Internal solver used by `catch_func()` to calculate the fishing effort needed to achieve target harvest rate or catch rate,
#' subject to partial retention due to mark-selective fishing.
#' Harvest rate is discounted by adult equivalents for preterminal fisheries.
#'
#' @param Eff Numeric, fishing effort
#' @param N Array `[ns, nage, n_r]`, total abundance (juvenile abundance for preterminal, return for terminal)
#' @param vul Array `[ns, nage, n_r]`, fishery vulnerability
#' @param ret Vector `[ns]`, retention rate
#' @param release_mort Vector `[ns]`, release mortality as a proportion, between 0-1. Only relevant if `ret < 1`.
#' @param type Character, either `"catch"`, or `"u"`, whether to solve for kept catch or harvest rate, respectively
#' @param u Numeric, harvest rate target
#' @param K Numeric, kept catch target
#' @param AEQ Array `[ns, nage, n_r]`, adult equivalents of catch
#' @param p_mature Array `[ns, nage, n_r]`, proportion mature by age class (used to calculate adult equivalent escapement)
#' @return Numeric. Returns the difference between the realized harvest rate (for the given value of `Eff`) and the target
#' @keywords internal
Effort_solver <- function(Eff, N, vul, ret, release_mort,
                          type = c("u", "catch"), u = 0, K = 0,
                          AEQ = array(1, dim(N)), p_mature = array(1, dim(N))) {
  type <- match.arg(type)

  F_ret <- vul * ret * Eff
  F_rel <- vul * (1 - ret) * release_mort * Eff
  F_total <- F_ret + F_rel
  dead_catch <- (1 - exp(-F_total)) * N
  catch_ret <- F_ret/F_total * dead_catch
  catch_ret[is.na(catch_ret)] <- 0

  if (type == "u") {
    Nsurv_mature <- N * exp(-F_total) * p_mature # Proportion that survive and then mature

    # This describes the aggregate harvest rate
    # In multi-population models, this does not mean that each population experiences the same exploitation rate
    fn <- sum(catch_ret * AEQ)/sum(Nsurv_mature + dead_catch * AEQ) - u
  } else {
    fn <- sum(catch_ret)/K - 1
  }
  return(fn)
}


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
#' @returns Named list:
#' - `K_NO` kept natural-origin catch, same dimension as `NO`
#' - `K_HO` kept hatchery-origin catch, same dimension as `HO`
#' - `D_NO` discarded (live + dead) natural-origin catch, same dimension as `NO`
#' - `D_HO` discarded (live + dead) hatchery-origin catch, same dimension as `HO`
#' - `DD_NO` dead discarded natural-origin catch, same dimension as `NO`
#' - `DD_HO` dead discarded hatchery-origin catch, same dimension as `HO`
#' - `U_NO` natural-origin harvest rate (ratio of kept catch and abundance), same dimension as `NO`
#' - `U_HO` hatchery-origin harvest rate (ratio of dead catch and abundance), same dimension as `HO`
#' - `Ex_NO` natural-origin exploitation rate (ratio of dead catch and abundance), same dimension as `NO`
#' - `Ex_HO` hatchery-origin exploitation rate (ratio of dead catch and abundance), same dimension as `HO`
#' @keywords internal
catch_func <- function(NO, HO, type = c("u", "catch"), U, K, V, MSF = FALSE, m = 1, release_mort = 0) {
  type <- match.arg(type)

  u_solve <- K_solve <- NA_real_

  NO_obs <- apply(NO, 2, sum)
  HO_obs <- apply(HO, 2, sum)

  ns <- dim(NO)[1]
  nage <- dim(NO)[2]

  if (type == "u") {
    if (is.numeric(U)) {
      u_solve <- U
    } else if (is.function(U)) {
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
      K_solve <- K(
        NO_obs,
        HO_obs,
        m
      )
    }
  }

  if (MSF) {
    N_solve <- apply(HO, 1:2, sum)
    ret_solve <- m
    ret_NO <- rep(0, ns)
    ret_HO <- m
  } else {
    N_solve <- apply(NO, 1:2, sum) + apply(HO, 1:2, sum)
    ret_solve <- ret_NO <- ret_HO <- rep(1, ns)
  }

  # Need to organize this into arrays by stock x age!!
  Effort <- get_F(
    u = u_solve, K = K_solve, type = type,
    M = array(0, dim(N_solve)), # ns x nage
    N = N_solve, vul = V,
    ret = matrix(ret_solve, ns, nage, byrow = TRUE),
    release_mort = matrix(release_mort, ns, nage, byrow = TRUE)
  )

  # K = kept catch, D = discards, DD = dead discards, U = harvest rate, Ex = exploitation rate
  K_NO <- D_NO <- DD_NO <- U_NO <- Ex_NO <- array(0, dim(NO))
  K_HO <- D_HO <- DD_HO <- U_HO <- Ex_HO <- array(0, dim(HO))

  # Fishing mortality
  Fret_NO <- array(V * ret_NO * Effort, dim(NO))
  Fret_HO <- array(V * ret_HO * Effort, dim(HO))

  Frel_NO <- array(V * (1 - ret_NO) * Effort, dim(NO))
  Frel_HO <- array(V * (1 - ret_HO) * Effort, dim(HO))

  K_NO[] <- (1 - exp(-Fret_NO)) * NO
  K_HO[] <- (1 - exp(-Fret_HO)) * HO

  D_NO[] <- (1 - exp(-Frel_NO)) * NO
  D_HO[] <- (1 - exp(-Frel_HO)) * HO

  DD_NO[] <- D_NO * release_mort
  DD_HO[] <- D_HO * release_mort

  U_NO[] <- K_NO/NO
  U_HO[] <- K_HO/HO

  Ex_NO[] <- (K_NO + DD_NO)/NO
  Ex_HO[] <- (K_HO + DD_HO)/HO

  # Divide by zero
  U_NO[is.na(U_NO)] <- U_HO[is.na(U_HO)] <- Ex_NO[is.na(Ex_NO)] <- Ex_HO[is.na(Ex_HO)] <- 0

  output <- list(
    K_NO = K_NO,
    K_HO = K_HO,
    D_NO = D_NO,
    D_HO = D_HO,
    DD_NO = DD_NO,
    DD_HO = DD_HO,
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
#' The apical F can be greater than the realized F, if retention < 1.
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
get_F <- function(u = 0, K = 0, type = c("u", "catch"), M, N = 1, vul = 1, ret = 1, release_mort = 0, Fmax = 20) {
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

F_solver <- function(.F, M, N = 1, vul = 1, ret = 1, release_mort = 0, u = 0, K = 0, type = c("u", "catch")) {
  type <- match.arg(type)

  F_ret <- vul * ret * .F
  F_rel <- vul * (1 - ret) * release_mort * .F
  Z <- F_ret + F_rel + M
  catch_ret <- F_ret/Z * (1 - exp(-Z)) * N
  catch_ret[is.na(catch_ret)] <- 0

  if (type == "u") {
    fn <- (1 - exp(-max(F_ret))) - u
  } else {
    fn <- sum(catch_ret)/K - 1
  }
  return(fn)

}


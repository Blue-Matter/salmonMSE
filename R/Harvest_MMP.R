
Harvest_MMP <- function(x = 1, DataList, reps = 1, u_terminal, u_preterminal,
                        p_terminal = c(2, 5), p_preterminal = c(1, 4), ...) {
  np <- length(DataList)
  nf <- length(DataList[[1]])

  multiRec <- lapply(1:np, function(p) {

    if (p %in% p_terminal) {
      Effort <- -log(1 - u_terminal)
    } else if (p %in% p_preterminal) {
      if (u_preterminal > 0) {
        y <- max(DataList[[1]][[1]]@Year) - DataList[[1]][[1]]@LHYear + 1
        nyears <- length(DataList[[1]][[1]]@Misc$FleetPars$Find[x, ])
        M <- DataList[[p]][[1]]@Misc$StockPars$M_ageArray[x, , nyears + y]
        F_preterminal <- get_F(u = u_preterminal, M = max(M))
      } else {
        F_preterminal <- 0
      }
      Effort <- F_preterminal
    } else {
      Effort <- 0
    }

    lapply(1:nf, function(f) {
      Rec <- new("Rec")
      HistE <- DataList[[p]][[f]]@OM$FinF[x] # Last historical fishing effort
      Rec@Effort <- rep(Effort/HistE, reps)
      return(Rec)
    })
  })
  return(multiRec)
}

#' Harvest component of salmon operating model
#'
#' A function that creates a multi-stock management procedure
#'
#' @param u_terminal Numeric between 0-1. Harvest rate of the terminal fishery.
#' @param u_preterminal Numeric between 0-1. Harvest rate of the preterminal fishery.
#'
#' @export
make_Harvest_MMP <- function(u_terminal = 0.1, u_preterminal = 0) {
  f <- Harvest_MMP
  formals(f)$u_terminal <- u_terminal
  formals(f)$u_preterminal <- u_preterminal
  class(f) <- "MMP"
  return(f)
}

#' @importFrom stats uniroot
get_F <- function(u = 0, M) {
  if (u > 0) {
    .F <- uniroot(F_solver, interval = c(1e-8, 3), M = M, u = u)
    .F$root
  } else {
    0
  }
}

F_solver <- function(.F, M, u = 0) {
  Z <- .F + M
  .F/Z * (1 - exp(-Z)) - u
}


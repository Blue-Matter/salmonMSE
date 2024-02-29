
Harvest_MMP <- function(x = 1, DataList, reps = 1, u, ...) {
  np <- length(DataList)
  nf <- length(DataList[[1]])

  u_total <- sum(u)
  E_total <- -log(1 - u_total)
  multiRec <- lapply(1:np, function(p) {
    lapply(1:nf, function(f) {
      Rec <- new("Rec")
      HistE <- DataList[[p]][[f]]@OM$FinF[x] # Last historical fishing effort
      Rec@Effort <- rep(E_total/HistE, reps)
      return(Rec)
    })
  })
  return(multiRec)
}

#' Harvest component of operating model
#'
#' A function that creates a multi-stock management procedure
#'
#' @param u Numeric. Harvest rate.
#'
#' @export
make_Harvest_MMP <- function(u) {
  f <- Harvest_MMP
  formals(f)$u <- u
  class(f) <- "MMP"
  return(f)
}

#' Flux computation
#'
#' Computes the flux of a metabolic pathway
#' 
#' 
#' @details 
#' Computes the flux of a metabolic pathway according to the Metabolic Control Theory (Kacser and Burns, 1973).
#' 
#' 
#' @usage flux(E_fun,A_fun,X_fun=1)
#' 
#' @param E_fun Numeric vector of concentrations   
#' @param A_fun Numeric vector of activities
#' @param X_fun Numeric value. Default is \code{1}
#' 
#'
#' @return \code{flux} returns a numeric value
#'
#' 
#' @seealso 
#' Use function \code{\link{activities}} to compute enzyme activities.
#' 
#' @references 
#' Kacser, H. and J. A. Burns, 1973. The control of flux. Symp. Soc. Exp. Biol. 27:65–104.
#' 
#' Kacser, H., J. A. Burns, H. Kacser, and D. A. Fell, 1995. The control of flux : 21 years on. Biochemical Society Transactions 23:341–366.
#'
#' @examples
#' E <- c(1,10,30)
#' A <- c(30,30,30)
#' J <- flux(E,A)
#' 
#' #result : J = 26.47059
#' 
#'
#' @export
#' 
#' 


# Flux
flux <- function(E_fun,A_fun,X_fun=1){
  if (length(A_fun)!=length(E_fun)) {
    stop("You forgot some enzymes for flux. A_fun and E_fun need to have the same length.")
  }
  J_fun <- X_fun/sum(1/(A_fun*E_fun))
  return(J_fun)
}

#' Expression of null response coefficient
#'
#' \emph{(Utilitary function)}. Gives result of a transformed expression of the response coefficient
#' 
#'
#' @details 
#' Null this expression corresponds to search the effective equilibrium in case of regulation, with or without competition.
#' 
#' 
#' @usage sol_eqeff(tau_fun,E_ini_fun,A_fun,B_fun,correl_fun)
#'     
#' @param tau_fun Numeric. Driving variable of the system
#' @param A_fun Numeric vector of activities
#' @inheritParams range_tau
#' @param correl_fun Character string indicating the constraint abbreviation. Three constraints are allowed here: \code{"RegNeg"}, \code{"CRPos"} and \code{"CRNeg"}.
#' 
#'
#'
#' @return A numeric value
#' 
#' @references Coton et al. (2021)
#' 
#' 
#' @seealso 
#' See function \code{\link{predict_eff}} to compute the effective equilibrium.
#'
#' 
#'
#' 





# Solution to null response coefficient to search the effective equilibrium in case of regulation
sol_eqeff <- function(tau_fun,E_ini_fun,A_fun,B_fun,correl_fun) {
  e0_fun <- E_ini_fun/sum(E_ini_fun)
  if (correl_fun=="RegNeg") {
    deriv_fun <- sum(1/(B_fun*A_fun*(tau_fun*(1/B_fun-e0_fun)+e0_fun)^2))
  }
  if (correl_fun=="CRPos"|correl_fun=="CRNeg") {
    deriv_fun <- sum((1/B_fun-e0_fun)/(A_fun*(tau_fun*(1/B_fun-e0_fun)+e0_fun)^2))
  }
  return(deriv_fun)
}

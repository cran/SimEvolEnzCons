#' Computation of mutation actual effect 
#'
#' Computes the actual effect \eqn{\delta} of a mutation
#' 
#'
#' @details 
#' Computes the actual effect \eqn{\delta} of a mutation depending on its canonical effect \eqn{\nu} and on the constraints.
#' 
#' 
#' 
#' @usage compute.delta(nu_fun,E_fun,correl_fun,B_fun=NULL)
#' 
#' 
#' @param nu_fun Numeric. Canonical effect of the mutation    
#' @param E_fun Numeric vector of concentrations
#' @inheritParams is.correl.authorized
#' @param B_fun Numeric vector of global co-regulation coefficients
#' 
#'
#'
#' @return Numeric vector. Each element \code{i} of the vector is the actual effect \eqn{\delta_i} for which enzyme \code{i} in \code{E_fun} is targeted by the mutation.
#' 
#' 
#' @seealso 
#' Use function \code{\link{is.correl.authorized}} to see allowed constraints for \code{correl_fun}.
#'
#' @examples
#' mu <- 1
#' E <- c(30,30,30)
#' beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
#' B <- apply(beta,1,sumbis)
#' correl <- "RegPos"
#' 
#' compute.delta(mu,E,correl,B)
#' 
#'
#' @export



compute.delta <- function(nu_fun,E_fun,correl_fun,B_fun=NULL) {
  #total concentration
  Etot_fun <- sum(E_fun)
  #number of enzymes
  n_fun <- length(E_fun)
  
  #verif param
  is.correl.authorized(correl_fun)
  is.B.accurate(B_fun,n_fun,correl_fun)
  
  
  #computes delta for each enzyme as target
  delta_fun <- rep(0,n_fun)
  if (correl_fun=="SC"|correl_fun=="RegPos"|correl_fun=="RegNeg") {
    delta_fun <- rep(nu_fun,n_fun)
  }
  if (correl_fun=="Comp") {
    delta_fun <- ((Etot_fun - E_fun)*nu_fun)/(Etot_fun + nu_fun)
  }
  if (correl_fun=="CRPos"|correl_fun=="CRNeg") {
    delta_fun <- ((Etot_fun - B_fun*E_fun)*nu_fun)/(Etot_fun + B_fun*nu_fun)
  }
  
  
  return(delta_fun)
}


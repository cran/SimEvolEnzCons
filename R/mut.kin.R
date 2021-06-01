#' Mutation of kinetic parameters
#'
#' Computes the mutant value of enzyme kinetic parameters
#' 
#'
#' @details 
#' This function used three mutation methods:
#' \itemize{
#'    \item Additive method (\code{typ_A=1}): mutant kinetic parameter is the sum of the resident one plus size of mutation \eqn{\nu}
#'    \item Multiplicative method (\code{typ_A=2}): mutant kinetic parameter is the product of the resident one and \eqn{1+\nu}
#'    \item Random method (\code{typ_A=3}): mutant kinetic parameter is equal to the mutation size
#'    }
#'    
#'  The method 1 is the default one.
#' 
#' 
#' @usage mut.kin(kin_fun,i_fun,nu_fun,typ_A=1)
#'
#' @param kin_fun Numeric vector of enzyme kinetic parameter, which are the resident values
#' @param i_fun Numeric value indicating the number of the enzyme targeted by the mutation
#' @param nu_fun Numeric value of mutation size
#' @param typ_A Numeric for mutation method. Default is 1. \emph{See details in \code{\link{mut.kin}}.}
#'
#'
#' @return Numeric vector of mutant values of enzyme kinetic parameters
#' 
#' 
#' @seealso 
#' See function \code{\link{mut.E.direct}} to compute mutation for enzyme concentrations.
#' 
#' See function \code{\link{activities}} to compute "activities" from enzyme kinetic parameters.
#'
#'
#'
#'
#' @examples
#' kin <- c(1,10,1000)
#' mu <- 1 #size of mutation
#' i <- 3 #enzyme directly targeted by mutation
#' 
#' mut.kin(kin,i,mu)
#'
#' @export




#Mutation for kinetic parameters
mut.kin <- function(kin_fun,i_fun,nu_fun,typ_A=1) {
  
  ######
  n_fun <- length(kin_fun)
  kinm <- kin_fun
  
  if (i_fun<1|i_fun>n_fun) {
    stop("You do not choose an enzyme in the sytem. i_fun need to be between 1 and n.")
  }
  
  ###### Mutation model
  # additive
  if (typ_A == 1) {
    kinm[i_fun]<-(kin_fun[i_fun]+nu_fun)
    }
  # multiplicative
  if (typ_A == 2) {
    kinm[i_fun]<-(kin_fun[i_fun]*(1+nu_fun))
    }
  # random
  if (typ_A == 3){
    kinm[i_fun] <- nu_fun
  }
  
  #flag to block negative kinetic parameters
  # if (kinm[k] <0){
  #   kinm[k] <- 0
  # }
  
  #NB : no modif for other kinetic parameters
  
  
  return(kinm)
}

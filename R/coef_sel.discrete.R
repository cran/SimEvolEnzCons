#' Selection coefficient computation
#'
#' Computes the selection coefficient using the discrete expression \eqn{s = (w^m - w^r)/w^r}
#' 
#'
#' @details 
#' Computes the selection coefficient between mutant and resident using the discrete expression \eqn{s = (w^m - w^r)/w^r},
#' assuming that fitness is proportional to flux.
#' Here, function \code{\link{flux}} is used to compute the flux, and therefore, the fitness.
#' 
#' All input vectors need to have the same length.
#' 
#' If there is no mutation affecting concentrations (resp. activities), 
#' mutant concentrations (resp. activities) are identical to resident one.
#' In this case, give to \code{E_mut} (resp. \code{A_fun}) the same value as \code{E_res} (resp. \code{A_res}), 
#' or put the default value \code{NULL}.
#' 
#' 
#' 
#' @usage coef_sel.discrete(E_res, A_res, E_mut=NULL, A_mut=NULL)
#'     
#'
#' @param E_res Numeric vector of concentrations for the resident
#' @param E_mut Numeric vector of concentrations for the mutant
#' @param A_res Numeric vector of activities for the resident.
#' Default value \code{NULL} corresponds to no mutation, \emph{i.e.} same value as \code{E_res}
#' @param A_mut Numeric vector of activities for the mutant.
#' Default value \code{NULL} corresponds to no mutation, \emph{i.e.} same value as \code{A_res}
#'  
#'
#'
#' @return Numeric value for selection coefficient
#' 
#' @seealso 
#' Use function \code{\link{activities}} to compute enzyme activities.
#' 
#'
#'
#'
#' @examples
#' 
#' ### Mutation of E
#' A <- c(1,10,30)
#' E <- c(30,30,30)
#' Em <- mut.E.direct(E,1,1,"SC")
#' 
#' coef_sel.discrete(E,A,Em)
#' 
#' 
#' ### Mutation of A
#' E <- c(30,30,30)
#' kin <- c(1,10,1000)
#' Keq <- c(10,1,100)
#' A <- activities(kin,Keq)
#' kin_m <- mut.kin(kin,3,1)
#' Am <- activities(kin_m,Keq) #equilibrium constant cannot be modified by mutation
#' 
#' coef_sel.discrete(E,A,E,Am)
#' 
#'
#' @export


# Discrete selection coefficient
coef_sel.discrete <- function(E_res, A_res, E_mut=NULL, A_mut=NULL) {
  
  ######## Set parameters
  #No mutation for concentrations
  if (length(E_mut) == 0) {
    E_mut <- E_res
  }
  
  #No mutation for activities
  if (length(A_mut) == 0) {
    A_mut <- A_res
  }
  
  
  ######### Test
  if (length(A_res)!=length(E_res)) {
    stop("An enzyme is missing before mutation. E_res and A_res need to have the same length.")
  }
  if (length(A_mut)!=length(E_mut)) {
    stop("An enzyme is missing after mutation. E_mut and A_mut need to have the same length.")
  }
  
  
  ########
  #Fitness = flux
  J_res <- flux(E_res,A_res)
  J_mut <- flux(E_mut,A_mut)
  
  #Selection coefficient
  s_fun <- J_mut/J_res - 1
  
  return(s_fun)
}


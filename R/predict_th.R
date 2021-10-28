#' Prediction of theoretical equilibrium
#'
#' Gives the theoretical equilibrium for relative concentrations
#' 
#'
#' @details 
#' Gives values at theoretical equilibrium for relative concentrations and response coefficients.
#' This equilibrium corresponds to null derivative for relative concentrations, without conditions on flux.
#' 
#' When there are regulation groups, preferably use \code{\link{predict_grp}}.
#' 
#' 
#' @usage predict_th(A_fun,correl_fun,B_fun=NULL)
#'     
#'
#' @param A_fun Numeric vector of activities
#' @inheritParams is.correl.authorized
#' @param B_fun Numeric vector of global co-regulation coefficients
#' 
#'
#'
#' @return List of two elements:
#' \itemize{
#' \item \code{$pred_e}: numeric vector of relative concentrations at theoretical equilibrium. Same length as \code{A_fun}
#' \item \code{$pred_r}: numeric vector of response coefficients at theoretical equilibrium. Same length as \code{A_fun}
#' }
#' 
#'
#' 
#' @section Special results:
#' In case of negative regulation (\code{correl_fun} = \code{"RegNeg"} or \code{"CRNeg"}), relative concentrations would be negative.
#' 
#' In case of competition plus regulation (\code{correl_fun} = \code{"CRPos"} or \code{"CRNeg"}), response coefficients is not defined and \code{$pred_r} returns \code{NaN}.
#' 
#' 
#' @seealso 
#' Use function \code{\link{activities}} to compute enzyme activities.
#' 
#' Use function \code{\link{is.correl.authorized}} to see allowed constraints for \code{correl_fun}.
#' 
#' Use function \code{\link{predict_grp}} to predict equilibria when there are regulation groups.
#'
#' @examples
#' #### For independancy "SC" or competition "Comp"
#' A <- c(1,10,30)
#' 
#' eq_th <- predict_th(A,"SC")
#' 
#' eq_th$pred_e
#' eq_th$pred_r
#' 
#' ###### In presence of regulation
#' A <- c(1,10,30)
#' beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
#' B <- apply(beta,1,sumbis) 
#' 
#' eq_th <- predict_th(A,"CRPos",B)
#' 
#' eq_th$pred_e
#' eq_th$pred_r
#' 
#'
#' @export



predict_th <- function(A_fun,correl_fun,B_fun=NULL){
  #Number of enzymes
  n_fun <- length(A_fun)
  
  ###### Verif parameters
  # cases allowed; if not, program stops
  is.correl.authorized(correl_fun)
  #test B_fun
  is.B.accurate(B_fun,n_fun,correl_fun)

  # #value of B in these cases is not important, and this writing avoid multiple test
  # if (correl_fun=="SC"|correl_fun=="Comp") {
  #   B_fun <- rep(1,n_fun)
  # }
  # #in other cases (presence of regulation), verif if a parameter is missing
  # if (length(B_fun)==0) {
  #   stop("B_fun is required.")
  # }
  # if (length(B_fun)!=n_fun) {
  #   stop("B_fun and A_fun need to have the same length. An enzyme is missing.")
  # }  

  
  
  #########
  
  #Independency, modele of mutation number 1
  if (correl_fun=="SC"){
    estar <- A_fun^(-1/3)/sum(A_fun^(-1/3))
    rstar <- A_fun^(-2/3)/sum(A_fun^(-2/3))
  }
  
  #Competition for ressources
  if (correl_fun=="Comp") {
    estar <- (A_fun^(-1/2)/sum(A_fun^(-1/2)))
    rstar <- rep(0,n_fun)
  }
  
  #Regulation between all enzymes
  if (correl_fun =="RegPos"|correl_fun=="RegNeg") {
    estar <- 1/B_fun
    rstar <- rep(1,n_fun)
  }
  
  #Competition and regulation between all enzymes
  if(correl_fun=="CRPos"|correl_fun=="CRNeg"){
    estar <- 1/B_fun
    #Response coefficients is not defined at this point
    rstar <- rep(NaN,n_fun)
  }
  
  return(list(pred_e=estar,pred_r=rstar))
}


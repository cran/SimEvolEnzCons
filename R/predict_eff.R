#' Prediction of effective equilibrium
#'
#' Gives the effective equilibrium for relative concentrations
#' 
#'
#' @details 
#' Gives values at effective equilibrium for relative concentrations and corresponding driving variable \eqn{\tau}.
#' This equilibrium corresponds to null derivative of relative concentrations, with a maximum for flux.
#' 
#' Effective equilibrium is found by searching the zero for response coefficients.
#' The R function uses in this objective is \code{\link[stats]{uniroot}}.
#' 
#' Note that sum of \code{1/B_fun} need to be equal to 1.
#' 
#' 
#' @importFrom stats uniroot
#' 
#' 
#' @usage predict_eff(E_ini_fun,B_fun,A_fun,correl_fun, tol=0.00000001)
#'     
#'
#' @param A_fun Numeric vector of activities
#' @inheritParams range_tau
#' @inheritParams is.correl.authorized
#' @param tol Tolerance for function \code{\link[stats]{uniroot}}
#' 
#'
#'
#' @return List of three elements:
#' \itemize{
#' \item \code{$pred_e}: numeric vector of relative concentrations at effective equilibrium. Same length as \code{A_fun}
#' \item \code{$pred_tau}: numeric value of driving variable \emph{tau} at effective equilibrium
#' \item \code{$pred_E}: numeric vector of absolute concentrations at effective equilibrium. Same length as \code{A_fun}
#' }
#' 
#'
#' 
#' @section Special results:
#' In case of independence (\code{correl_fun="SC"}) or positive regulation (\code{correl_fun="RegPos"}), there is no effective equilibrium, and function \code{predict_eff} stops.
#' 
#' In case of competition (\code{correl_fun="Comp"}), effective and theoretical equilibria are confounded. Function \code{predict_eff} also stops, so use preferably function \code{\link{predict_th}} to compute equilibrium.
#' 
#' If \code{E_ini_fun} is a multiple of \code{1/B_fun}, effective equilibrium is confounded with theoretical equilibrium and initial point (see \code{\link{droites}} for details).
#' Function \code{predict_eff} returns \code{E_ini_fun} for \code{$pred_E} and 0 for \code{$pred_tau}, with a warning message.
#' 
#' 
#' 
#' @seealso 
#' Use function \code{\link{activities}} to compute enzyme activities.
#' 
#' Use function \code{\link{is.correl.authorized}} to see allowed constraints for \code{correl_fun}.
#' 
#' Use function \code{\link{predict_th}} to compute theoretical equilibrium.
#' 
#' 
#' @references 
#' Coton et al. (2021)
#'
#'
#'
#'
#' @examples
#' ###### In presence of competition plus regulation
#' A <- c(1,10,30)
#' E0 <- c(30,30,30)
#' beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
#' B <- apply(beta,1,sumbis) 
#' 
#' eq_eff <- predict_eff(E0,B,A,"CRPos")
#' 
#' eq_eff$pred_e
#' eq_eff$pred_tau
#' eq_eff$pred_E
#' 
#'
#' @export




predict_eff <- function(E_ini_fun,B_fun,A_fun,correl_fun, tol=0.00000001) {
  ###" Config
  #Number of enzymes
  n_fun <- length(A_fun)
  #To avoid problem due to use of dataframe
  E_ini_fun <- as.matrix(E_ini_fun)
  # Computes initial relative concentrations
  e0_fun <- E_ini_fun/sum(E_ini_fun)
  
  
  ######### Verif parameters
  #authorized parameters
  is.correl.authorized(correl_fun)
  # missing parameters
  if (length(E_ini_fun)!=n_fun) {
    stop("An enzyme are missing. A_fun and E_ini_fun have not the same length.")
  }
  #adequaty of B_fun
  is.B.accurate(B_fun,n_fun,correl_fun)
  
  
  
  ########## Special cases
  # No solution
  if (correl_fun=="SC"|correl_fun=="RegPos") {
    pred_e <- rep(NA,n_fun)
    stop(paste("Effective equilibrium does not exists in the case '",correl_fun,"'.",sep=""))
  }
  #Same solution as theoretical equilibrium
  if (correl_fun=="Comp") {
    pred_e <- predict_th(A_fun,correl_fun)$pred_e
    stop("In case of competition, effective and theoretical equilibria are identical. Preferably use function predict_th(). ")
  }
  #If theoretical equilibrium confounded with initial relative concentration, nothing moves
  if (all(e0_fun==(1/B_fun))) {
    pred_e <- predict_th(A_fun,correl_fun,B_fun)$pred_e
  }


  
  
  #######
  #Computes bounds of tau
  #Such as all relative concentrations are between 0 and 1, because dJ/dtau has several extrema
  range_fun <- range_tau(E_ini_fun,B_fun)

  
  # ########## Part coding before range_tau(): don't need anymore
  # # For the chosen interval of tau, test if relative concentrations are between 0 and 1
  # 
  # # Computes interval of tested tau for searching the zero of the function
  # to<- seq(min(range_fun),max(range_fun), by=0.01)
  # #Matrix of relative concentrations on line D (nb rows =length of tested tau, nb col = nb enzymes)
  # e_test_fun <- matrix(rep(0,length(to)*n_fun),nrow=length(to),ncol=n_fun)
  # e_fun <- NULL
  # tu <- rep(0,length(to))
  # tau_fun <- NULL
  # #for all tested values of tau
  # for (i in 1:length(to)){
  #   #Testing if all associated relative concentrations associated to to[i] are between 0 and 1
  #   if ((sum(sapply(to[i],droite_e,E_ini_fun,B_fun)>=0)==n_fun)&(sum(sapply(to[i],droite_e,E_ini_fun,B_fun)<=1)==n_fun)) {
  #     #only if is true, values of 
  #     e_test_fun[i,] <- t(sapply(to[i],droite_e,E_ini_fun,B_fun)) #on ne garde que ces valeurs
  #     tu[i] <- to[i]
  #   }
  # }
  # #conserving only values such as sum of relative concentrations are different from 0, i.e. the values not retained in matrix e_test_fun
  # e_fun <- e_test_fun[apply(e_test_fun,1,sum)!=0,]
  # tau_fun <- tu[apply(e_test_fun,1,sum)!=0]
  # #I could use rbind() instead of e_test_fun and tu...
  
  
  
  ########"
  
  if (any(e0_fun!=1/B_fun)) {
    #Search of zero for response coefficient
    zero_ini <-uniroot(f=sol_eqeff, interval=c(min(range_fun)+tol,max(range_fun)-tol),E_ini_fun,A_fun,B_fun,correl_fun, tol=tol)
    #zero_ini <-uniroot(f=sol_eqeff, interval=c(min(tau_fun+tol),max(tau_fun-tol)),E_ini_fun,A_fun,B_fun,correl_fun, tol=tol)
    #Recuperation of tau
    tau_Jmax_fun <- zero_ini$root
    #recuperation of associated relative concentrations
    e_Jmax_fun <- droite_e(zero_ini$root,E_ini_fun,B_fun)
  } else {
    #if all e0=1/B, line becomes a point, so we put initial value
    tau_Jmax_fun <- 0
    e_Jmax_fun <- e0_fun
  }
  
  #recuperation of associated absolute concentrations
  if (correl_fun=="RegNeg") {
    E_Jm_fun <- droite_E.Reg(tau_Jmax_fun,E_ini_fun,B_fun)
  } else {
    #in other cases, ie with competition, Etot is fixed
    E_Jm_fun <- e_Jmax_fun*sum(E_ini_fun)
  }
  
  

  return(list(pred_e=e_Jmax_fun,pred_tau=tau_Jmax_fun,pred_E=E_Jm_fun))
}


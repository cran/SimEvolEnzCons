#' Prediction of effective equilibrium for all possible initial relative concentrations
#'
#' Gives the effective equilibrium for relative concentrations for various initial concentrations
#' 
#'
#' @details 
#' Effective equilibrium is computed with function \code{\link{predict_eff}}.
#' 
#' \emph{\bold{WARNING: Function \code{predict_eff_allE0} is only available for three enzymes!} Length of \code{A_fun} and \code{B_fun} need to be 3.}
#' 
#' Each relative concentration is taken between 0 and 1 by 0.01, then triplet of relative concentrations are sorted to have a sum equal to 1.
#' Then relative concentrations are multiplied by \code{Etot_fun} to have initial concentrations.
#' 
#' For parameter \code{correl_fun}, authorized input are \code{"RegNeg"}, \code{"CRPos"} and \code{"CRNeg"}.
#' 
#' @usage predict_eff_allE0(B_fun,A_fun,correl_fun,Etot_fun=100,X_fun=1, tol=0.00000001)
#'     
#'
#' @inheritParams predict_eff
#' @inheritParams flux.dome.graphics
#' 
#'
#'
#' @return Invisible list of 3 elements:
#' \itemize{
#'    \item \code{$all_eq_eff}: Dataframe of 4 851 rows and eight columns (named \code{e1,e2,e3,tau,E1,E2,E3,J}) for effective equilibrium from possible initial concentrations.
#'    Each row corresponds to a set of initial concentrations, and columns are respectively relative concentrations (\code{$e1,$e2,$e3}), driving variable \eqn{\tau} (\code{$tau}), absolute concentrations \code{$E1,$E2,$E3} and flux \code{$J} at effective equilibrium; 
#'    \item \code{$all_E0}: Dataframe of 4 851 rows and three columns corresponding to initial concentrations. Each row is a triplet of initial concentrations;
#'    \item \code{$param}: List of input parameters
#' }
#'
#'  
#' 
#' 
#'
#' @export




predict_eff_allE0 <- function(B_fun,A_fun,correl_fun,Etot_fun=100,X_fun=1, tol=0.00000001) {
  ###" Config
  #Number of enzymes
  n_fun <- length(A_fun)
  

  
  
  ######### Verif parameters
  #authorized parameters
  is.correl.authorized(correl_fun)
  #adequaty of B_fun
  is.B.accurate(B_fun,n_fun,correl_fun)
  
  
  #######" Initialization
    #Computes all possible e0 (ei0 no-null and sum(e0)=1)
  e1<-rep(seq(0,1, by=0.01), 100) #nb pts in x
  e2<-rep(seq(0,1, by=0.01), each=100) #nb pts in y
  e3<-rep(seq(0,1, by=0.01), each=10100) #nb pts in z
  e_sup<-cbind(e1,e2,e3)
  #only triplet such as sum(e0)=1
  e_tri <- Etot_fun*e_sup[apply(e_sup,1,sum)==1&(e_sup[,1]!=0&e_sup[,2]!=0&e_sup[,3]!=0),]

  
  #######
  #vector of all tau
  tau_Jmax <- rep(0,nrow(e_tri)) #NULL 
  e_Jmax <- matrix(0,ncol=3,nrow=nrow(e_tri))
  E_Jmax <- matrix(0,ncol=3,nrow=nrow(e_tri))
  Jmax <- rep(0,nrow(e_tri))
  E0_all <- matrix(0,ncol=3,nrow=nrow(e_tri))
  
  # for all possible e0
  for (k in 1:nrow(e_tri)){
    #define initial value used
    E0_fun <- e_tri[k,]
    #compute effective equilibrium
    pred <- predict_eff(E0_fun,B_fun,A_fun,correl_fun,tol=tol)
    J_fun <- flux(pred$pred_E,A_fun,X_fun)
    #recuperation of results
    tau_Jmax[k] <- pred$pred_tau #rbind(tau_Jmax,pred$pred_tau)
    e_Jmax[k,] <- pred$pred_e #rbind(e_Jmax,t(pred$pred_e))
    E_Jmax[k,] <- pred$pred_E #rbind(E_Jmax,t(pred$pred_E))
    Jmax[k] <- J_fun #rbind(Jmax,J_fun)
    E0_all[k,] <- E0_fun #rbind(E0_all,E0_fun)
  }
  #Assembling data
  tau_Jmax <- as.data.frame(tau_Jmax)
  colnames(tau_Jmax) <- c("tau")
  e_Jmax <- as.data.frame(e_Jmax)
  colnames(e_Jmax) <- c("e1","e2","e3")
  E_Jmax <- as.data.frame(E_Jmax)
  colnames(E_Jmax) <- c("E1","E2","E3")
  Jmax <- as.data.frame(Jmax)
  colnames(Jmax) <- c("J")
  ens_Jmax <- cbind.data.frame(e_Jmax,tau_Jmax,E_Jmax,Jmax)
  E0_all <- as.data.frame(E0_all)
 

  return(invisible(list(all_eq_eff=ens_Jmax,all_E0=E0_all,param=list(A=A_fun,B=B_fun,correl=correl_fun,X=X_fun,Etot=Etot_fun))))
}


#' Simulation of enzyme evolution
#'
#' This function simulates evolution of enzyme concentrations under constraints
#' 
#'
#' @details 
#' Time step is defined between appearance of two mutation.
#' 
#' Time step is simulated with function \code{\link{simul.next.resident}}, which gives values of resident after a time step.
#' 
#' To reduce size of result matrix, certain results only are conserved.
#' State of simulation is observed every \code{pasobs} time steps.
#' There is \code{npt} observations, so result matrix has \code{npt} rows.
#' In total, a simulation includes \code{pasobs*npt} time steps. By default, there is \code{125 000} time steps.
#' 
#' Chosen equilibrium for \code{pred_enzsim} and \code{pred_contsim} are the theoretical equilibrium for constraints \code{"SC"}, \code{"Comp"} and \code{"RegPos"}, 
#' and the effective one for constraints \code{"RegNeg"}, \code{"CRPos"} and \code{"CRNeg"}. If there is regulation groups (\code{sum(1/B)>1}), the the theoretical equilibrium \emph{within groups} is given.
#' 
#' 
#' 
#' @usage simul.evol.enz.one(E_ini_fun, kin_fun, Keq_fun, N_fun, correl_fun, 
#' beta_fun=NULL, X_fun=1, pasobs=250, npt=500, max_mut_size_E=1, max_mut_size_A=1,
#' pmutA=0, typ_E=1, typ_A=1, use.old.mut=FALSE)
#'
#' @param pasobs Numeric. Number of time steps between two successive observations of the system. Default is \code{250}
#' @param npt Numeric. Number of observations. Default is \code{500}
#' @param E_ini_fun Numeric vector of the initial concentrations
#' @param kin_fun Numeric vector of the initial kinetic parameters
#' @inheritParams simul.next.resident
#'  
#'
#'
#' @return Invisible list of 4 elements:
#' \itemize{
#'    \item \code{$res_sim}: numeric matrix of \code{npt} rows and \code{3*n+3} columns.
#'    Each row corresponds to state at each observation step (i.e. after \code{pasobs} mutations),
#'    and columns are respectively concentrations (\code{1:n}), kinetic parameters (\code{n+1:2n}), total concentration (\code{2n+1}), total kinetic (\code{2n+2}), flux/fitness (\code{2n+3}), activities (\code{2n+4:3n+3});
#'    \item \code{$pred_enzsim}: numeric matrix of \code{npt} rows and \code{n} columns, corresponding to relative concentrations at equilibrium, for each observation step (in rows);
#'    \item \code{$pred_contsim}: same as \code{$pred_enzsim}, but for response coefficients
#'    \item \code{$param}: list of input parameters: \itemize{
#'        \item \code{n}: number of enzymes,
#'        \item \code{E0}: numeric vector of initial concentrations,
#'        \item \code{kin}: numeric vector of initial kinetic parameters,
#'        \item \code{Keq}: numeric vector of constant equilibrium,
#'        \item \code{beta}: matrix of co-regulation coefficients,
#'        \item \code{B}: numeric vector of global co-regulation coefficients,
#'        \item \code{correl}: character string indicating the constraint abbreviation,
#'        \item \code{N}: population size,
#'        \item \code{pasobs}: number of steps between two system observations,
#'        \item \code{npt}: number of system observations,
#'        \item \code{X}: parameter for flux computation,
#'        \item \code{pmutA}: probability for activity mutation,
#'        \item other input parameters}
#' }
#' 
#' Note that \code{n} is the number of enzymes, which is the length of \code{E_ini_fun}.
#' 
#' 
#' @seealso 
#' See function \code{\link{simul.next.resident}} to see how works each time steps. 
#' 
#' Use function \code{\link{simul.evol.enz.multiple}} to compute several simulations.
#'
#'
#'
#'
#'
#' @export




# Simulations of enzyme evolution
simul.evol.enz.one <- function(E_ini_fun, kin_fun, Keq_fun, N_fun, correl_fun, beta_fun=NULL, 
                           X_fun=1, pasobs=250, npt=500, max_mut_size_E=1, max_mut_size_A=1,
                           pmutA=0, typ_E=1, typ_A=1, use.old.mut=FALSE) {
  
  
  ########### Settings
  n_fun <- length(E_ini_fun)
  if (length(beta_fun)==0) {
    beta_fun <- diag(rep(1,n_fun))
  }
  B_fun <- apply(beta_fun,1,sumbis)

  #Number of time steps/generations = observation steps (i.e. number of time step between two observations) x number of observations
  ngene <- pasobs*npt
  
  ########## Test
  is.correl.authorized(correl_fun)
  is.beta.accurate(beta_fun,n_fun,correl_fun)
  if (any(E_ini_fun<0)|any(kin_fun<0)|N_fun<0) {
    stop("Enter only positive numeric for E_ini_fun, A_fun and N_fun")
  }
  if (pasobs<0|npt<0) {
    stop("Simulations go only to the future. Enter positive value for pasobs and npt.")
  }
  
  ########## Initialization
  
  #Matrix for results of simulations
  #Each row of resobs = one observation, after pasobs time steps
  #E, kin, sum(E), sum(kin), J, A (npt lignes x 3*n+3 colonnes)
  resobs<-matrix(0,nrow=npt,ncol=(3*n_fun+3))
  
  #First resident
  E_res <- E_ini_fun
  kin_res <- kin_fun
  
  #Activities and flux
  A_fun <- activities(kin_fun,Keq_fun)
  J_fun <- flux(E_ini_fun,A_fun,X_fun)
  
  # First data
  resobs[1,]<-c(E_ini_fun,kin_fun,sum(E_ini_fun),sum(kin_fun),J_fun,A_fun)
  
  
  
  ## Predictions for theoretical equilibrium
  pred_enz<-matrix(0,nrow=npt,ncol=n_fun) # for enzyme concentrations e* or tilde{e}
  pred_cont<-matrix(0,nrow=npt,ncol=n_fun) # for response coefficients R*
  #First data
  star <- predict_th(A_fun,correl_fun,B_fun)
  pred_cont[1,] <- star$pred_r
  if (correl_fun=="RegNeg"|correl_fun=="CRPos"|correl_fun=="CRNeg") {
    if (sum(1/B_fun)==1) {
      star <- predict_eff(E_ini_fun,B_fun,A_fun,correl_fun)
    }
  }
  pred_enz[1,] <- star$pred_e
  
  
  
  
  ############## Simulation
  
  #For each time step
  for (l in (pasobs+1):ngene) {
    #beginning from pasobs+1, because if l=pasobs, (l/pasobs)=1 and initial values are modified
    #so there is (ngene-pasobs) mutation time steps
    
    
    #Computes next resident values
    next_gen <- simul.next.resident(E_res,kin_res,Keq_fun,N_fun,correl_fun,beta_fun,
                                    X_fun,max_mut_size_E,max_mut_size_A,pmutA,typ_E,typ_A,use.old.mut)
    
    A_next <- activities(next_gen$kin_next,Keq_fun)
    
    #Set new resident values
    E_res <- next_gen$E_next
    kin_res <- next_gen$kin_next
    
    
    ###### Registration values every 'pasobs' time steps, named observation steps
    
    # => accumulates mutations until (l/pasobs) is an integer
    #test if it is an observation step, by verifying if (l/pasobs) is an integer
    if(l/pasobs==floor(l/pasobs)) {
      #register resident values
      resobs[l/pasobs,] <- c(next_gen$E_next,next_gen$kin_next,next_gen$Etot_next,next_gen$kintot_next,next_gen$w_next,A_next)
      #register resident predictions
      star <- predict_th(A_next,correl_fun,B_fun)
      pred_cont[l/pasobs,] <- star$pred_r
      if (correl_fun=="RegNeg"|correl_fun=="CRPos"|correl_fun=="CRNeg") {
        if (sum(1/B_fun)==1) {
          star <- predict_eff(E_ini_fun,B_fun,A_next,correl_fun)
        }
      }
      pred_enz[l/pasobs,] <- star$pred_e
      
      
      # Note: probably simpler but longer with rbind()...
    }
    
  }
  
  
  return(invisible(list(res_sim=resobs,pred_enzsim=pred_enz,pred_contsim=pred_cont,
              param=list(E0=E_ini_fun,kin=kin_fun,Keq=Keq_fun,n=n_fun,beta=beta_fun,B=B_fun,
                         correl=correl_fun,N=N_fun,pasobs=pasobs,npt=npt,
                         X=X_fun,max_mut_size_E=max_mut_size_E,max_mut_size_A=max_mut_size_A,
                         pmutA=pmutA,typmut_E=typ_E,typmut_A=typ_A,nsim=1))))
}



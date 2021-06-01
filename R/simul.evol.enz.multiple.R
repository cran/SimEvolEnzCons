#' Multiple simulations of enzyme evolution
#'
#' This function gives multiple simulations of evolution of enzyme concentrations under constraints
#' 
#'
#' @details 
#' Details about how does simulations work are in function \code{\link{simul.evol.enz.one}} documentation. 
#' 
#' \emph{Apply would also work, but only for multiple values of \code{E_ini_fun}.
#' \code{simul.evol.enz.multiple} gives all results in a same table, and has the possibility to change kinetic parameters between simulation.}
#' 
#' \bold{Initial parameters}
#' 
#' You can choose if initial concentrations and initial kinetic parameters are identical between all simulations by setting \code{same.E0} and \code{same.kin0}.
#' 
#' You can choose randomly initial concentrations and initial kinetic parameters by setting \code{is.random.E0} and \code{is.random.kin0}.
#' 
#' Four cases are possible. Examples are given for initial concentrations, but are the same for kinetic parameters, where \code{same.E0}, \code{is.random.E0} and \code{E_ini_fun} are respectively \code{same.kin0}, \code{is.random.kin0} and \code{kin_fun}.
#' \itemize{
#'    \item \code{same.E0} is \code{TRUE} and \code{is.random.E0} is \code{FALSE} is default parameters. In this case, input \code{E_ini_fun} of length \code{n} is used for all simulations.
#'    \item If \code{same.E0} and \code{is.random.E0} are \code{TRUE}, \code{E_ini_fun} is chosen randomly and used for all \code{nsim} simulations.
#'    \item If \code{same.E0} is \code{FALSE} and \code{is.random.E0} \code{TRUE}, each simulation has different random initial concentrations.
#'    \item If \code{same.E0} and \code{is.random.E0} are \code{FALSE}, a matrix of \code{n} columns and \code{nsim} rows is needed for \code{E_ini_fun}, each row corresponding to one simulation.
#' }
#' 
#' \code{$list_init} is the return for initial concentrations, kinetic parameters and activities.
#' 
#' 
#' \bold{Other details}
#' 
#' Note that \code{n} is the number of enzyme and is determined by length of \code{Keq_fun}.
#'
#' 
#'   
#' 
#' @usage simul.evol.enz.multiple(E_ini_fun, kin_fun, Keq_fun, nsim, 
#' N_fun, correl_fun, beta_fun=NULL, X_fun=1, pasobs=250, npt=500, 
#' same.E0=TRUE, is.random.E0=FALSE, same.kin0=TRUE, is.random.kin0=FALSE,
#' Etot_fun=100, kin_max=1000, max_mut_size_E=1, max_mut_size_A=1,
#' pmutA=0, typ_E=1, typ_A=1, use.old.mut=FALSE)
#'
#' 
#' @inheritParams simul.evol.enz.one
#' @inheritParams flux.dome.graphics
#' @param Etot_fun Numeric. Total concentration if \code{is.random.E0} is \code{TRUE}. Default is \code{100}.
#' @param nsim Numeric. Number of simulations
#' @param kin_max Numeric. Maximal value for random \code{kin_fun}. Default is \code{1000}.
#' @param same.E0,same.kin0 Logical. Uses same \code{E_ini_fun} / \code{kin_fun} for all simulations ? Default is \code{TRUE}.
#' @param is.random.E0,is.random.kin0 Logical. Is \code{E_ini_fun} / \code{kin_fun} chosen randomly ? Default is \code{FALSE}.
#'  
#'
#'
#' @return Invisible list of 6 elements:
#' \itemize{
#'    \item \code{$tabR}: dataframe of \code{nsim x npt} rows and \code{3*n+4} columns.
#'    Each row corresponds to state at each observation step (i.e. after \code{pasobs} mutations),
#'    and columns are respectively concentrations (\code{1:n}), kinetic parameters (\code{n+1:2n}), total concentration (\code{2n+1}), total kinetic (\code{2n+2}), flux/fitness (\code{2n+3}), activities (\code{2n+4:3n+3}), and simulation number (\code{$sim} or \code{3n+4});
#'    \item \code{$tabP_e}: numeric matrix of \code{npt} rows and \code{n+1} columns, corresponding to relative concentrations at equilibrium (column 1 to \code{n}) for each observation step (in rows), plus column \code{$sim};
#'    \item \code{$tabP_r}: same as \code{$tabP_e}, but for response coefficients
#'    \item \code{$list_init}: list of 3 elements, containing initial values of concentrations in \code{$E0}, kinetic parameters in \code{$kin0} and activities in \code{$A0} for each simulation. Each element is a numeric matrix of \code{nsim} rows (one by simulation) and \code{n} columns (one by enzyme);
#'    \item \code{$list_final}: list of 3 elements, containing final values of concentrations in \code{$E_f}, kinetic parameters in \code{$kin_f} and activities in \code{$A_f} for each simulation. Each element is a numeric matrix of \code{nsim} rows and \code{n} columns;
#'    \item \code{$param}: list of input parameters.
#' }
#' 
#' Note that \code{n} is the number of enzymes, which is the length of \code{E_ini_fun}.
#' 
#' 
#' @seealso 
#' See function \code{\link{simul.evol.enz.one}} to see how works a simulation.
#' 
#'
#' @examples 
#' 
#' \donttest{
#' # case for 3 enzymes
#' n <- 3
#' E0 <- c(30,30,30)
#' kin <- c(1,10,30)
#' Keq <- c(1,1,1)
#' nsim <- 2 # 2 simulations
#' N <- 1000
#' beta <- diag(1,n)
#' beta[upper.tri(beta)] <- c(0.32,0.32*(-0.43),-0.43)
#' #beta_12 = 0.32, beta_13 = beta_12 x beta_23, beta_23 = -0.43
#' t_beta <- t(beta) #because R fills matrix column by column
#' beta[lower.tri(beta)] <- 1/t_beta[lower.tri(t_beta)] #beta_ji = 1/beta_ij
#' if (n==3) {beta[lower.tri(beta)] <- 1/beta[upper.tri(beta)]} #only available if n=3
#' correl <- "RegNeg"
#' 
#' evol_sim <- simul.evol.enz.multiple(E0,kin,Keq,nsim,N,correl,beta)
#' evol_sim$tabR[,1:n] #concentrations
#' evol_sim$tabR[,(n+1):(2*n)] #kinetic parameters
#' evol_sim$tabR[,(2*n+1)] #total concentration
#' evol_sim$tabR[,(2*n+2)] #total kinetic
#' evol_sim$tabR[,(2*n+3)] #flux
#' evol_sim$tabR[,(2*n+4):(3*n+3)] #activities
#' evol_sim$tabR$sim #simulation number
#' evol_sim$tabR[evol_sim$tabR$sim==2,] #results for 2nd simulation
#' }
#'
#' @export



# Simulations of enzyme evolution
simul.evol.enz.multiple <- function(E_ini_fun, kin_fun, Keq_fun, nsim, N_fun, correl_fun, beta_fun=NULL, 
                           X_fun=1, pasobs=250, npt=500, same.E0=TRUE, is.random.E0=FALSE, same.kin0=TRUE, is.random.kin0=FALSE,
                           Etot_fun=100, kin_max=1000, max_mut_size_E=1, max_mut_size_A=1,
                           pmutA=0, typ_E=1, typ_A=1, use.old.mut=FALSE) {
  
  
  ########### Settings
  n_fun <- length(Keq_fun)
  if (length(beta_fun)==0) {
    beta_fun <- diag(rep(1,n_fun))
  }
  B_fun <- apply(beta_fun,1,sumbis)
  
  
  ########## Test
  is.correl.authorized(correl_fun)
  is.beta.accurate(beta_fun,n_fun,correl_fun)
  if (any(E_ini_fun<0)|any(kin_fun<0)|N_fun<0) {
    stop("Enter only positive numeric for E_ini_fun, A_fun and N_fun")
  }
  if (pasobs<0|npt<0) {
    stop("Simulations go only to the future. Enter positive value for pasobs and npt.")
  }
  #Length of input E_ini_fun
  #if E0 is not random
  if (is.random.E0==FALSE) {
    #if use same E0 for all simulations
    if (same.E0==TRUE) {
      if (length(E_ini_fun)!=n_fun) {
        stop("To use same initial concentrations, please put a vector of length 'n' for 'E_ini_fun'.")
      }
      
      #if use different E0 for each simulations
    } else {
      # if E0 is a vector and there more than one simulation
      if (class(E_ini_fun)=="numeric") {
        if (nsim!=1) {
          stop("To use different initial concentrations, please put a matrix of 'nsim' rows and 'n' columns for 'E_ini_fun'.")
        }
      }
      #if E0 is a matrix but do not have n_fun columns or nsim rows
      if (class(E_ini_fun)=="matrix") {
        if (ncol(E_ini_fun)!=n_fun|nrow(E_ini_fun)!=nsim) {
          stop("To use different initial concentrations, please put a matrix of 'nsim' rows and 'n' columns for 'E_ini_fun'.")
        }
      }
      
    }
  }
  #Length of input kin_fun
  #if kin0 is not random
  if (is.random.kin0==FALSE) {
    #if use same kin0 for all simulations
    if (same.kin0==TRUE) {
      if (length(kin_fun)!=n_fun) {
        stop("To use same initial kinetic parameters, please put a vector of length 'n' for 'kin_fun'.")
      }
      
      #if use different kin0 for each simulations
    } else {
      # if kin0 is a vector and there more than one simulation
      if (class(kin_fun)=="numeric") {
        if (nsim!=1) {
          stop("To use different kinetic parameters, please put a matrix of 'nsim' rows and 'n' columns for 'kin_fun'.")
        }
      }
      #if kin0 is a matrix but do not have n_fun columns or nsim rows
      if (class(kin_fun)=="matrix") {
        if (ncol(kin_fun)!=n_fun|nrow(kin_fun)!=nsim) {
          stop("To use different kinetic parameters, please put a matrix of 'nsim' rows and 'n' columns for 'kin_fun'.")
        }
      }
    }
  }
  
  
  ########## Initialization
  
  ## Matrix for results of simulations
  #Each row = one observation, and $sim for simulation number
  #E, kin, sum(E), sum(kin), J, A, n sim (npt lignes x 3*n+4 colonnes)
  tabR <- NULL
  
  ## Predictions for theoretical equilibrium
  tabP_e <- NULL
  tabP_r <- NULL
  
  ## List of initial and final values for each simulations
  #concentrations, kinetic parameters, activities
  #Each element of the list is a matrix of n_fun columns and nsim rows
  initial_values <- list(E0=NULL,kin0=NULL,A0=NULL)
  final_values <- list(E_f=NULL,kin_f=NULL,A_f=NULL)
  
  #if use same E0 for all simulations => set E0_sim now
  if (same.E0==TRUE) {
    #set E0_sim, which be used for all simulations
    
    #if E0 is not random
    if (is.random.E0==FALSE) {
      #take input value
        E0_sim <- E_ini_fun
        
        #if E0 is random
      } else {
        #take random values
        E0_alea <- runif(n_fun,0,Etot_fun)
        #recalibrate to have summ(E0_sim) = Etot_fun
        E0_sim <- E0_alea*Etot_fun/sum(E0_alea)
      }
    
    # if different E0 for each simulations => set E0_sim in each simulation
  }
  #idem for kin0
  if (same.kin0==TRUE) {
    if (is.random.kin0==FALSE) {
      kin0_sim <- kin_fun
    } else {
      #take random values
      kin0_sim <- runif(n_fun,0,kin_max)
    }
  }

  
  
  ############## Simulations
  
  #For each simulation
  for (i_sim in 1:nsim) {
    
    ## Initial values
    
    #if different E0 for simulations
    if (same.E0==FALSE) {
      if (is.random.E0==FALSE) {
        E0_sim <- E_ini_fun[i_sim,]
      } else {
        #take random values
        E0_alea <- runif(n_fun,0,Etot_fun)
        #recalibrate to have sum(E0_sim) = Etot_fun
        E0_sim <- E0_alea*Etot_fun/sum(E0_alea)
      }
    }
    #if different kin0 for simulations
    if (same.kin0==FALSE) {
      if (is.random.kin0==FALSE) {
        kin0_sim <- kin_fun[i_sim,]
      } else {
        #take random values
        kin0_sim <- runif(n_fun,0,kin_max)
      }
    }
    A0_sim <- activities(kin0_sim,Keq_fun)
    
    #keep initial values
    initial_values$E0 <- rbind(initial_values$E0,E0_sim)
    initial_values$kin0 <- rbind(initial_values$kin0,kin0_sim)
    initial_values$A0 <- rbind(initial_values$A0,A0_sim)
    
    
    ## Start simulation
    simu <- simul.evol.enz.one(E0_sim, kin0_sim, Keq_fun, N_fun, correl_fun, 
                               beta_fun, X_fun, pasobs, npt, max_mut_size_E, max_mut_size_A,
                               pmutA, typ_E, typ_A, use.old.mut)
    res <- simu$res_sim
    enzstar <- simu$pred_enzsim
    contstar <- simu$pred_contsim
    
    ## Final values
    final_values$E_f <- rbind(final_values$E_f,res[npt,1:n_fun])
    final_values$kin_f <- rbind(final_values$kin_f,res[npt,(n_fun+1):(2*n_fun)])
    final_A <- activities(res[npt,(n_fun+1):(2*n_fun)],Keq_fun)
    final_values$A_f <- rbind(final_values$A_f,final_A)
    #final_value$A_f <- rbind(final_values$A_f,res[npt,(2*n+4):(3*n+3)])
    
    
    ## Keep data
    #transfo in dataframe
    res <- data.frame(res)
    #add simulation number
    res$sim <- i_sim
    #put results of this simulation under the other
    tabR <- rbind.data.frame(tabR,res)
    
    ## Idem for predictions
    enzstar <- data.frame(enzstar)
    enzstar$sim <- i_sim
    tabP_e <- rbind.data.frame(tabP_e,enzstar)
    contstar <- data.frame(contstar)
    contstar$sim <- i_sim
    tabP_r <- rbind.data.frame(tabP_r,contstar)
    
    ## Progress tracking
    #print(i_sim)
    message(paste0("Simulation number : ", i_sim))
    
    ## End of simulation
  }

  
  return(invisible(list(tabR=tabR, list_init=initial_values, list_final=final_values, tabP_e=tabP_e, tabP_r=tabP_r,
              param=list(n=n_fun,nsim=nsim,E0=initial_values$E0,kin0=initial_values$kin0,Keq=Keq_fun,beta=beta_fun,B=B_fun,
                         correl=correl_fun,N=N_fun,pasobs=pasobs,npt=npt,ngene=npt*pasobs,
                         X=X_fun,Etot0=Etot_fun,kin_max=kin_max,
                         max_mut_size_E=max_mut_size_E,max_mut_size_A=max_mut_size_A,
                         pmutA=pmutA,typmut_E=typ_E,typmut_A=typ_A,
                         same.E0=same.E0, is.random.E0=is.random.E0, same.kin0=same.kin0, is.random.kin0=is.random.kin0))))
}



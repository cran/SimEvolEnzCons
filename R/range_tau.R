#' Bounds of \emph{tau}
#'
#' Computes the bounds of the driving variable (or position) \eqn{\tau} on the straight line of relative concentrations such as thereof are between 0 and 1 
#' 
#'
#' @details 
#' Relative concentrations is defined in [0,1].
#' This function \code{range_tau} computes the bounds of \eqn{\tau} such as all relative concentrations are between 0 and 1.
#' 
#' The inferior (resp. superior) bound of \eqn{\tau} corresponds to minimal (resp. maximal) value of \eqn{\tau}
#' such as all relative concentrations are superior or equal to 0 \bold{and} inferior or equal to 1,
#' with at least one relative concentration equal to 0 or 1.
#' 
#' 
#' 
#' 
#' @usage range_tau(E_ini_fun,B_fun)
#'     
#'
#' @param E_ini_fun Numeric vector of initial concentrations
#' @param B_fun Numeric vector of global co-regulation coefficients. Same length as \code{E_ini_fun}.
#'
#'
#' @return Numeric vector of the inferior and the superior bounds of driving variable \eqn{\tau}
#' 
#' 
#' @seealso 
#' See details of function \code{\link{droite_e}} for detailed explanations on \eqn{\tau}.
#' 
#'
#' @examples
#' beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
#' B <- apply(beta,1,sumbis)
#' E0 <- c(30,30,30)
#' 
#' range_tau(E0,B)
#'
#' @export





#Bounds of tau
range_tau <- function(E_ini_fun,B_fun) {
  
  ####
  #Number of enzymes
  n_fun <- length(E_ini_fun) 
  #To avoid problem due to use of dataframe
  E_ini_fun <- as.matrix(E_ini_fun)
  # Computes initial relative concentrations
  e0_fun <- E_ini_fun/sum(E_ini_fun)
  
  
  ######### Verif parameters
  if (length(E_ini_fun)!=length(B_fun)) {
    stop("An enzyme are missing. B_fun and E_ini_fun have not the same length.")
  }
  
  ##### Special cases
  if (all(e0_fun==(1/B_fun))) {
    warning("E_ini_fun is a multiple of 1/B_fun. Line becomes a point and tau does not exist.")
  }
  
  
  
  
  
  #########" Bornes de tau
  # tau is between -e0/(1/B-e0) for (e_i=0) and (1-e0)/(1/B-e0) for (e_i=1) for all i
  
  if (any(e0_fun!=1/B_fun)) {
    #inf_t (resp. sup_t) is vector of all possible inferior (resp. superior) bounds of tau
    inf_t <- NULL
    sup_t <- NULL
    
    #for all enzymes, computation of bounds
    for (i in 1:n_fun) {
      
      #to avoid division by 0
      if (1/B_fun[i]-e0_fun[i]!=0) {
        
        # computes -e0/(1/B-e0)
        t1 <- -e0_fun[i]/(1/B_fun[i]-e0_fun[i])
        
        # computes (1-e0)/(1/B-e0)
        t2 <- (1-e0_fun[i])/(1/B_fun[i]-e0_fun[i])
        
        # if denominator is positive
        if (1/B_fun[i]-e0_fun[i]>0) {
          # t1 < tau < t2
          inf_t <- c(inf_t,t1)
          sup_t <- c(sup_t,t2)
          
          #if denominator is negative, bounds are inverted
        } else {
          # t1 > tau > t2
          inf_t <- c(inf_t,t2)
          sup_t <- c(sup_t,t1)
        }
      }
    }
    
    # choose the maximal inferior (resp. minimal superior) bound of tau
    mi_t <- max(inf_t)
    ma_t <- min(sup_t)
    
    # inferior and superior bounds of tau
    borne_fun <- c(mi_t,ma_t)
    
  }
  
  ## if e0=1/B, line does not exist
  if (all(e0_fun==(1/B_fun))) {
    #to simplify further function, put 0 for e0 and 1 for e*, but does not change droite_e results
    borne_fun <- c(0,1)
  }
  
  return(borne_fun)
}



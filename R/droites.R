#' @name droites
#' @aliases droite_e
#' @aliases droite_E.Reg
#' @aliases droite_E.CR
#' @aliases droite_tau
#' 
#' 
#' @title Line of relative enzyme concentrations (co-regulations cases)
#'
#' @description Computes the position of the point of relative concentrations on the line on which they move on in the cases of co-regulations
#' 
#' 
#' @details 
#' In the cases of co-regulations, relative enzymes concentrations evolve along a straight line.
#' This line is determined by to factors: initial enzyme concentrations and global co-regulation factors.
#' The driving variable \eqn{\tau} is a parameter indicating the position of the relative enzyme concentrations \emph{e} on this line.
#' 
#' 
#' Function \code{droite_e} gives the relative concentrations corresponding to the input position \code{tau_fun}.
#' 
#' Function \code{droite_E.Reg} gives the absolute concentrations corresponding to the input position \code{tau_fun} in case of \bold{regulation only} (constraints abbreviation \code{"RegPos"} or \code{"RegNeg"}).     
#' 
#' Function \code{droite_E.CR} gives the absolute concentrations corresponding to the input position \code{tau_fun} in case of \bold{competition and regulation} (constraints abbreviation \code{"CRPos"} or \code{"CRNeg"}).
#' 
#' Function \code{droite_tau} gives the position \eqn{\tau} corresponding to the input enzyme concentrations \code{E_fun}.
#' 
#' 
#' Note that if initial relative concentrations \code{E_ini_fun} is a multiple of \code{1/B_fun}, the line becomes a point and position \eqn{\tau} does not exist.
#' 
#' 
#'     
#' @param tau_fun Numeric value of the position of relative enzyme concentrations on the line
#' @param E_ini_fun Numeric vector of initial concentrations
#' @param B_fun Numeric vector of global co-regulation coefficients. Same length as \code{E_ini_fun}.
#' @param E_fun numeric vector of current concentrations on the line. Same length as \code{E_ini_fun}.
#' 
#'
#' @section Special results:
#' \bold{\emph{Initial point} }
#' 
#' If \code{tau_fun} is equal to 0, \code{droite_e} returns the value of initial relative concentrations, \emph{i.e.} value of \code{E_ini_fun} divided by \code{sum(E_ini_fun)};
#' \code{droite_E.Reg} and \code{droite_E.CR} returns the value of \code{E_ini_fun}.
#' 
#' If \code{E_fun} is a multiple of \code{E_ini_fun}, function \code{droite_tau} returns \code{0} \emph{(initial point)}.
#' 
#' \bold{\emph{End point} }
#' 
#' If \code{tau_fun} is equal to 1, \code{droite_e} returns the reverse value of \code{B_fun};
#' \code{droite_E.CR} returns a multiple of the reverse value of \code{B_fun};
#' \code{droite_E.Reg} returns \code{Inf}.
#' 
#' If \code{E_fun} is a multiple of \code{1/B_fun}, function \code{droite_tau} returns \code{1} \emph{(end point)}.
#' 
#' \bold{\emph{Line becomes a point} }
#' 
#' If \code{E_ini_fun} is a multiple of \code{1/B_fun},\code{droite_e} returns the value of \code{1/B_fun};
#' \code{droite_E.Reg} returns an error, because concentrations \code{E} can be any multiple of \code{1/B_fun} without variation of relative concentrations;
#' \code{droite_E.CR} returns \code{E_ini_fun};
#' \code{droite_tau} returns an error because \eqn{\tau} does not exist in this case.
#' 
#' \bold{\emph{Line does not exist} }
#' 
#' If there only one enzyme (length of \code{E_ini_fun} is equal to 1), relative concentrations is always equal to 1.
#' \code{droite_e} should return 1;
#' \code{droite_E.Reg}, \code{droite_E.CR} and \code{droite_tau} should return an error.
#'
#' 
#' @seealso 
#' To compute global co-regulation coefficients \code{B_fun} from co-regulation matrix \code{beta_fun}, see the example or use function \code{\link{compute.B.from.beta}}.
#'
#'
#' @examples
#' beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
#' B <- apply(beta,1,sumbis)
#' E0 <- c(30,30,30)
#' tau <- 0.5
#' 
#' droite_e(tau,E0,B)
#' 
#' E <- droite_E.Reg(tau,E0,B)
#' droite_tau(E,E0,B)
#' 
#' E <- droite_E.CR(tau,E0,B)
#' droite_tau(E,E0,B)
#'
#' 
#' 
#' 
NULL





#' @rdname droites
#' @usage droite_e(tau_fun,E_ini_fun,B_fun)
#' @return \code{droite_e} returns a numeric vector of relative concentrations
#' @export

######### Line to deduce initial concentration from position (cases Reg et CR)
droite_e <- function(tau_fun,E_ini_fun,B_fun) {
  #computes initial relative concentrations
  e0_fun <- E_ini_fun/sum(E_ini_fun)
  
  ### Verif parameters
  if (length(E_ini_fun)!=length(B_fun)) {
    stop("An enzyme are missing. B_fun and E_ini_fun have not the same length.")
  }
  
  ####
  #if there is only one enzyme
  if (length(E_ini_fun)==1) {
    #relative concentration is equal to 1
    e_fun <- 1
  } else {
    #deduce current relative concentrations from line expression
  e_fun <- tau_fun*(1/B_fun-e0_fun) + e0_fun
  return(e_fun)
  }
  
}





#' @rdname droites
#' @usage droite_E.Reg(tau_fun,E_ini_fun,B_fun)
#' @return \code{droite_E.Reg} returns a numeric vector of absolute concentrations
#' @export

############ Line to deduce absolute concentration from position (case Reg)
droite_E.Reg <- function(tau_fun,E_ini_fun,B_fun) {
  ### Verif parameters
  if (length(E_ini_fun)!=length(B_fun)) {
    stop("An enzyme are missing. B_fun and E_ini_fun have not the same length.")
  }
  # #computes initial relative concentrations
  # e0_fun <- E_ini_fun/sum(E_ini_fun)
  # if (all(e0_fun==(1/B_fun))) {
  #   stop("E_ini_fun is a multiple of 1/B_fun. Line becomes a point and tau does not exist.")
  # }
  
  ####
  #compute initial total concentration
  Etot_fun <- sum(E_ini_fun)
  #deduce current absolute concentrations from line expression
  E_fun <- E_ini_fun-Etot_fun/B_fun+Etot_fun/((1-tau_fun)*B_fun)
  return(E_fun)
}





#' @rdname droites
#' @usage droite_E.CR(tau_fun,E_ini_fun,B_fun)
#' @return \code{droite_E.CR} returns a numeric vector of absolute concentrations
#' @export

############" Line to deduce absolute concentration from position (case CR)
droite_E.CR <- function(tau_fun,E_ini_fun,B_fun) {
  ### Verif parameters
  if (length(E_ini_fun)!=length(B_fun)) {
    stop("An enzyme are missing. B_fun and E_ini_fun have not the same length.")
  }
  
  ####
  #deduce current absolute concentrations from line expression
  E_fun <- tau_fun*(sum(E_ini_fun)/B_fun-E_ini_fun)+E_ini_fun
  return(E_fun)
}




#' @rdname droites
#' @usage droite_tau(E_fun,E_ini_fun,B_fun)
#' @return \code{droite_tau} returns a numeric value giving the position on the line
#' @export

############ Line to deduce position tau from current absolute concentrations
droite_tau <- function(E_fun,E_ini_fun,B_fun) {
  
  #####
  #Number of enzymes
  n_fun <- length(E_ini_fun)
  #computes initial relative concentrations
  e0_fun <- E_ini_fun/sum(E_ini_fun)
  #computes current relative concentrations
  e_fun <- E_fun/sum(E_fun)
  
  
  ### Verif parameters
  if (length(E_ini_fun)!=length(E_fun)) {
    stop("An enzyme have disapeared. E_fun and E_ini_fun have not the same length.")
  }
  if (length(E_ini_fun)!=length(B_fun)) {
    stop("An enzyme are missing. B_fun and E_ini_fun have not the same length.")
  }
  if (all(e0_fun==(1/B_fun))) {
    stop("E_ini_fun is a multiple of 1/B_fun. Line becomes a point and tau does not exist.")
  }

  
  ####
  
  #deduces position tau from line expression
  tau_fun <- (e_fun - e0_fun)/(1/B_fun - e0_fun)
  #to avoid problem with data.frame in further functions and all.equal()
  tau_fun <- as.numeric(tau_fun)
  
  ###Test if tau is different between enzymes
  #Use of all.equal because of small differences
  #replication n times (=nb enzymes) of the first value of tau to test against all values of tau
  #if class is logical, isTRUE(all.equal()) is TRUE, else, all.equal() is "character" indicating mean differences between tested vectors
  if (class(all.equal(rep(tau_fun[1],n_fun),tau_fun))=="character") {
    stop("E_fun is not on the line.")
  }
  
  #only first value is preserved, because all values are identical
  return(tau_fun[1])
}


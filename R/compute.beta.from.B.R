#' Co-regulation coefficient computation
#'
#' Computes the matrix of co-regulation coefficients \eqn{M_\beta} from a vector of global co-regulation coefficients \code{B}
#' 
#'
#' @details 
#' Default value \code{L_Phi_fun = NULL} is appropriate only if enzymes are supposed to be all co-regulated, i.e. no \eqn{\beta} value is null.
#' Return the same result as as a list of one element such as \code{L_Phi_fun = list(1:n)}, where \eqn{n} is the number of enzymes, which is also the length of \code{B_fun}.
#' 
#' \code{L_Phi_fun} must be a list of \eqn{p} elements (the number of regulation groups) containing numbers between 1 and \code{n} (number of enzymes), where each list element contains the numbers of the enzyme in the group.
#' Each enzyme only occurs in one group.
#' See function \code{\link{class_group}} to have an idea of \code{L_Phi_fun} structure.
#' For \code{compute.beta.from.B}, independent enzymes can be not contained in \code{L_Phi_fun}, and thus, \code{L_Phi_fun} may be smaller than \code{B_fun}.
#' 
#' 
#' @usage compute.beta.from.B(B_fun,L_Phi_fun=NULL)
#'     
#' 
#' @param B_fun Numeric vector of global co-regulation coefficients
#' @param L_Phi_fun List of regulation groups. Default is \code{NULL}.
#' 
#' 
#'
#'
#' @return Numeric matrix \code{n*n} of the co-regulation coefficients, where \code{n} is the length of \code{B_fun}.
#' 
#' If \code{beta_fun} is \code{NULL}, \code{compute.beta.from.B} returns \code{NULL}.
#' 
#' 
#' @seealso 
#' Use function \code{\link{is.B.accurate}} to verify \code{B_fun} conformity.
#' 
#' See function \code{\link{class_group}} to know format of \code{L_Phi_fun}.
#'
#'
#'
#' @examples
#' B <- 1/c(0.5,0.2,0.3)
#' correl <- "RegPos"
#' 
#' is.B.accurate(B,3,correl)
#' 
#' beta <- compute.beta.from.B(B)
#' 
#' #Seven enzymes and three groups
#' n <- 7
#' p <- 3
#' B <- c(1.1824,  3.695, -8.593023,  1.3, 13,  6.5, 1)
#' L_Phi <- list(1:3) #firt three enzymes in first group
#' L_Phi[[2]] <- 4:6
#' L_Phi[[3]] <- 7 #last enzyme independent
#' 
#' beta <- compute.beta.from.B(B,L_Phi)
#' 
#' 
#'
#' @export




compute.beta.from.B <- function(B_fun,L_Phi_fun=NULL) {
  
  #add is.B.accurate() ?
  
  # if beta_fun is NULL
  if (length(B_fun)==0) {
    beta_fun <- NULL
  }
  
  # in other cases
  if (length(B_fun)>1) {
    #B*t(1/B)
    beta_fun <- B_fun %*% t(1/B_fun)
  }
  
  if (length(L_Phi_fun)>0) {
    #diagonal matrix of same length as B
    beta_fun <- diag(1,length(B_fun))
    #number of regulation group
    p_fun <- length(L_Phi_fun)
    #check number of enzymes = same totale elements in B and L_Phi
    # if (length(B_fun)!=length(unlist(L_Phi_fun))) {
    #   stop("An enzyme is missing. B_fun and L_Phi_fun do not have the same enzyme number.")
    #   #not necessary if we consider that missing enzymes are independent
    # }
    for(q in 1:p_fun) {
      #vector of enzymes in the group
      phi_q <- unlist(L_Phi_fun[[q]])
      #compute matrix of beta ONLY for enzymes of the group (= regulation submatrix of the group)
      beta_fun[phi_q,phi_q] <- B_fun[phi_q] %*% t(1/B_fun[phi_q])
    }
  }
  
  return(beta_fun)
}


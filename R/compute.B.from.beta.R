#' Global co-regulation coefficient computation
#'
#' Computes the global co-regulation coefficients \code{B} from a matrix of co-regulation coefficients \code{beta}
#' 
#'
#' @details 
#' 
#' \code{beta_fun} have same number of rows and columns.
#' 
#' 
#' @usage compute.B.from.beta(beta_fun)
#'     
#' 
#' @param beta_fun Matrix of co-regulation coefficients
#' 
#' 
#'
#' @return Numeric vector of the \code{n} global co-regulation coefficients.
#' 
#' If \code{beta_fun} is \code{NULL}, \code{compute.B.from.beta} returns \code{NULL}.
#' 
#' @seealso 
#' Use function \code{\link{is.beta.accurate}} to verify \code{beta_fun} conformity.
#'
#'
#'
#' @examples
#' beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
#' correl <- "RegPos"
#' 
#' is.beta.accurate(beta,3,correl)
#' 
#' B <- compute.B.from.beta(beta)
#' 
#'
#' @export




compute.B.from.beta <- function(beta_fun) {
  
  #add is.beta.accurate() ?
  
  # if beta_fun is NULL
  if (length(beta_fun)==0) {
    B_fun <- NULL
  }
  
  # if only one enzyme
  if (length(beta_fun)==1) {
    B_fun <- 1
  }
  
  # in other cases
  if (length(beta_fun)>1) {
    #sum by row
    B_fun <- apply(beta_fun,1,sumbis)
  }
  
  return(B_fun)
}


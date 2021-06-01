#' Help to name parameter 'correl'
#'
#' Give the correct abbreviation of the applied constraint, used for parameter \code{correl_fun}
#'  
#'
#' @details 
#' Explored constraints are competition and/or regulation.
#' 
#' If you choose to put regulation, matrix of co-regulation coefficient \code{beta} \bold{or} vector of global co-regulation coefficients \code{B} are needed.
#' 
#' See function \code{\link{is.correl.authorized}} to know possible value of parameter \code{correl_fun}.
#' 
#' 
#' 
#' @usage name.correl(is.comp,is.reg,beta_fun=NULL)
#'     
#' @param is.comp Logical. Is there competition for resources ?
#' @param is.reg Logical. Is there regulation between enzymes ?
#' @param beta_fun Numeric matrix of co-regulation coefficients. Default is NULL, but needed if there is regulation (\code{is.reg==TRUE}).
#' 
#'
#'
#' @return A character string, used for parameter \code{correl_fun}. See possible values in \code{\link{is.correl.authorized}}.
#' 
#'
#'
#' @examples
#' 
#' #Independence 
#' name.correl(is.comp=FALSE,is.reg=FALSE)
#' #returns "SC"
#' 
#' #Regulation
#' beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
#' name.correl(is.comp=FALSE,is.reg=TRUE,beta_fun=beta)
#' #returns "RegPos"
#' 
#' 
#'
#' @export




name.correl <- function(is.comp,is.reg,beta_fun=NULL) {
  
  # Independence
  if (is.comp==FALSE&is.reg==FALSE) {
    correl_fun <- "SC"
  }
  
  # Competition only
  if (is.comp==TRUE&is.reg==FALSE) {
    correl_fun <- "Comp"
  }
  
  # Regulation
  if (is.reg==TRUE) {
    
    # verif if a parameter is missing
    if (length(beta_fun)==0) {
      stop("Co-regulation coefficients 'beta_fun' or 'B_fun' are required.")
    }
    
    
    # Kind of regulation
    if (all(beta_fun>=0)==TRUE) {
      # positive regulation
      which.reg <- "Pos"
    } else {
      #negative regulation
      which.reg <- "Neg"
    }
    
    
    # and competition ?
    if (is.comp==FALSE) {
      # regulation only
      correl_fun <- paste("Reg",which.reg,sep="")
    } else {
      # with competition
      correl_fun <- paste("CR",which.reg,sep="")
    }
    
  }

  return(correl_fun)
}

  
  
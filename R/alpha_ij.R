#' Redistribution coefficient computation
#'
#' Computes the redistribution coefficients between all enzymes
#' 
#'
#' 
#' @references 
#' Lion, S., F. Gabriel, B. Bost, J. Fiévet, C. Dillmann, and D. De Vienne, 2004.
#' An extension to the metabolic control theory taking into account correlations between enzyme concentrations.
#' European Journal of Biochemistry 271:4375–4391.
#' 
#' 
#' @usage alpha_ij(E_fun,correl_fun,beta_fun=NULL)
#'     
#' @param E_fun Numeric vector of concentrations
#' @inheritParams is.correl.authorized
#' @param beta_fun Matrix of co-regulation coefficients
#' 
#' 
#' 
#' @details 
#' The redistribution coefficients are influenced by the applied constraint.
#'     
#' For further details on how the redistribution coefficient is calculated, see Lion \emph{et al.} (2004).
#' 
#' For cases with co-regulations (i.e. \code{correl_fun} value is \code{"RegPos"}, \code{"RegNeg"}, \code{"CRPos"} or \code{"CRNeg"}), \code{beta_fun} is obligatory. In other cases (i.e. \code{correl_fun} value is \code{"SC"} or \code{"Comp"}), \code{beta_fun} is ignored, that is why default is \code{NULL}.
#' 
#' When \code{beta_fun} is obligatory, it has to be a square matrix of size \code{n}, where this number \code{n} is the length of \code{E_fun}.
#' 
#'
#'
#' @return Square matrix of size \code{n} (where \code{n} is the length of \code{E_fun}, which is the number of enzymes) of the redistribution coefficients between all enzymes.
#' 
#' 
#' @seealso 
#' See function \code{\link{is.correl.authorized}} to have more information about constraints and on the usage of argument \code{correl_fun}.
#'
#' @examples
#' E <- c(30,30,30)
#' beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
#' correl <- "SC"
#' 
#' alph <- alpha_ij(E,correl,beta)
#' 
#'
#' @export


# Calcul du coef de redistri
alpha_ij <- function(E_fun,correl_fun,beta_fun=NULL) {
  # Number of enzymes
  n_fun <- length(E_fun)
  
  ##### Verif parameters
    # cases allowed; if not, program stops
  is.correl.authorized(correl_fun)
  
    # verif beta accuracy
  is.beta.accurate(beta_fun,n_fun,correl_fun)
  # if (nrow(beta_fun)!=n_fun|ncol(beta_fun)!=n_fun) {
  #   stop("Dimensions of beta_fun matrix is not equal to E_fun length.")
  # }

  
  
  ### Need further
  #Computation of global co-regulation coefficients in case of regulation
  if (correl_fun!="SC"&correl_fun!="Comp") {
     B_fun <- apply(beta_fun,1,sumbis)
  }
 
  
  
  #### Computation of redistribution coefficients
  Alpha_fun <- diag(rep(1,n_fun))
  #Alpha_fun <- matrix(rep(0,n_fun*n_fun),nrow=n_fun)
  #diag(Alpha_fun) <- 1
  # for (i_fun in 1:n_fun) {
  #   Alpha_fun[i_fun,i_fun] <- 1
  # }
  
  if (correl_fun =="Comp") {
    for (i_fun in 1:n_fun) {
      for (j_fun in 1:n_fun) {
        if (j_fun!=i_fun) {
          Alpha_fun[i_fun,j_fun] <- -E_fun[j_fun]/(sum(E_fun)-E_fun[i_fun])
        }
      }
    }
  }
  
  if (correl_fun == "RegNeg"|correl_fun=="RegPos") {
    Alpha_fun <- beta_fun
  }
  
  if (correl_fun =="CRPos"|correl_fun=="CRNeg") {
    for (i_fun in 1:n_fun) {
      #to avoid division by zero
      if (E_fun[i_fun]/sum(E_fun) != 1/B_fun[i_fun]) {
        for (j_fun in 1:n_fun) {
          if (j_fun != i_fun) {
            Alpha_fun[i_fun,j_fun] <- (beta_fun[i_fun,j_fun]*sum(E_fun)-B_fun[i_fun]*E_fun[j_fun])/(sum(E_fun)-B_fun[i_fun]*E_fun[i_fun])
          }
        }
      }
    }
  }
  return(Alpha_fun)
}

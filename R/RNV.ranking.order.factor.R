#' Name and value of RNV-ranking-order factor
#'
#' Gives the name and values of the RNV-ranking-order factor
#' 
#'
#' @details 
#' \bold{Factors governing ranking order of RNV}
#' 
#' RNV-ranking-order factor depends on constraint. It would be:
#' \itemize{
#'    \item in case of independence (\code{correl_fun="SC"}): activities \code{A} power 1/3;
#'    \item in case of competition only (\code{correl_fun="Comp"}): activities \code{A} power 1/4;
#'    \item in case of regulation only (\code{correl_fun="RegPos"} or \code{"RegNeg"}): absolute value of inverse of global co-regulation coefficients \code{|1/B|};
#'    \item in case of competition and regulation (\code{correl_fun="CRPos"} or \code{"CRNeg"}): absolute value of inverse of global co-regulation coefficients minus initial relative concentrations \code{|1/B - e0|}.
#' }
#' 
#' 
#' 
#' 
#' 
#' @usage RNV.ranking.order.factor(A_fun,correl_fun,E_ini_fun,B_fun=NULL)
#'
#' @inheritParams predict_th
#' @param E_ini_fun Numeric vector of initial concentrations.
#' 
#'
#'
#' @return Invisible list of 2 elements:
#' \itemize{
#'    \item \code{$value}: numeric vector (length \code{n}) of the values of the RNV-ranking-order factor for each enzymes
#'    \item \code{$name}: character sting indicating the name of the RNV-ranking-order variable
#' }
#' 
#' 
#' @examples 
#' A <- c(1,20,30)
#' E0 <- c(30,30,30)
#' 
#' #Independence
#' rank_var <- RNV.ranking.order.factor(A,"SC",E0)
#' all.equal(A^(1/3),rank_var$value) #TRUE
#' 
#' #Positive regulation
#' B <- 1/c(0.2,0.5,0.3)
#' rank_var <- RNV.ranking.order.factor(A,"RegPos",E0,B)
#' all.equal(1/B,rank_var$value) #TRUE
#'
#'
#' @export



RNV.ranking.order.factor <- function(A_fun,correl_fun,E_ini_fun,B_fun=NULL) {
  
  ### Config
  #Number of enzymes
  n_fun <- length(A_fun)
  #To avoid problem due to use of dataframe
  E_ini_fun <- as.matrix(E_ini_fun)
  # Computes initial relative concentrations
  e0_fun <- E_ini_fun/sum(E_ini_fun)
  
  ### Tests (see function predict.eff)
  is.correl.authorized(correl_fun)
  is.B.accurate(B_fun,n_fun,correl_fun)
  if (length(E_ini_fun)!=n_fun) {
    stop("An enzyme are missing. A_fun and E_ini_fun have not the same length.")
  }
  
  
  ### Name and values of RNV-ranking-order variable
  
  #rank_var = value of RNV-ranking-order factor, vector of length n
  #rank_var_name = name of this factor, character
  
  #SC
  if (correl_fun=="SC") {
    rank_var <- A_fun^(1/3)
    rank_var_name <- "Activity^(1/3)"
  }
  
  #Comp
  if (correl_fun=="Comp") {
    rank_var <- A_fun^(-1/4)
    rank_var_name <- "Activity^(-1/4)"
  }
  
  #RegPos and RegNeg
  if (correl_fun=="RegPos"|correl_fun=="RegNeg") {
    rank_var <- abs(1/B_fun)
    rank_var_name <- "Absolute value of 1/B"
  }
  
  #CrPos and CRNeg
  if (correl_fun=="CRPos"|correl_fun=="CRNeg") {
    
    
    
    rank_var <- abs(1/B_fun-e0_fun)
    rank_var_name <- "Absolute value of (1/B - e0)"
  }
  
  
  return(list(value=rank_var,name=rank_var_name))
}




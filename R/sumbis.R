#' Sum of vector elements
#' 
#' \code{sumbis} computes the sum of a vector
#' 
#' 
#' @details 
#' \code{sumbis} computes the sum of a vector,
#'     to avoid some problems in comparison test due to R \code{\link{sum}} function.
#' 
#' This problem can also be avoided by using \code{all.equal} rather than \code{all} for example.
#'
#' \code{sumbis} is used for compute global co-regulation coefficients in \code{\link{compute.B.from.beta}}.
#'
#' @usage sumbis(v)
#'
#' @param v Numeric vector
#'
#' @return Numeric corresponding to the sum of all values present in the input vector
#' 
#' 
#' @examples
#'   # Computation of co-regulation coefficients
#' #Number of enzymes
#' n <- 3
#' #Matrix of co-regulations (n column and n rows)
#' beta<-matrix(rep(0,n*n),nrow=n)
#' #Vector of global co-regulation coefficients (length n)
#' B <- rep(0,n)
#' 
#' #For each enzyme
#' for (j in 1:n){
#'   beta[j,j]<-1
#' }
#' 
#' #Set of co-regulation coefficients
#' beta[1,2]<- 0.1
#' beta[2,3]<- 2.0 #-0.43
#' beta[1,3]<- beta[1,2]*beta[2,3]
#' beta[2,1]<- 1/beta[1,2]
#' beta[3,1]<- 1/beta[1,3]
#' beta[3,2]<- 1/beta[2,3]
#' 
#' #Computation of global co-regulation coefficients
#' for (j in 1:n){
#'   B[j]<-sumbis(beta[j,])
#' }
#' #or
#' apply(beta,1,sumbis)
#' 
#' # result : B = c(1.3, 13.0,  6.5)
#'
#' @export


sumbis=function(v) {
  S=0
  for (i in 1:length(v)){
    #add to previous sum the current element of vector v
    S=S+v[i]
  }
  return(S)
}

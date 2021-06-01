#' Activity computation
#'
#' Computes the enzyme pseudo-activities
#' 
#'
#' @details 
#' Computes the pseudo-activities of enzymes,
#'     i.e. for each enzyme, the product of its kcat/Km (kinetic parameters)
#'     by the product of the upstream reactions equilibrium constants.
#'     
#' In other functions, pseudo-activity is also named activity.
#' 
#' \code{kin_fun} and \code{Keq_fun} need to have the same length.
#' 
#' @references Lion, S., F. Gabriel, B. Bost, J. Fiévet, C. Dillmann, and D. De Vienne, 2004. 
#' An extension to the metabolic control theory taking into account correlations between enzyme concentrations.
#' European Journal of Biochemistry 271:4375–4391.
#' 
#' 
#' @usage activities(kin_fun,Keq_fun)
#'     
#'
#' @param kin_fun Numeric vector of kinetic parameters (catalytic constant kcat divided by Michaelis constant Km)
#' @param Keq_fun Numeric vector of equilibrium constants
#'
#'
#' @return Numeric vector of activities of each enzyme, of same length as \code{kin_fun} and \code{Keq_fun}.
#' 
#'
#' @examples
#' #Values from CoA metabolism
#' kin <- c(53/0.29,50/0.78,29,6.22) #kinetic parameters kcat/Km
#' Keq <- c(1.1e+8,4.9e+3,1.1e+3,0.228) #equilibrium constants
#' A <- activities(kin,Keq) #activities
#' 
#' # results : A = c(1.827586e+02, 7.051282e+09, 1.563100e+13, 3.687838e+15)
#'
#' @export


#Activities computation = upstream Keq x kcat/KM
activities <- function(kin_fun,Keq_fun){
  #number of enzymes
  n_fun <- length(kin_fun)
  
  if (length(Keq_fun)!=n_fun) {
    stop("You forgot some reactions. kin_fun and Keq_fun need to have the same length")
  }
  
  #Keq_upstream is the vector of product of upstream Keq for each enzymes
  Keq_upstream <- NULL #valeur neutre de la multiplication, car la premiere reaction de la chaine n'a pas d'amont
  
  #up_fun is a numeric value used to compute product upstream Keq
  #up_fun is set to 1, neutral value of multiplication, because there is no upstream Keq to the first enzyme
  up_fun <- 1
  
  #for each enzyme
  for (i_fun in 1:n_fun) {
    
    #add the current value of product of upstream Keq
    Keq_upstream <- c(Keq_upstream,up_fun)
    
    #multiply the Keq of current reaction to the upstream Keq product of this reaction => it becomes the upstream product to the next enzyme
    up_fun <- up_fun * Keq_fun[i_fun]
  }
  
  #activities = upstream Keq x kcat/KM
  A_fun <- Keq_upstream*kin_fun
  
  return(A_fun)
}


# 
# #Activities computation = upstream Keq x kcat/KM
# activities_bis <- function(kin_fun,Keq_fun){
#   #number of enzymes
#   n_fun <- length(kin_fun)
#   
#   #Keq_upstream is the vector of product of upstream Keq for each enzymes
#   #first 
#   Keq_upstream <- c(1) #valeur neutre de la multiplication, car la premiere reaction de la chaine n'a pas d'amont
#   up_fun <- 1 #valeur neutre
#   for (i_fun in 2:n_fun) { #première valeur n'existe pas, car pas d'amont
#     up_fun <- up_fun * Keq_fun[i_fun-1] #on prend la valeur de la reaction juste avant i et on l' "ajoute" aux precedentes
#     Keq_upstream <- c(Keq_upstream,up_fun)
#   }
#   A_fun <- Keq_upstream*kin_fun #Keq en amont x kcat/KM
#   return(A_fun)
# }



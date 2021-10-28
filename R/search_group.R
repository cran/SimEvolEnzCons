#' Search the regulation group for an enzyme
#'
#' Give the number of the regulation group \eqn{Phi_q} where is the interest enzyme giving the list of regulation group
#'  
#'
#' @details 
#' Enzymes are classified in regulation groups depending on the co-regulation matrix (see function \code{link{class_group}}).
#'
#' Function \code{search_group} also allows to find the group type from the list of group types.
#' In that case, \code{i} is the interest group, and \code{Lv} is the list of group type, which the output of function \code{\link{group_types}}.
#' 
#' More largely, the function gives the number of a list element that contains a particular number. 
#' 
#' @usage search_group(i,Lv)
#'     
#' @param i Integer which is the number of the interest enzyme
#' @param Lv List of regulation group, preferably the output of function \code{\link{class_group}}
#' 
#'
#' @return Return an integer which is the number \eqn{q} of the group \eqn{Phi_q} that contains the interest enzyme
#' 
#' @seealso 
#' Function \code{\link{class_group}} to compute the list of regulation groups
#' 
#'
#' @examples
#' 
#' ## One group
#' beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
#' L_Phi <- class_group(beta)
#' 
#' #enzyme 2 is in group 1
#' search_group(2,L_Phi)
#' 
#' ## Two groups
#' n <- 3
#' beta <- diag(1,n) 
#' beta[1,2] <- -0.32 
#' beta[2,1] <- 1/beta[1,2]
#' 
#' L_Phi <- class_group(beta)
#' 
#' search_group(2,L_Phi) #gives 1 
#' search_group(3,L_Phi) #gives 2
#' 
#' @export


search_group <- function(i,Lv) {
  ## Verif parameters
  if (i != round(i)) {
    stop("Parameter i needs to be an integer.")
  }
  
  #search if i is present in Lv, and gives the EXACT position of i in the complete list
  #LIST of same size as Lv that indicates if i is in each possible position of Lv (true if it is i, false if not)
  Lpresent <- lapply(Lv, function(x) x%in%i)
  
  #search where is i, and indicates the group that contains i
  #numeric VECTOR of same length as Lv, that indicates the element of Lv where is i (true if element contains i, false if not)
  Vpresent <- unlist(lapply(Lpresent,any))
  
  #search which element of Vpresent is true, and gives the position of this element
  #number of the group that contains i
  ki <- which(Vpresent)
  
  return(ki)
}

#fonction qui renvoie le numéro du noyau groupe appartient l'enzyme i
  #i=numéro de l'enzyme
  #Lv = liste des vecteurs de snoyaux contenant les numéros des enzymes appartenant à chaque noyau
  #Lpresent = liste de même dimension que Lv, indiquant si i se trouve dans un élément de Lv ou non
  #Vpresent = vecteur de même longueur que Lv indiquant dans quel vecteur de Lv se touve i
  #ki = numéro du vecteur contenant i


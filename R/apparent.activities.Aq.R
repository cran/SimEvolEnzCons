#' Apparent activity computation
#'
#' Computes the group apparent pseudo-activities at equilibrium
#' 
#'
#' @details 
#' Computes the apparent pseudo-activities of regulation groups at theoretical equilibrium (no competition) or effective equilibrium (competition).
#'  
#' When there is competition, \code{eiq_eff} is needed.
#' 
#' If there is only one regulation group (all enzymes are co-regulated), apparent activity of the group is not useful.
#' When the enzyme is in a singleton (enzyme independent from all other enzymes), the apparent activity of the group is identical to the activity of the enzyme.
#' 
#'     
#' In other functions, pseudo-activity is also named activity.
#' 
#' 
#' 
#' @usage apparent.activities.Aq(phi_q,A_fun,B_fun,correl_fun,eiq_eff=NULL)
#'     
#'
#' @inheritParams predict_eff
#' @param phi_q Numeric vector containing the numbers of enzymes in the group (ideally, element of the list returns by \code{\link{class_group}})
#' @param eiq_eff Numeric vector of the intra-group relative concentrations at effective equilibrium. 
#' 
#'
#' @return Numeric value of the apparent activity of chosen group.
#' 
#' @seealso 
#' Function \code{\link{class_group}} to classify enzymes in groups.
#' 
#'
#' @examples
#' #Two groups
#' n <- 3
#' A <- c(1,10,30)
#' beta <- diag(1,n) 
#' beta[1,2] <- -0.32 
#' beta[2,1] <- 1/beta[1,2]
#' B <- compute.B.from.beta(beta)
#' 
#' L_Phi <- class_group(beta)
#' 
#' correl <- "RegNeg"
#' #apply function of all groups
#' A_app <- unlist(lapply(L_Phi,apparent.activities.Aq,A,B,correl)) 
#'
#' @export


#fonction calculant l'activité apparente d'un vecteur contenant le numéro de chaque enzyme du noyau de rég
apparent.activities.Aq <- function(phi_q,A_fun,B_fun,correl_fun,eiq_eff=NULL) {
  #number of enzymes
  n_fun <- length(A_fun)
  
  #verif parameters
  is.correl.authorized(correl_fun)
  is.B.accurate(B_fun,n_fun,correl_fun)
  
  
  #### Without competition
  if (any(correl_fun==c("SC","RegPos","RegNeg"))) {
    #calculation matrix for all enzymes of the group
    Mcal <- matrix(0,nrow=n_fun,ncol=n_fun)
    #double-sum : for each enzyme in the group
    for (i in phi_q) {
      for (j in  phi_q) {
        Mcal[i,j] <- B_fun[i]^2*B_fun[j]/A_fun[j]
      }
    }
    #reverse of double-sum on i and j
    Aq_fun <- 1/sum(Mcal)
  }
  
  
  
  #### With competition
  if (any(correl_fun==c("Comp","CRPos","CRNeg"))) {
    #verif presence of eiq_eff
    if (length(eiq_eff)!=n_fun) {
      stop("An enzyme is missing in eiq_eff.")
    }
    
    #sum on enzymes in the group
    Aq_fun <- 1/sum(1/(A_fun[phi_q] * B_fun[phi_q] * (eiq_eff[phi_q])^2 ))
  }
  
  
  return(Aq_fun)
}

#phi_q = vecteur contenant le numéro des enzymes du noyau de régulation q, de longueur m_q
    #Afun = vecteur des activités pour la function, de longueur n
    #Bfun = vecteur des B pour la function, de longueur n
    #Mcal = matrice de calcul pour toutes les enzymes, de taille n*n
    #Msum = somme de la matrice de calcul



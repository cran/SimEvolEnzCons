#' Flux shape computing from one point
#'
#' \code{flux.shape.from.one.point} computes flux in every dimension (corresponding to each enzyme) from a given point (vector of concentrations)
#' 
#'  
#' 
#'
#' @details 
#' Every enzyme correspond to one dimension in a \emph{n}-dimensional graph.
#' 
#' From a given resident point \code{E_fun}, each value on dimension \emph{i} is considered as a possible mutant of enzyme concentration \eqn{E_i}.
#' Every "mutants" are taken between 0 and \code{Etot_fun} by 0.01 step.
#' In every dimension, function \code{flux.shape.from.one.point} computes flux and selection coefficient (discrete \code{\link{coef_sel.discrete}} and continuous \code{\link{coef_sel.continue}}) from this point.
#' 
#' \code{E_fun} (resp. \code{E_ini_fun}) is rescaled by a cross product to have sum of \code{E_fun} (resp. \code{E_ini_fun}) equal to \code{Etot_fun}.
#' 
#' If \code{from.eq=TRUE}, analyzed point is:
#' \itemize{
#'    \item the theoretical equilibrium in case of independence \code{"SC"} and competition only \code{"Comp"};
#'    \item near the theoretical equilibrium (tau=0.95) in case of positive regulation only \code{"RegPos"} (due to infinite possible values for concentrations at this point);
#'    \item the effective equilibrium in other cases.
#' } 
#' 
#' Default of \code{E_ini_fun} is \code{NULL} and corresponds to \code{correl_fun} equal to \code{"SC"} or \code{"Comp"},
#' but in other cases (due to presence of regulation, \code{E_ini_fun} is obligatory and needs to have the same length as \code{A_fun}).  
#' 
#' @importFrom RColorBrewer brewer.pal		
#' 
#' @usage flux.shape.from.one.point(Etot_fun, A_fun, correl_fun, beta_fun=NULL, 
#' E_ini_fun=NULL, from.eq=TRUE, E_fun=NULL, X_fun=1, with.alpha=FALSE, grp.reg=FALSE)
#'     
#' @param Etot_fun Numeric. The total concentration
#' @inheritParams coef_sel.continue
#' @param E_ini_fun Numeric vector corresponding to initial concentrations.
#' @param from.eq Logical. Is the analyzed point is the equilibrium point ?
#' If \code{TRUE}, flux and selection coefficients are computed from the equilibrium, else \code{E_fun} is required.
#' Default is \code{TRUE}.
#' @param E_fun Numeric vector of the concentrations at analysed point.
#' If \code{from.eq=TRUE}, \code{E_fun} is ignored.
#' @inheritParams flux
#' @param with.alpha Logical. For case \code{CR}, is the computing method use \code{alpha} formula (TRUE) or pass by \code{tau} computing (FALSE) ?
#' @param grp.reg Logical. Is there is some regulation groups in beta matrix ? If \code{TRUE}, \code{tau} will not be computed and give 0.
#'
#'
#' @return Invisible list of 6 elements:
#' \itemize{
#'    \item \code{$x} : Numeric vector of all values that mutated enzymes can take, between 0 and \code{Etot_fun}, by 0.01.
#'    Length of \code{(Etot_fun-0)*100}.
#'    \item \code{$J} : Numeric matrix of \code{n} columns and \code{(Etot_fun-0)*100} rows.
#'    Each column correspond to one direction (\emph{i.e. which enzyme is "mutated"}) and each row to each value of flux in this direction.
#'    \item \code{$sel_disc} : Numeric matrix corresponding to discrete selection coefficient.
#'    Same properties.
#'    \item \code{$sel_cont} : Numeric matrix corresponding to continuous selection coefficient.
#'    Same properties.
#'    \item \code{$tau} : Numeric matrix corresponding to position \eqn{\tau} in case of regulation.
#'    Same properties.
#'    \item \code{$param} : List of input parameters
#'    }
#' 
#'
#'
#' @seealso 
#' To understand differences between discrete and continuous selection coefficients, see function \code{\link{coef_sel.discrete}} and \code{\link{coef_sel.continue}}.
#' 
#'   
#' @examples
#' \donttest{
#' fsfop <- flux.shape.from.one.point(100,c(1,10,30),"SC")
#' }
#'
#' @export







flux.shape.from.one.point <- function(Etot_fun, A_fun, correl_fun, beta_fun=NULL, E_ini_fun=NULL, from.eq=TRUE, E_fun=NULL, X_fun=1, with.alpha=FALSE, grp.reg=FALSE) {
  
  ###### Setting
  n_fun <- length(A_fun)
  if (correl_fun=="SC"|correl_fun=="Comp") {beta_fun <- diag(1,n_fun)}
  B_fun <- compute.B.from.beta(beta_fun)
  
  ##### Test
  is.correl.authorized(correl_fun)
  is.beta.accurate(beta_fun,n_fun,correl_fun)
  #Regulation : need E0
  if (length(E_ini_fun)==0&correl_fun!="SC"&correl_fun!="Comp") {
    stop("In case of regulation, initial point 'E_ini_fun' is needed.")
  }
  #Not from equilibrium : E needed
  if (from.eq==FALSE&length(E_fun)==0) {
    stop("If analysed point is not equilibrium, enter a point 'E_fun'.")
  }
  
  #### Initialization
  
  # rescale E0 and E to have sum(E0) = Etot_fun = sum(E)
  E_ini_fun <- E_ini_fun*Etot_fun/sum(E_ini_fun)
  E_fun <- E_fun*Etot_fun/sum(E_fun)
  
  
  # if analysed point is chosen at equilibrium
  if (from.eq==TRUE) {
    
    #Independency or competition only : from theoretical equilibrium
    if (correl_fun=="SC"|correl_fun=="Comp") {
      e_fun <- predict_th(A_fun,correl_fun)$pred_e
    }
    
    #Positive co-regulation
    if (correl_fun=="RegPos") {
      #do not chose directly E, because of infinite value of E at theoretical equilibrium in this case
      tau_r <- 0.95
      e_fun <- droite_e(tau_r,E_ini_fun,B_fun)
    }
    
    #In other cases, effective equilibrium
    if (correl_fun=="RegNeg"|correl_fun=="CRPos"|correl_fun=="CRNeg") {
      eq_eff <- predict_eff(E_ini_fun,B_fun,A_fun,correl_fun)
      e_fun <- eq_eff$pred_e
    }
    
    #to have sum(E) = Etot
    E_fun <- e_fun*Etot_fun
    
    if (correl_fun=="RegNeg") {
      #because in case of RegNeg, Etot0 != Etot at eq
      E_fun <- eq_eff$pred_E
    }
  }
  
  # Flux at analysed point
  J_fun <- flux(E_fun,A_fun,X_fun)
  
  # Range of enzyme variation
  #= all possible value of "mutated" enzyme
  x <- seq(0,Etot_fun,by=0.01)
  #x <- seq(-10,10,by=0.1) ## Gamme de taille de mutation
  

  
  
  ########### Variables
  # i = "mutated" enzyme
  # y = concentration after mutation of i (value between 0 and Etot, taken in x)
  # j = other enzymes that do not mutate
  # Em = value of Ej after mutation of i
  # vec = vectors contening all Em for every enzymes when i is mutated (vector of length n)
  # Jm = value of J corresponding to each value of y
  # vecJ = vetor contening all values of Jm for each value of y for one enzyme i (same length as x)
  #  /!\ vecJ is reset for each mutated enzyme i
  # J_all = matrix contening all Jm, with mutated enzyme i in column and vecJ corresponding to i in row (one row = one value of x)
  
  
  # Matrix of n col (for each "mutated" enzyme) and ... row for each possible value for "mutated" enzyme
  J_all <- matrix(0,nrow=length(x),ncol=n_fun)
  sel_all <- matrix(0,nrow=length(x),ncol=n_fun)
  sel_all_tru <- matrix(0,nrow=length(x),ncol=n_fun)
  tau_all <- NULL
  if (correl_fun!="SC"&correl_fun!="Comp"){
    #in cases of regulation only
    tau_all <- matrix(0,nrow=length(x),ncol=n_fun)
  }
  
  ############
  
  #for each enzyme. i is the mutated one
  for (i in 1:n_fun) {
    
    # vector for this enzyme
    vecJ <- NULL
    vecsel <- NULL
    vecseltru <- NULL
    vectau <- NULL
    
    # for each mutated value of i : y is concentration of enzyme i after its mutation
    for (y in x) {
      #for (z in x) { # z = mutation size of i
      
      # value of all enzyme after mutation of i
      #take only positive concentrations
      vecE <- NULL
      #conserve negative concentration
      find_tau <- NULL
      
      # computes value of all enzyme concentrations
      for (j in 1:n_fun) {
        
        # if j is the mutated enzyme
        if (j==i) {
          Em <- y
          #Em <- E[j]+z
          
          # if j is not the mutated enzyme
        } else {
          
          # Independency
          if (correl_fun =="SC") {
            Em <- E_fun[j]
          }
          
          #Competition
          if (correl_fun=="Comp") {
            #proportions of non-mutated enzyme are constant before and after mutation
            Em <- (E_fun[j]*(Etot_fun-y)/(Etot_fun-E_fun[i]))
            #Em <- (E[j]*(Etot-E[i]-z)/(Etot-E[i]))
          }
          
          # Regulation only
          if (correl_fun=="RegPos"|correl_fun=="RegNeg") {
            Em <- E_fun[j]+beta_fun[i,j]*(y-E_fun[i])
            #Em <- E[j]+(B[i]/B[j])*(y-E[i])
          }
          
          # Competition and regulation
          if (correl_fun=="CRPos"|correl_fun=="CRNeg") {
            #computes tau position of the mutant
            thor <- (y - E_fun[i])/(Etot_fun/B_fun[i]-E_fun[i])
            #computes Em_j value from tau position
            Em <- E_fun[j] + thor*(Etot_fun/B_fun[j]-E_fun[j])
            #indirectly, it is equal to Ej=Ej^r+alpha_ij * delta_i (where delta_i = y-Ei^r)
            if(with.alpha==TRUE) {
              Em <- E_fun[j] + ((Etot_fun*beta_fun[i,j]/B_fun[i] - E_fun[j])/(Etot_fun/B_fun[i] - E_fun[i]))*(y - E_fun[i])
            }
            
            #   thor <- (y - E[i])/(Etot/B[i]-E[i]) #calcul du tau Ã  partir de l'enzyme i mutante
            #   vec <- droiteE_CR(thor,E,B)
          }
        }
        
        #to avoid negative concentrations, but need Em to find tau (store Em in find_tau)
        find_tau <- c(find_tau,Em)
        if (Em<0){
          Em <- 0
        }
        
        # store value of each enzyme after mutation
        vecE <- c(vecE,Em)
      }
      
      
      ### Computes interest variable
      
      # flux value
      Jm <- flux(vecE,A_fun,X_fun)
      
      # apparent mutation effect
      delta <- y - E_fun[i]
      #delta <- z
      
      # selection coefficient
      sel <- coef_sel.continue(i,vecE,A_fun,delta,correl_fun,beta_fun)
      sel_tru <- (Jm-J_fun)/J_fun
      
      #store to every of mutated enzyme i
      vecJ <- c(vecJ, Jm)
      vecsel <- c(vecsel,sel)
      vecseltru <- c(vecseltru,sel_tru)
      if (correl_fun!="SC"&correl_fun!="Comp"){
        if (grp.reg==FALSE) {
          tau <- droite_tau(find_tau,E_ini_fun,B_fun)
        } else {tau <- 0} #if there is regulation group, the line does not have sense between all enzymes
        vectau <- c(vectau,tau)
      }
      #end of y loop 
    }
    
    # store all data
    J_all[,i] <- vecJ
    sel_all[,i] <- vecsel
    sel_all_tru[,i] <- vecseltru
    if (correl_fun!="SC"&correl_fun!="Comp"){
      tau_all[,i] <- vectau
    }
  }

  return(invisible(list(x=x,J=J_all,sel_disc=sel_all_tru,sel_cont=sel_all,tau=tau_all,
                        param=list(correl=correl_fun,A=A_fun,beta=beta_fun,B=B_fun,E0=E_ini_fun,E=E_fun,Etot=Etot_fun,X=X_fun,n=n_fun))))
}





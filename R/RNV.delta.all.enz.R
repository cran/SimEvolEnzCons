#' Delta at RNV for all enzymes
#'
#' This function computes actual mutation effect \eqn{\delta} at RNV for each enzyme 
#' 
#'
#' @details 
#' The \emph{Range of Neutral Variations} (RNV) are mutant concentration values such as coefficient selection is between \eqn{1/(2N)} and \eqn{-1/(2N)}.
#' 
#' Inferior (resp. superior) bound of RNV corresponds to selection coefficient equal to \eqn{-1/(2N)} (resp. \eqn{1/(2N)}).
#' 
#' Function \code{RNV.delta.all.enz} computes the actual mutation effect \eqn{delta_i} at RNV bounds, where \eqn{i} is the enzyme targeted by the mutation.
#' 
#' 
#' Depending on applied constraints \code{correl_fun}, it exists 1 or 2 RNV.
#' In case of independence (\code{"SC"}) or positive regulation between all enzymes (\code{"RegPos"}), flux has no limit and there is only one RNV.
#' In other cases (competition and/or negative regulation), because flux can reach a maximum, there is two RNV:
#' a "near" one, for small mutations, and a "far" one for big mutations that put mutants in the other side of flux dome.
#' 
#' 
#' \bold{Known bug}
#' 
#' Due to use of \code{\link{range_delta}} to limit search area, output \eqn{\delta} is computed to have mutant concentration between 0 and \code{sum(E_res_fun)} (resident total concentration).
#' If \eqn{\delta} is too high (mutant concentration over resident total concentration), output \eqn{\delta} could be \code{NA} rather than a numeric value. This case might happen when \code{N_fun} is too low and \code{correl_fun="SC"} or \code{"RegPos"} (no limit on total concentration), 
#' 
#' 
#' @importFrom stats optimise
#' @importFrom stats uniroot
#' 
#' 
#' 
#' @usage RNV.delta.all.enz(E_res_fun,A_fun,N_fun,correl_fun,beta_fun=NULL,tol_fun=0.0001)
#'
#' @inheritParams simul.next.resident
#' @param A_fun Numeric vector of activities
#' @inheritParams range_delta
#'  
#'
#'
#' @return List of \emph{n} elements, one for each enzyme \emph{i} considering that \eqn{i} is targeted by the mutation.
#' 
#' For each element of the list, we have a vector of length 2 or 4, depending on applied constraint \code{correl_fun}. Values of this vector are:
#' \enumerate{
#'    \item \emph{delta_i} value for inferior bounds of the near RNV
#'    \item \emph{delta_i} value for superior bounds of the near RNV
#'    \item \emph{delta_i} value for inferior bounds of the far RNV (if exists)
#'    \item \emph{delta_i} value for superior bounds of the far RNV (if exists)
#' }
#' 
#' If superior bound is not accessible, value is \code{NA}.
#' 
#' Note that \code{n} is the number of enzymes, which is the length of \code{E_ini_fun}.
#' 
#' 
#' @seealso 
#' \eqn{delta} at RNV bounds is obtained by nullify the expression in \code{\link{odd.discrete.sel.coef}}.
#' 
#'
#' @references 
#' Coton et al. (2021)
#'
#' @examples
#' Er <- c(30,30,30)
#' A <- c(1,10,30)
#' beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
#' correl <- "CRPos"
#' N <- 1000
#'
#' RNV.delta.all.enz(Er,A,N,correl,beta)
#'
#' correl <- "SC"
#' RNV.delta.all.enz(Er,A,N,correl)
#'
#'
#'
#' @export







RNV.delta.all.enz <- function(E_res_fun,A_fun,N_fun,correl_fun,beta_fun=NULL,tol_fun=0.0001) {
  
  ###### Setting
  n_fun <- length(E_res_fun)
  
  
  
  ####### Tests
  if (length(A_fun)!=n_fun) {
    stop("An enzyme is missing. A_fun and E_res_fun need to have the same length.")
  }
  is.correl.authorized(correl_fun)
  is.beta.accurate(beta_fun,n_fun,correl_fun)
  
  
  
  
  #######" Initialization
  
  #compute limits of Neutral Zone in terms of selection coefficient
  limits_NZ <- c(-1/(2*N_fun),1/(2*N_fun))  #c("inf","sup")
  nb_limits_NZ <- length(limits_NZ)
  
  #compute all redistribution coefficients
  Alpha_fun <- alpha_ij(E_res_fun,correl_fun,beta_fun)
  
    
    
    
  ###### Computation
  
   #All computed delta are putted in a list.
  #Each element of the list corresponds to one enzyme, and in each element, a vector of inf and sup bounds for "small" mutations then for "big" mutations
  #ens_delta_NZ <- vector("list",length=n_fun)
  ens_delta_NZ_bis <- vector("list",length=n_fun) #idem, mais pour la reprÃ©sentation graphique
  
  # for each target enzyme by the mutation
  for (i_fun in 1:n_fun) {
    
    ### definition of exploration zone and its limits for delta depending on constraints
    
    #compute limits of delta such as 0 < Em < Etot <=> -Er/alpha < delta_i < (Etot-Er)/alpha
    lim_delta <- range_delta(E_res_fun,Alpha_fun,i_fun,tol_fun)
    
    # in case of independency or regulation between all enzymes
    if (correl_fun=="SC"|correl_fun=="RegPos") {
      # flux increase indefinitely => only 1 RNV
      nb_RNV_fun <- 1
      # one RNV = no partition of the search area
      lim_explo <- lim_delta
      #set with a random value because not used
      delta_smax <- TRUE
      smax <- TRUE
      
      # in other case
    } else {
      # flux can reach a maximum => 2 RNV ("small" and "big" mutations)
      nb_RNV_fun <- 2
      #search maximum of studied function to partionate search area of delta ; limits_NZ[1], because form of the studied function do not change with the use limits of NZ
      delta_smax <- optimise(f=odd.discrete.sel.coef, interval=lim_delta,maximum=TRUE,i_fun, E_res_fun, A_fun, Alpha_fun,limits_NZ[1], tol=0.00000001)
      #computes corresponding selection coefficient ?
      smax <- coef_sel.discrete(E_res_fun,A_fun,mut.E.indirect(delta_smax$maximum,E_res_fun,Alpha_fun,i_fun),A_fun) #calcul du coef de sÃ©lection correspondant au max de la fonction
      #exploration in two parts: between -Er and max then between max and Etot0-Er
      lim_explo <- c(lim_delta[1],delta_smax$maximum,lim_delta[2])
      
      #if delta corresponding to maximum of the studied function is negative
      # if (delta_smax$maximum<0) {
      #   reverse bounds, because in this case, "small" mutation are between max and Etot-Er, and "big" mutations between -Er and max
      #   lim_explo <- rev(lim_explo)
      # }
      
    }
    
    
    
    
    ### search delta_i such as s_i= -1/2N then +1/2N 
    
    #for each possible RNV ("small" and "big" mutations)
    for (t_fun in 1:nb_RNV_fun) {
      
      # foe each limits of the neutral zone in terms of selection coefficient (inf then sup)
      for (l_fun in 1:nb_limits_NZ) {
        
        #which the limit of NZ in this case
        this_limit_NZ <- limits_NZ[l_fun]
        
        #values of the studied function at the exploration bounds of delta
        value_lim_explo <- sapply(lim_explo[t_fun:(t_fun+1)], odd.discrete.sel.coef, i_fun, E_res_fun, A_fun, Alpha_fun, this_limit_NZ)
          #sapply(lim_explo[t_fun:(t_fun+1)],other.sel.dis,E_fun,A_fun,N_fun,correl_fun,Alpha_fun,i_fun,this_limit_NZ) #valeur du coef de selection Ã  cette borne et pour les limites d'explo choisies
        
        ###compute delta at RNV bounds
        
        #if we searched at the superior bounds of RNV (l_fun == 2) but the maximal selection coefficient < to the superior bounds of the NZ (1/(2N))
        if (smax<=this_limit_NZ&l_fun==2) {
          #this RNV bounds cannot exist for delta
          ens_delta_NZ_bis[[i_fun]] <- c(ens_delta_NZ_bis[[i_fun]],NA)
          #ens_delta_NZ[[i_fun]] <- c(ens_delta_NZ[[i_fun]],tail(ens_delta_NZ[[i_fun]],1)) #on rÃ©plique la borne inf

          # in other cases (superior bounds exits)
        } else {
          
          #if for the studied function, values at the exploration limits have the sign
          if (sign(value_lim_explo[1])==sign(value_lim_explo[2])) { 
            #there is no zero for the studied function between these two limits, so no bounds for the RNV
            ens_delta_NZ_bis[[i_fun]] <- c(ens_delta_NZ_bis[[i_fun]],NA)
            #ens_delta_NZ[[i_fun]] <- c(ens_delta_NZ[[i_fun]],NA)
            
            # in other cases: cases SC&&RegPos, or not at max, or in inf bounds
          } else {
            #all RNV
            lim_ZN_i <-uniroot(f=odd.discrete.sel.coef, interval=lim_explo[t_fun:(t_fun+1)],i_fun,E_res_fun,A_fun,Alpha_fun,this_limit_NZ, tol=0.00000001)
            #ens_delta_NZ[[i_fun]] <- c(ens_delta_NZ[[i_fun]],lim_ZN_i$root)
            ens_delta_NZ_bis[[i_fun]] <- c(ens_delta_NZ_bis[[i_fun]],lim_ZN_i$root) 
          }
        }
        
      }
    }
    
    
    
    ###for each enzyme, shaping results
    
    #only in case of 2 RNV
    if (nb_RNV_fun==2) {
      
      # if (smax<=limits_NZ[2]) { #si la borne sup n'est pas accessible
      #   ens_delta_NZ[[i_fun]] <- c(ens_delta_NZ[[i_fun]][c(1,3,2,4)]) #on inverse les deux du milieu pour dupliquer la mÃªme RNV
      # }
      
      #if delta corresponding to maximum of the studied function is negative
      if (delta_smax$maximum<0) { #dans le cas oÃ¹ le max est nÃ©gatif, les petites mutations sont Ã  droite (entre le max et Etot-E_i) et non Ã  gauche (entre -Ei et le max)

        #reverse bounds, because in this case, "small" mutation are between max and Etot-Er, and "big" mutations between -Er and max
        ens_delta_NZ_bis[[i_fun]] <- c(ens_delta_NZ_bis[[i_fun]][c(3,4,1,2)])
        
        
        #ens_delta_NZ[[i_fun]] <- c(ens_delta_NZ[[i_fun]][3:4],ens_delta_NZ[[i_fun]][1:2]) #on inverse le placement
        #ens_delta_NZ[[i_fun]] <- c(ens_delta_NZ[[i_fun]][c(3,4,1,2)])
        #NB : si le delta_max est nÃ©gatif mais que la borne sup n'est pas accessible, cela ne change pas vu que les deux sont identiques
        
      }
    }
  }
  
  #return list, with for each element of the liste, one enzyme
  #and for eache enzyme: nearest "small" RNV, inf bounds - nearest "small" RNV, sup bounds - farthest "big" RNV, inf bounds - farthest "big" RNV, sup bounds
  # nearest RNV correspond to "small" mutations, and farthest RNV to "big" mutations
  #return(c(ens_delta_NZ,ens_delta_NZ_bis))
  return(ens_delta_NZ_bis)
}



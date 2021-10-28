#' Prediction of equilibrium with regulation groups
#'
#' Gives the equilibrium for intra-group, inter-group and total relative concentrations at equilibrium
#' 
#'
#' @details 
#' Gives values at effective equilibrium for intra-group \eqn{e_i^q}, inter-group \eqn{e^q} and total \eqn{e_i} relative concentrations, and group driving variable \eqn{\tau^q},
#' and also for absolute concentrations \eqn{E^i} and concentrations sum in groups \eqn{E^q}.
#' The equilibrium corresponds to null derivative for relative concentrations.
#' 
#' However, does not compute the theoretical intra-group equilibrium when there is competition, which is \eqn{e_i^q = 1/B_i}.
#' 
#' 
#' @usage predict_grp(E_ini_fun,beta_fun,A_fun,correl_fun, tol=0.00000001)
#'     
#'
#' @inheritParams predict_eff
#' @inheritParams class_group
#' 
#'
#'
#' @return List of seven elements:
#' \itemize{
#' \item \code{$pred_eiq}: numeric vector of intra-group relative concentrations \eqn{e_i^q} at equilibrium. Same length as \code{A_fun}.
#' \item \code{$pred_eq}: numeric vector of inter-group relative concentrations \eqn{e^q} at equilibrium. The length is the number of regulation groups.
#' \item \code{$pred_ei}: numeric vector of total relative concentrations \eqn{e_i} at equilibrium. Same length as \code{A_fun}.
#' \item \code{$pred_tau}: numeric vector of driving variable \eqn{tau^q} at equilibrium. The length is the number of regulation groups.
#' \item \code{$pred_Ei}: numeric vector of enzyme absolute concentrations \eqn{E_i} at equilibrium. Same length as \code{A_fun}.
#' \item \code{$pred_Eq}: numeric vector of sum of absolute concentrations in groups \eqn{E^q} at equilibrium. The length is the number of regulation groups.
#' \item \code{$pred_Etot}: numeric value of total concentration at equilibrium.
#' }
#' 
#' @section Special results:
#' 
#' When there are more than one positive or negative group and singletons with competition (\code{"CRPos"} or \code{"CRNeg"}), the equilibria are not predictable.
#' 
#'
#' 
#' @seealso 
#' Use function \code{\link{activities}} to compute enzyme activities.
#' 
#' Use function \code{\link{is.correl.authorized}} to see allowed constraints for \code{correl_fun}.
#' 
#' Use function \code{\link{predict_th}} (resp. \code{\link{predict_eff}}) to compute theoretical (resp. effective) equilibrium when there is no regulation groups (enzymes are all independent or all co-regulated).
#'
#' @examples
#' #### For independancy "SC"
#' A <- c(1,10,30)
#' E0 <- c(30,30,30)
#' beta <- diag(1,3)
#' 
#' eq <- predict_grp(E0,beta,A,"SC")
#' #same results for pred_e and pred_ei
#' eq_th <- predict_th(A,"SC")
#' 
#' ###### In presence of regulation, all enzyme co-regulated
#' A <- c(1,10,30)
#' beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
#' B <- apply(beta,1,sumbis) 
#' 
#' eq_grp <- predict_grp(E0,beta,A,"CRPos")
#' #same results for pred_e and pred_ei
#' eq_eff <- predict_eff(E0,B,A,"CRPos")
#' 
#' 
#' #Two groups: one negative group + one singleton
#' n <- 3
#' beta <- diag(1,n) 
#' beta[1,2] <- -0.32 
#' beta[2,1] <- 1/beta[1,2]
#' 
#' eq_grp <- predict_grp(E0,beta,A,"RegNeg")
#' eq_grp <- predict_grp(E0,beta,A,"CRNeg")
#' 
#' 
#' #Two groups: one positive group + one singleton 
#' n <- 3
#' beta <- diag(1,n) 
#' beta[1,2] <- 0.43 
#' beta[2,1] <- 1/beta[1,2]
#' 
#' eq_grp <- predict_grp(E0,beta,A,"RegPos")
#' eq_grp <- predict_grp(E0,beta,A,"CRPos")
#' 
#' 
#' \donttest{
#' #With saved simulation
#' data(data_sim_RegPos)
#' n <- data_sim_RegPos$param$n
#' num_s <- 1
#' pred_eq <- predict_grp(data_sim_RegPos$list_init$E0[num_s,1:n],
#' data_sim_RegPos$param$beta,data_sim_RegPos$list_init$A0[num_s,1:n],data_sim_RegPos$param$correl)
#' 
#' data(data_sim_RegNeg_1grpNeg1grpPos)
#' pred_eq <- predict_grp(data_sim_RegNeg_1grpNeg1grpPos$list_init$E0[num_s,],
#' data_sim_RegNeg_1grpNeg1grpPos$param$beta,c(1,10,30,50),"RegNeg")
#' 
#' }
#'
#' @export



predict_grp <- function(E_ini_fun,beta_fun,A_fun,correl_fun, tol=0.00000001) {
  #Number of enzymes
  n_fun <- length(A_fun)
  #To avoid problem due to use of dataframe
  E_ini_fun <- as.matrix(E_ini_fun)
  #beta without reg
  if (correl_fun=="Comp"|correl_fun=="SC") {
    beta_fun <- diag(1,n_fun)
  }
  #B
  B_fun <- compute.B.from.beta(beta_fun)
  #list of groups
  L_Phi_fun <- class_group(beta_fun)
  p_fun <- length(L_Phi_fun)
  #types of groups
  grp_typ <- group_types(beta_fun)
  
  
  # initial sum of concentration in groups : Eq0
  # Eini_btw <- rep(NA,p_fun)
  # for (q in 1:p_fun) {
  #   phi_q <- unlist(L_Phi_fun[[q]])
  #   Eini_btw[q] <- sum(E_ini_fun[phi_q])
  # }
  
  
  
  ###### Verif parameters
  # cases allowed; if not, program stops
  is.correl.authorized(correl_fun)
  # adequacy of beta_fun
  is.beta.accurate(beta_fun,n_fun,correl_fun)
  # missing parameters
  if (length(E_ini_fun)!=n_fun) {
    stop("An enzyme are missing. A_fun and E_ini_fun have not the same length.")
  }
  
  
  
  
  #set intra, inter and total relative concentrations + driving variable of groups
  eq_ingrp <- rep(NA,n_fun) #eiq
  eq_btwgrp <- rep(NA,p_fun) #eq
  eq_total <- rep(NA,n_fun) #ei
  eq_tau <- rep(NA,p_fun) #tauq
  Eq_btwgrp <- rep(NA,p_fun) #Eq
  E_fun <- rep(NA,n_fun) #Ei
  
  
  ######## Without competition
  
  if (any(correl_fun==c("SC","RegPos","RegNeg"))) {
    
    ## Within groups : eiq
    
    for (q in 1:p_fun) {
      phi_q <- L_Phi_fun[[q]]
      
      #which kind of group
      which.typ.q <- search_group(q,grp_typ)
      
      #if positive group
      if(names(which.typ.q)=="grp_pos") {
        # eq th = 1/B
        eq_ingrp[phi_q] <- predict_th(A_fun[phi_q],"RegPos",B_fun[phi_q])$pred_e
        eq_tau[q] <- 1
        E_fun[phi_q] <- Inf
      }
      
      #if singleton
      if(names(which.typ.q)=="grp_single") {
        # eiq = 1
        eq_ingrp[phi_q] <- 1
        #tau_q = NA
        E_fun[phi_q] <- Inf
      }
      
      #if negative group
      if(names(which.typ.q)=="grp_neg") {
        # eq eff
        eq_values <- predict_eff(E_ini_fun[phi_q],B_fun[phi_q],A_fun[phi_q],"RegNeg")
        # eiq
        eq_ingrp[phi_q] <- eq_values$pred_e
        # tau
        eq_tau[q] <- eq_values$pred_tau
        # Ei
        #because relationship between tau_q and Ei is the same as droite_E.Reg but with Eq0 (which is the input in predict_eff) rather than Etot0 
        E_fun[phi_q] <- eq_values$pred_E
        # Eq : sum of concentrations of the group
        #Eq_btwgrp[q] <- sum(eq_values$pred_E) 
      }
      
      # Eq : sum of concentrations of the group
      Eq_btwgrp[q] <- sum(E_fun[phi_q])
      
    }
    #end of within groups
    
    
    ## Between groups :  eq
    #apparent activities of groups
    A_app <- unlist(lapply(L_Phi_fun,apparent.activities.Aq,A_fun,B_fun,correl_fun))
    
    
    #only positive groups and singletons  
    if (length(grp_typ$grp_neg)==0) {
      #eq th
      eq_btwgrp <- predict_th(A_app,"SC")$pred_e
    }
    
    #only negative groups
    if (length(grp_typ$grp_neg)==p_fun) {
      #fixed concentrations at eq eff
      eq_btwgrp <- Eq_btwgrp/sum(Eq_btwgrp)
    }
    
    #positive +/- singletons AND negative groups
    if (length(grp_typ$grp_neg)>0 & length(grp_typ$grp_neg)<p_fun) {
      #for each group
      for (q in 1:p_fun) {
        #which kind of group
        which.typ.q <- search_group(q,grp_typ)
        #if negative group
        if(names(which.typ.q)=="grp_neg") {
          #because of positive groups, Etot -> Inf
          eq_btwgrp[q] <- 0
        } else {
          #recup numbers of positive and negative groups
          grp_pos_sgl <- c(grp_typ$grp_pos,grp_typ$grp_single)
          #eq th with apparent activities only for these groups
          eq_btwgrp[grp_pos_sgl] <- predict_th(A_app[grp_pos_sgl],"SC")$pred_e
        }
      }
    }
    #end of between groups
    
    
    #total concentration = sum of concentration
    Etot_fun <- sum(E_fun)
    
    
    # ## only positive groups and singletons  
    # if (length(grp_typ$grp_pos)==0) {
    #   
    #   # Within group : eiq*
    #   if (correl_fun=="SC") {
    #     #1/B = 1 for all enzymes
    #     eq_ingrp <- rep(1,n_fun)
    #   } else {
    #     #1/B
    #     eq_ingrp <- predict_th(A_fun,correl_fun,B_fun)$pred_e
    #     eq_tau <- rep(1,p_fun) ######## /!\ NA si singleton /!\
    #   }
    #   
    #   
    #   # Between groups : eq*
    #   #apparent activities of groups
    #   A_app <- unlist(lapply(L_Phi_fun,apparent.activities.Aq,A_fun,B_fun,correl_fun))
    #   
    #   eq_btwgrp <- predict_th(A_app,"SC")$pred_e
    #   
    # }
    # 
    # ## only negative groups
    # if (length(grp_typ$grp_neg)==p_fun) {
    #   # Within group : eiq~
    #   
    #   #for each group
    #   for (q in 1:p) {
    #     #enzyme in group q
    #     phi_q <- L_Phi_fun[[q]]
    #     #compute effective equilibrium with enzymes in group q
    #     eq_values <- predict_eff(E_ini_fun[phi_q],B_fun[phi_q],A_fun[phi_q],correl_fun)
    #     # eiq
    #     eq_ingrp[phi_q] <- eq_values$pred_e
    #     # tau
    #     eq_tau[q] <- eq_values$pred_tau
    #     # Eq : sum of concentrations of the group
    #     Eq_btwgrp[q] <- sum(eq_values$pred_E) 
    #     
    #   }
    #   
    #   # Between groups : eq
    #   eq_btwgrp <- Eq_btwgrp/sum(Eq_btwgrp)
    #   
    #   
    # }
    # 
    # ## positive and negative groups
    # if (length(grp_typ$grp_neg)>0 & length(grp_typ$grp_neg)<p_fun) {
    #   
    # }
    
    
    
  
  }
  #end of without competition
  
  
  
  
  
  ####### With competition
  
  if (any(correl_fun==c("Comp","CRPos","CRNeg"))) {
    
    
    ## Within group : eiq
    for (q in 1:p_fun) {
      #take the enzyme numbers of this group q
      phi_q <- unlist(L_Phi_fun[[q]])
      
      if (correl_fun=="Comp"|length(phi_q)==1) {
        #if enzymes are all independent or if there only one enzyme in the group, always Ei/Eq=1
        eq_ingrp[phi_q] <- 1
      } else {
        #how many groups are not singleton
        nb_not_single <- sum(length(grp_typ$grp_pos),length(grp_typ$grp_neg))
        
        #if there only one positive or negative group and any number of singletons
        if (nb_not_single==1) {
          #compute effective equilibrium with enzymes in group q
          eq_values <- predict_eff(E_ini_fun[phi_q],B_fun[phi_q],A_fun[phi_q],correl_fun)
          # eiq
          eq_ingrp[phi_q] <- eq_values$pred_e
          # tau
          eq_tau <- eq_values$pred_tau
        }
        #else we cannot compute the equilibrium
        
      }
    
    #end of loop on q
    }
    
    
    ## Between groups : eq
    #apparent activities of groups
    A_app <- unlist(lapply(L_Phi_fun,apparent.activities.Aq,A_fun,B_fun,correl_fun,eiq_eff=eq_ingrp))
    
    #same as theoretical eq without regulation with apparent activities
    eq_btwgrp <- predict_th(A_app,"Comp")$pred_e
    
  }
  #end of with competition
  
  
  
  ##### With or without competition
  
  ## Total : eiq
  
  for (j in 1:n_fun){#pour toutes les enzymes
    #search in which group is the enzyme j with search_group, then ei = eiq*eq
    eq_total[j] <- eq_ingrp[j]*eq_btwgrp[search_group(j,L_Phi_fun)]
  }
  
  
  ## absolute concentrations in case of competition
  if (any(correl_fun==c("Comp","CRPos","CRNeg"))) {
    #fixed total concentration
    Etot_fun <- sum(E_ini_fun)
    # Eq = eq*Etot
    Eq_btwgrp <- eq_btwgrp*Etot_fun
    # Ei = ei*Etot
    E_fun <- eq_total*Etot_fun
  }
  
  
  ##### Output
  
  return(list(pred_eiq=eq_ingrp, pred_eq=eq_btwgrp, pred_ei=eq_total, pred_tau=eq_tau, pred_Ei=E_fun, pred_Eq=Eq_btwgrp, pred_Etot=Etot_fun))
}


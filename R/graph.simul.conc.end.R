#' Graphics of simulation ends
#' 
#' Gives graphics of enzyme concentrations two-by-two at end of simulation
#' 
#' @details 
#' This function shows the concentration of one enzyme against the concentration of another enzyme, for selected enzymes.
#' 
#' Simulation ends are supposed to be near to equilibrium, given the population size and therefore the neutral zone.
#' 
#' 
#' 
#' @usage graph.simul.conc.end(all_res_sim,new.window=FALSE,add.eq=TRUE,
#' which.sim=NULL,which.enz=NULL,...)
#' 
#' @inheritParams graph.simul.by.time.by.sim
#' @param which.enz Numeric vector containing integer numbers between 1 and \code{n}. Which enzymes would you represent? If \code{NULL} (default), all enzymes would be represented.
#'
#'
#' @import graphics
#' 
#' @seealso 
#' Use function \code{\link{simul.evol.enz.multiple}} to simulate enzyme evolution.
#'
#' Use function \code{\link{graph.simul.by.time.by.sim}} to see enzyme concentrations through time.
#' 
#' 
#' @examples
#'
#' data(data_sim_CRNeg_1grpNeg1sgl)
#' graph.simul.conc.end(data_sim_CRNeg_1grpNeg1sgl)
#'
#'
#' @export


graph.simul.conc.end <- function(all_res_sim,new.window=FALSE,add.eq=TRUE,
                                 which.sim=NULL,which.enz=NULL,...) {
  
  #input parameters
  n <- all_res_sim$param$n #number of enzymes
  nsim <- all_res_sim$param$nsim #number of simulations
  Keq_fun <- all_res_sim$param$Keq #equilibrium constants
  X_fun <- all_res_sim$param$X #numerator of flux
  beta_fun <- all_res_sim$param$beta #co-regulation coefficients
  B_fun <- all_res_sim$param$B #global co-regulation coefficients
  correl_fun <- all_res_sim$param$correl #constraints on system
  same.kin0 <- all_res_sim$param$same.kin0
  pmutA <- all_res_sim$param$pmutA #proba mutation of kin
  
  #simulation results
  tabR <- all_res_sim$tabR
  
  #begin of simulations
  E_ini <- all_res_sim$list_init$E0
  
  #end of simulations
  E_f <- all_res_sim$list_final$E_f
  A_f <- all_res_sim$list_final$A_f
  
  #which simulations will be represented. If NULL, all
  if (length(which.sim)==0) {which.sim <- 1:nsim}
  #same for enzymes
  if (length(which.enz)==0) {which.enz <- 1:n}
  
  
  ##### Equilibrium #####
  
  ### list of equilibrium for all initial concentration and activities
  all_eq_grp <- vector("list",length=nsim)
  eq_E <- NULL
  for (i in 1:nsim) {
    all_eq_grp[[i]] <- predict_grp(E_ini[i,1:n],beta_fun,A_f[i,1:n],correl_fun)
    #take equilibrium for absolute concentrations only
    eq_E <- rbind(eq_E,all_eq_grp[[i]]$pred_Ei)
  }
  
  #computes equilibrium for competition alone
  eq_comp_E <- NULL
  #if there is competition
  if (any(correl_fun==c("CRPos","CRNeg"))) {
    #if no variation of A between and within simulations
    if (pmutA==0&same.kin0==TRUE) {
      #if Etot0 is identical between simulations (compared with first simulation)
      if (all.equal(as.numeric(apply(E_ini,1,sum)), rep(sum(E_ini[1,]),nsim))) {
        Etot0 <- sum(E_ini[1,])
        eq_comp_E <- predict_th(A_f[1,],"Comp")$pred_e*Etot0
      }
    }
  }
  
  
  ##### Graphics #######
  
  #ordered Ej
  #Ej (row 1) in x-axis and Ek (row 2) in y-axis. One column by figure.
  
  ordE <- NULL #matrix(0,nrow=2,ncol=factorial(n-1))
  #for each possible case 2-by-2
  #n-1: last enzyme not in x-axis
  # for (j in 1:(n-1)) {
  #   #triangle sup
  #   for (k in (j+1):n) {
  #     ordE <- cbind(ordE,c(j,k))
  #   }
  # }
  #number of selected enzymes
  nb.sel.enz <- length(which.enz)
  #last enzyme not in x-axis
  for (j in which.enz[1:(nb.sel.enz-1)]) {
    #superior triangle: avoid to repeat graph -> already represented enzymes not pass in y-axis
    #which(which.enz==j): because loop on enzyme number, but selection on position in which.enz, so we take position of enzyme j in which.enz
    for (k in which.enz[(which(which.enz==j)+1):nb.sel.enz]) {
      ordE <- cbind(ordE,c(j,k))
    }
  }
  
  #for each enz vs enz
  for (i in 1:ncol(ordE)) {
    
    #takes enzyme j in x-axis and enzyme k in y-axis
    j <- ordE[1,i] ; k <- ordE[2,i]
    
    #Ek fct Ej
    plot(E_f[,j],E_f[,k],xlab=paste0("E",j),ylab=paste0("E",k),pch=20,...)#,cex.main=2,cex.lab=2,cex.axis=1.6)
    
    # equilibrium in red
    if (pmutA==0&add.eq==TRUE) {
      points(eq_E[,j],eq_E[,k],col=2,pch=20,...)
    }
    
    #le max du max
    if (add.eq==TRUE) {
      points(eq_comp_E[j],eq_comp_E[k],bg='green',pch=17,col='green3',...)#bg='violet',pch=21, cex=1.5
    }
    
  }
  
  
}


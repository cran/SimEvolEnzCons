#' @name simul.evol.graph.methods
#' @aliases graph.simul.by.time.by.sim
#' @aliases graph.simul.by.time.by.enz
#' @aliases graph.simul.others.by.sim
#' @aliases graph.simul.by.time.RNV
#' 
#' @title Graphic methods for simulations of enzyme evolution
#'
#' @description Graphics illustrating enzyme evolution simulations obtained by function \code{\link{simul.evol.enz.multiple}}.
#' 
#' Function \code{graph.simul.by.time.by.sim} gives graphics depending on time, colored by simulations.
#' 
#' Function \code{graph.simul.by.time.by.enz} gives graphics depending on time, colored by enzymes, with series of graphics for each simulation.
#' 
#' Function \code{graph.simul.others.by.sim} gives different graphics depending on other variables than time in x-axis.
#' 
#' Function \code{graph.simul.by.time.RNV} gives graphics of \emph{Range of Neutral Variation} (RNV) for each enzyme.
#' 
#'
#' @details 
#' \emph{If only one simulation may be represented, use preferably function \code{graph.simul.by.time.by.enz}.}
#' 
#' Colors for simulations are taken in palette \code{rainbow}.
#' Colors for enzymes correspond to their number plus one.
#' 
#' \bold{Function \code{graph.simul.by.time.by.sim}} gives graphs of flux, relative concentrations, absolute concentrations, total concentration, kinetic parameters and activities through time.
#' In addition, if exists, gives also driving variable \eqn{\tau} in relation to time.
#' \emph{Lines are colored according to simulation.} There is one graph by enzyme if necessary.
#' Dashed lines correspond to theoretical equilibrium, and dotted lines to effective equilibrium.
#' Every graph follow same scheme:
#' \enumerate{
#'    \item empty graph with time in x-axis and interesting variable in y-axis
#'    \item for each simulation \emph{i}
#'    \item add connected points for current variable for simulation \emph{i}
#'    \item add text for simulation number \emph{i} at end of x-axis
#'    \item eventually, add predicted values
#' }
#' 
#' 
#' \bold{Function \code{graph.simul.by.time.by.enz}} gives graphs of flux, relative concentrations, absolute concentrations, total concentration, kinetic parameters, activities and response coefficients through time. 
#' In addition, if exists, gives also driving variable \eqn{\tau} in relation to time and flux in relation to \eqn{\tau}.
#' \emph{Lines are colored according to enzyme.} There is one graph by simulation. An heading with parameters of current simulation can be added with \code{gr.sim.heading}
#' Dashed lines correspond to theoretical equilibrium, and dotted lines to effective equilibrium.
#' \enumerate{
#'    \item for each simulation \emph{i}
#'    \item line graph with time in x-axis and interesting variable in y-axis
#'    \item eventually, add predicted values
#'    \item add legend
#' }
#' 
#' 
#' \bold{Function \code{graph.simul.others.by.sim}} gives graphs of:
#' \itemize{
#'    \item final concentrations in relation to initial concentrations
#'    \item final relative concentrations in relation to initial relative concentrations
#'    \item final kinetic parameters in relation to initial kinetic parameters
#'    \item final activities in relation to initial activities
#'    \item final concentrations in relation to final activities
#'    \item final relative concentrations in relation to final activities
#'    \item flux in relation to concentrations and activities (3D-graph)
#'    \item flux in relation to absolute or relative concentrations (one graph by enzyme, colored by simulation)
#' }
#' One color by enzyme. The colored numbers correspond to the simulations.
#' 
#' 
#' \bold{Function \code{graph.simul.by.time.RNV}} gives graphs of, for each simulations:
#' \itemize{
#'    \item concentrations with RNV bounds
#'    \item apparent mutation effects \eqn{\delta} at RNV bounds, for each enzyme considering at mutant
#'    \item RNV size
#'    \item RNV size divided by total concentration
#'    \item flux with neutral zone bounds, in relation to time and in relation to enzyme concentrations of each enzyme (where data are ordered to facilitate view)
#' }
#' Lines are colored by enzymes. Bounds of RNV is colored depending on \code{col_RNV}.
#' 
#' 
#' \bold{Graphical parameters}
#' 
#' To modify line width, input both \code{lwd} and \code{lwd.eq}.
#' Input \code{lwd} without input \code{lwd.eq} modifies only equilibrium line width, and not all line width.
#' 
#' 
#' 
#' @param all_res_sim List, the output of function \code{\link{simul.evol.enz.multiple}} (results of evolution simulation).
#' @param new.window Logical. Do graphics appear in a new window? 
#' @param add.eq Logical. Do equilibrium appear on graph? 
#' @param which.sim Numeric vector containing integer numbers between 1 and \code{nsim}. Which simulations would you represent? If \code{NULL} (default), all simulations would be represented.
#' @param ... Arguments to be passed in \code{plot} function, such as \code{lwd} or \code{cex}.
#' @param lwd.eq Numeric. Line width for equilibrium line only.
#'
#' @import graphics
#' @importFrom grDevices rainbow
#' @importFrom scatterplot3d scatterplot3d
#'
#'  
#' @seealso 
#' Use function \code{\link{simul.evol.enz.multiple}} to simulate enzyme evolution.
#' 
#' Function \code{\link[scatterplot3d]{scatterplot3d}} is used to make the 3D-graph in function \code{graph.simul.others.by.sim}.
#' 
#' 
#'
#'
#'
#' @examples 
#'  
#'   # With saved simulation
#' data(data_sim_RegNeg)
#' 
#' graph.simul.by.time.by.sim(data_sim_RegNeg,new.window=TRUE)
#' graph.simul.by.time.by.enz(data_sim_RegNeg,new.window=TRUE,which.sim=c(1))
#' graph.simul.others.by.sim(data_sim_RegNeg,new.window=TRUE)
#' graph.simul.by.time.RNV(data_sim_RegNeg,new.window=TRUE,which.sim=c(1))
#'  
#'  
#'  
#'  \donttest{
#'  #New simulation
#' # case for 3 enzymes
#' n <- 3
#' E0 <- c(30,30,30)
#' kin <- c(1,10,30)
#' Keq <- c(1,1,1)
#' nsim <- 2 # 2 simulations
#' N <- 1000
#' beta <- diag(1,n)
#' beta[upper.tri(beta)] <- c(0.32,0.32*(-0.43),-0.43)
#' #put : beta_12 = 0.32, beta_13 = beta_12 x beta_23, beta_23 = -0.43
#' t_beta <- t(beta) #because R fills matrix column by column
#' beta[lower.tri(beta)] <- 1/t_beta[lower.tri(t_beta)] #beta_ji = 1/beta_ij
#' if (n==3) {beta[lower.tri(beta)] <- 1/beta[upper.tri(beta)]} #only available if n=3
#' correl <- "RegNeg"
#' 
#' evol_sim <- simul.evol.enz.multiple(E0,kin,Keq,nsim,N,correl,beta,npt=250)
#' graph.simul.by.time.by.sim(evol_sim,new.window=TRUE)
#' graph.simul.by.time.by.enz(evol_sim,new.window=TRUE,which.sim=c(1))
#' graph.simul.others.by.sim(evol_sim,new.window=TRUE)
#' graph.simul.by.time.RNV(evol_sim,new.window=TRUE,which.sim=c(1))
#' 
#' }
#' 
#'
#'
#'
#'
#' 
NULL



#Through time, colored by simulation
#################################
#################################
#################################
#' @rdname simul.evol.graph.methods
#' @usage graph.simul.by.time.by.sim(all_res_sim,new.window=FALSE,add.eq=TRUE,which.sim=NULL,
#' gr.J.time=FALSE,gr.e.time=TRUE,gr.E.time=FALSE,gr.Etot.time=FALSE,
#' gr.kin.time=FALSE,gr.A.time=FALSE,gr.tau.time=FALSE,
#' lwd.eq=1.5,...)
#' 
#' @param gr.J.time,gr.e.time,gr.E.time,gr.Etot.time,gr.kin.time,gr.A.time Logical.
#'  Add graph flux \code{J} / relative concentrations \code{e} / absolute concentrations \code{E} / total concentration \code{Etot} / 
#'  kinetic parameters \code{kin} / activities \code{A} in relation to time?
#'  
#'  
#' 
#' @return Function \code{graph.simul.by.time.by.sim} returns invisible list of 5 elements:
#' \itemize{
#'    \item \code{$eq_th_e}: Numeric matrix of \code{n} columns and \code{nsim} rows. Every row corresponds to relative concentrations at theoretical equilibrium computed from initial values of current simulation;
#'    \item \code{$eq_th_r}: Same structure, for response coefficients at theoretical equilibrium;
#'    \item \code{$eq_eff_e}: Same structure, for relative concentrations at effective equilibrium (if exists), else \code{NA};
#'    \item \code{$eq_eff_E}: Same structure, for absolute concentrations at effective equilibrium (if exists), else \code{NA};
#'    \item \code{$eq_eff_tau}: Numeric matrix of one \code{column} and \code{nsim} rows, corresponding to driving variable \eqn{\tau} at effective equilibrium in case of regulation, else \code{NULL}.
#' }
#' 
#' @export

graph.simul.by.time.by.sim <- function(all_res_sim,new.window=FALSE,add.eq=TRUE,which.sim=NULL,
                                          gr.J.time=FALSE,gr.e.time=TRUE,gr.E.time=FALSE,gr.Etot.time=FALSE,gr.kin.time=FALSE,gr.A.time=FALSE,gr.tau.time=FALSE,
                                          lwd.eq=1.5,...) {
  #here, variable 'i' is always used to indicate simulation number and 'j' for enzyme number
  
  ####### Take interessing parameters
  
  #input parameters
  n <- all_res_sim$param$n #number of enzymes
  nsim <- all_res_sim$param$nsim #number of simulations
  Keq_fun <- all_res_sim$param$Keq #equilibrium constants
  X_fun <- all_res_sim$param$X #numerator of flux
  beta_fun <- all_res_sim$param$beta #co-regulation coefficients
  B_fun <- all_res_sim$param$B #global co-regulation coefficients
  correl_fun <- all_res_sim$param$correl #constraints on system
  N_fun <- all_res_sim$param$N #population size
  pasobs <- all_res_sim$param$pasobs
  npt <- all_res_sim$param$npt #number of observations
  ngene <- npt*pasobs #number of generations
  pmutA <- all_res_sim$param$pmutA #proba mutation of kin
  same.E0 <- all_res_sim$param$same.E0
  is.random.E0 <- all_res_sim$param$is.random.E0
  same.kin0 <- all_res_sim$param$same.kin0
  is.random.kin0 <- all_res_sim$param$is.random.kin0
  
  
  #simulation results
  tabR <- all_res_sim$tabR
  tabP_e <- all_res_sim$tabP_e
  tabP_c <- all_res_sim$tabP_r
  
  #initial and final values of E, kin, A and relative concentrations e
  #matrix such as: each simulation in row, initial values from column 1 to n and final values from column n+1 to 2n
  leg_E <- cbind(all_res_sim$list_init$E0,all_res_sim$list_final$E_f)
  leg_kin <- cbind(all_res_sim$list_init$kin0,all_res_sim$list_final$kin_f)
  leg_A <- cbind(all_res_sim$list_init$A0,all_res_sim$list_final$A_f)
  leg_e <- leg_E
  for (i in 1:nsim){
    #initial relative concentrations e0
    leg_e[i,1:n] <- leg_E[i,1:n]/sum(leg_E[i,1:n])
    #final relative concentrations e_f
    leg_e[i,((n+1):(2*n))] <- leg_E[i,((n+1):(2*n))]/sum(leg_E[i,((n+1):(2*n))])
  }
  
  #theoretical equilibrium
  # A_base <- activities(kin_base,Keq)
  # prediction <- predict(A_base,typmut_E,correl)
  all_eq_th <- apply(all_res_sim$list_init$A0,1,predict_th,correl_fun,B_fun)#due to R transformation of matrix in numeric class when there is only one row (nsim=1)...
  #all_eq_th <- apply(leg_A[,1:n],1,predict_th,correl_fun,B_fun)
  eq_th_e <- NULL
  eq_th_r <- NULL
  for (i in 1:nsim) {
    eq_th_e <- rbind(eq_th_e,all_eq_th[[i]]$pred_e)
    eq_th_r <- rbind(eq_th_r,all_eq_th[[i]]$pred_r)
  }
  # enzpred <- eq_th$pred_e
  # cpred <- eq_th$pred_c
  lty_eq <- 2 #dashed lines
  
  #effective equilibrium
  if (correl_fun=="RegNeg"|correl_fun=="CRPos"|correl_fun=="CRNeg") {
    lty_eq <- 3 #dotted lines
    eq_eff_e <- NULL
    eq_eff_E <- NULL
    eq_eff_tau <- NULL
    for (i in 1:nsim) {
      eq_eff_by_sim <- predict_eff(leg_E[i,1:n],B_fun,leg_A[i,1:n],correl_fun)
      eq_eff_e <- rbind(eq_eff_e,as.vector(eq_eff_by_sim$pred_e))
      eq_eff_E <- rbind(eq_eff_E,as.vector(eq_eff_by_sim$pred_E))
      eq_eff_tau <- rbind(eq_eff_tau,as.vector(eq_eff_by_sim$pred_tau))
    }
  } else {
    #effective equilibrium does not exist
    eq_eff_e <- matrix(NA,nrow=nsim,ncol=(2*n))
    eq_eff_E <- matrix(NA,nrow=nsim,ncol=(2*n))
    eq_eff_tau <- matrix(NA,nrow=nsim,ncol=1)
  }
  
  
  
  
  #########################
  #######""" Graphics
  #############################
  
  #which simulations will be represented. If NULL, all
  if (length(which.sim)==0) {which.sim <- 1:nsim}
  
  #color vector for simulations
  mycol_sim <- rainbow(nsim)
  #color for theoretical equilibrium (identical for all simulation if A does not change)
  #if (pmutA==0) {col_th <- 1} else {col_th <- mycol_sim}
  #color for effective equilibrium: its changes if different A0 or E0
  if (same.kin0==TRUE&same.E0==TRUE&pmutA==0) {col_eff <- 1} else {col_eff <- mycol_sim}
  
  #time axis: number of generations, but only 'npt' observations every 'pasobs' generations
  time.axis <- seq(1,ngene,by=pasobs)
    
  
  #####"
  #Every graph follow same scheme:
  #1. empty graph with time in x axis and interesting variable in y axis
  #2. for each simulation
  #3. add connected points for simulation i for current variable
  #4. add text for simulation number at end of x axis
  #5. eventually, add predicted values
  
  
  
  ###### Flux
  if (gr.J.time==TRUE) {
    if (new.window==TRUE) {dev.new()}
    #par(mar=c(5,5,5,5))
    plot(seq(1,ngene,by=pasobs),tabR[1:npt,(2*n+3)],xlim=c(0,ngene+1000),ylim=c(0,max(tabR[,(2*n+3)])),type="n",xlab="Generation",ylab="Flux",...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
    #empty graph: x axis = nb of generations, y axis = flux, ylim = [0,max value of J]
    
    for (i in which.sim){ #for each simulations
      #add connected points on graph for flux for simulation i
      points(seq(1,ngene,by=pasobs),tabR[tabR$sim==i,(2*n+3)],type="l",col=mycol_sim[i],...)#,lwd=1.5)
      #add simulation number at end of x axis
      text(x=(ngene+1000),y=tabR[i*npt,(2*n+3)],labels=i,col=mycol_sim[i])
    }
  }
  
  
  
  
  ##### Relative concentrations
  if (gr.e.time==TRUE) {
    for (j in 1:n){ #for each enzyme
      
      #one graph by enzyme
      if (new.window==TRUE) {dev.new()}
      
      #par(mar=c(5,5,5,5))
      plot(seq(1,ngene,by=pasobs),tabR[1:npt,j],ylim=c(0,1), xlim=c(0,ngene+1000),type="n",xlab="Generation",ylab=paste("Relative enzyme concentration of enzyme ",j, sep=""),...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
      #empty graph: x axis = nb of generations, y axis = current enzyme j, ylim = [0,1] because it's relative concentrations
      
      for (i in which.sim){ #for each simulation
        #add connected points on graph for relative concentrations of enzyme j for simulation i
        points(seq(1,ngene,by=pasobs),tabR[tabR$sim==i,j]/tabR[tabR$sim==i,(2*n+1)],type="l",col=mycol_sim[i],...)#,lwd=1.5)
        #add simulation number at end of x axis
        text(x=(ngene+1000),y=tabR[i*npt,j]/tabR[i*npt,(2*n+1)],labels=i,col=mycol_sim[i])
        
        if (pmutA!=0&add.eq==TRUE) {
          #add prediction of relative concentration during simulation, depending on kin mutations
          points(seq(1,ngene,by=pasobs),y=tabP_e[(tabP_e[,n+1]==i),j],type="l",lty=lty_eq,col=mycol_sim[i],lwd=lwd.eq)#,lwd=1.5)
        }
      }
      
      #predict value (theoretical & effective equilibrium)
      if (pmutA==0&add.eq==TRUE) { 
        abline(h=eq_th_e[,j],lty=2,col=1,lwd=lwd.eq)
        abline(h=eq_eff_e[,j],lty=3,col=col_eff,lwd=lwd.eq)
      }
    }
  }
  
   
   
  ##### Absolute concentrations
  if (gr.E.time==TRUE) {
    for (j in 1:n){
      if (new.window==TRUE) {dev.new()}
      #par(mar=c(5,5,5,5))
      plot(seq(1,ngene,by=pasobs),tabR[1:npt,j],ylim=c(0,max(tabR[,1:n])), xlim=c(0,ngene+1000),type="n",xlab="Generation",ylab=paste("Enzyme concentration of enzyme ",j, sep=""),...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
      
      for (i in which.sim){
        points(seq(1,ngene,by=pasobs),tabR[tabR$sim==i,j],type="l",col=mycol_sim[i],...)#,lwd=1.5)
        text(x=(ngene+1000),y=tabR[i*npt,j],labels=i,col=mycol_sim[i])
        #prédiction e* x Etot selon mutation de A
        #points(seq(1,ngene,by=pasobs),y=tabP_e[(tabP_e[,n+1]==i),j]*tabR[tabR$sim==i,2*n+1],type="l",lty=2,lwd=1.5,col=mycol_sim[i])
      }
      if (pmutA==0&add.eq==TRUE) { 
        abline(h=eq_eff_E[,j],lty=3,col=col_eff,lwd=lwd.eq)
      }
    }
  }
  
  
  
  
  ##### Etot
  if (gr.Etot.time==TRUE) {
    if (new.window==TRUE) {dev.new()}
    #par(mar=c(5,5,5,5))
    plot(seq(1,ngene,by=pasobs),tabR[1:npt,j],ylim=c(0,max(tabR[,2*n+1])), xlim=c(0,ngene+1000),type="n",xlab="Generation",ylab="Total concentration",...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
    
    for (i in which.sim){
      points(seq(1,ngene,by=pasobs),tabR[tabR$sim==i,2*n+1],type="l",col=mycol_sim[i],...)#,lwd=1.5)
      text(x=(ngene+1000),y=tabR[i*npt,2*n+1],labels=i,col=mycol_sim[i])
    }
  }
  

  
  ##### Kinetic parameters
  if (gr.kin.time==TRUE) {
    for (j in 1:n){
      if (new.window==TRUE) {dev.new()}
      #par(mar=c(5,5,5,5))
      plot(seq(1,ngene,by=pasobs),tabR[1:npt,(n+j)],xlim=c(0,ngene+1000),ylim=c(0,max(tabR[,(n+1):(2*n)])),type="n",xlab="Generation",ylab=paste("Kinetic parameters of enzyme ",j, sep=""),...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
      
      for (i in which.sim){
        points(seq(1,ngene,by=pasobs),tabR[tabR$sim==i,n+j],type="l",col=mycol_sim[i],...)#,lwd=1.5)
        text(x=(ngene+1000),y=tabR[i*npt,n+j],labels=i,col=mycol_sim[i])
      }
    }
  }
  
  
  
  ##### Activities
  if (gr.A.time==TRUE) {
    for (j in 1:n){
      if (new.window==TRUE) {dev.new()}
      #par(mar=c(5,5,5,5))
      plot(seq(1,ngene,by=pasobs),tabR[1:npt,(2*n+3+j)],xlim=c(0,ngene+1000),ylim=c(0,max(tabR[,(2*n+3+1):(2*n+3+n)])),type="n",xlab="Generation",ylab=paste("Activity of enzyme ",j, sep=""),...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
      
      for (i in which.sim){
        points(seq(1,ngene,by=pasobs),tabR[tabR$sim==i,2*n+3+j],type="l",col=mycol_sim[i],...)#,lwd=1.5)
        text(x=(ngene+1000),y=tabR[i*npt,2*n+3+j],labels=i,col=mycol_sim[i])
      }
    }
  }
  

  
  # #Log de la concentrations
  # for (j in 1:n){ #pour chaque enzyme
  #   x11()
  #   #par(mar=c(5,5,5,5))
  #   plot(seq(1,ngene,by=pasobs),log(tabR[1:npt,j]),ylim=c(0,max(log(tabR[,1:n]))), xlim=c(0,ngene+1000),type="n",xlab="Generation",ylab=paste("Enzyme concentration of enzyme ",j, sep=""),...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
  #   #graphe vide de : abscisse = nb generations, ordonnee = concentration de l'enzyme j en cours
  # 
  #   for (i in which.sim){ #un graphe par simulation
  #     points(seq(1,ngene,by=pasobs),log(tabR[tabR$sim==i,j]),type="l",col=mycol_sim[i],...)#)
  #     #graphe type points liés : abscisse = nb generation, ordonnee = concentration de l'enzyme j en cours
  #     text(x=(ngene+1000),y=log(tabR[i*npt,j]),labels=i,col=mycol_sim[i]) #numéro de la simulation affichée sur la queue du graphe, de la couleur de la simulation
  #     #prediction e* x Etot selon mutation de A
  #     #points(seq(1,ngene,by=pasobs),y=tabP_e[(tabP_e[,n+1]==i),j]*tabR[tabR$sim==i,2*n+1],type="l",lty=2,lwd=lwd.eq,col=mycol_sim[i]) #concentration prédite en fonction des mutations de A en pointillés, de la couleur de la simulation
  #   }
  # }
  
  
  ##### Position tau
  if (gr.tau.time==TRUE) {
    
    #if tau exists
    if (correl_fun=="RegPos"|correl_fun=="RegNeg"|correl_fun=="CRPos"|correl_fun=="CRNeg") {
      #compute bounds of tau such as 0<e<1
      bounds_tau <- apply(leg_E[,1:n],1,range_tau,B_fun)
      #compute tau values for each simulation : simul in column
      tabtau <- matrix(0,nrow=npt,ncol=nsim)
      for (i in 1:nsim) {
        last_res <- tabR[tabR$sim==i,]
        #compute all tau for current simulation i
        tabtau[,i] <- apply(last_res[,1:n],1,droite_tau,leg_E[i,1:n],B_fun)
      }

      
      # tau = f(t)
      if (new.window==TRUE) {dev.new()}
      #par(mar=c(5,5,5,5))
      plot(seq(1,ngene,by=pasobs),tabtau[,1],xlim=c(0,ngene+1000),ylim=c(0,1),type="n",xlab="Generation",ylab="Tau",...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
      for (i in which.sim){
        points(seq(1,ngene,by=pasobs),tabtau[,i],type="l",col=mycol_sim[i],...)#,lwd=1.5)
        text(x=(ngene+1000),y=tabtau[npt,i],labels=i,col=mycol_sim[i])
        
        #bounds of tau
        #abline(h=bounds_tau,lwd=2)
        #theoretical equilibrium
        #abline(h=1,lty=2,lwd=2)
        #effective equilibrium
        if (pmutA==0&add.eq==TRUE) { 
          abline(h=eq_eff_tau[i],col=mycol_sim[i],lty=3,lwd=lwd.eq)
        }
        
      }
      
      # #J = f(tau)
      # if (new.window==TRUE) {dev.new()}
      # #par(mar=c(5,5,5,5))
      # plot(last_tau,last_res[,2*n+3],"l",lty=1,lwd=1.5,xlab="Position tau",xlim=c(min(bounds_tau),max(bounds_tau,1)),ylim=c(0,max(last_res[,2*n+3])),ylab="Flux",...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
      # abline(v=bounds_tau,lwd=2)
      # abline(v=1,lty=2,lwd=2)
      # if (pmutA==0&add.eq==TRUE) { 
      #   abline(v=eq_eff_tau[i],col=1,lty=3,lwd=1.5)
      # }
    }
    
  }
  
  
  return(invisible(list(eq_th_e=eq_th_e,eq_th_r=eq_th_r,eq_eff_e=eq_eff_e,eq_eff_E=eq_eff_E,eq_eff_tau=eq_eff_tau)))
}







#Through time, colored by enzyme
#################################
#################################
#################################
#' @rdname simul.evol.graph.methods
#' @usage graph.simul.by.time.by.enz(all_res_sim,new.window=FALSE,add.eq=TRUE,which.sim=NULL,
#' gr.J.time=TRUE,gr.e.time=TRUE,gr.E.time=FALSE,gr.Etot.time=FALSE,
#' gr.kin.time=FALSE,gr.A.time=FALSE,gr.tau.time=FALSE,
#' gr.rep.time=FALSE,gr.sim.heading=FALSE,lwd.eq=1.5,...)
#' 
#' @param gr.tau.time Logical. Add graph depending on driving variable \eqn{\tau} if exists? 
#' @param gr.rep.time Logical. Add graph response coefficients in relation to time? 
#' @param gr.sim.heading Logical. Add an heading before each series of graphics corresponding to current simulation? 
#' 
#' @return Function \code{graph.simul.by.time.by.enz} returns nothing.
#' @export

graph.simul.by.time.by.enz <- function(all_res_sim,new.window=FALSE,add.eq=TRUE,which.sim=NULL,
                                       gr.J.time=TRUE,gr.e.time=TRUE,gr.E.time=FALSE,gr.Etot.time=FALSE,gr.kin.time=FALSE,gr.A.time=FALSE,
                                       gr.tau.time=FALSE,gr.rep.time=FALSE,gr.sim.heading=FALSE,lwd.eq=1.5,...) {
  #here, variable 'i' is always used to indicate simulation number and 'j' for enzyme number
  
  ####### Take interessing parameters
  
  #input parameters
  n <- all_res_sim$param$n #number of enzymes
  nsim <- all_res_sim$param$nsim #number of simulations
  Keq_fun <- all_res_sim$param$Keq #equilibrium constants
  X_fun <- all_res_sim$param$X #numerator of flux
  beta_fun <- all_res_sim$param$beta #co-regulation coefficients
  B_fun <- all_res_sim$param$B #global co-regulation coefficients
  correl_fun <- all_res_sim$param$correl #constraints on system
  N_fun <- all_res_sim$param$N #population size
  Etot_fun <- all_res_sim$param$Etot0 #initial total concentration
  pasobs <- all_res_sim$param$pasobs
  npt <- all_res_sim$param$npt #number of observations
  ngene <- npt*pasobs #number of generations
  pmutA <- all_res_sim$param$pmutA #proba mutation of kin
  same.E0 <- all_res_sim$param$same.E0
  is.random.E0 <- all_res_sim$param$is.random.E0
  same.kin0 <- all_res_sim$param$same.kin0
  is.random.kin0 <- all_res_sim$param$is.random.kin0
  
  
  #simulation results
  tabR <- all_res_sim$tabR
  tabP_e <- all_res_sim$tabP_e
  tabP_c <- all_res_sim$tabP_r
  
  #initial and final values of E, kin, A and relative concentrations e
  #matrix such as: each simulation in row, initial values from column 1 to n and final values from column n+1 to 2n
  leg_E <- cbind(all_res_sim$list_init$E0,all_res_sim$list_final$E_f)
  leg_kin <- cbind(all_res_sim$list_init$kin0,all_res_sim$list_final$kin_f)
  leg_A <- cbind(all_res_sim$list_init$A0,all_res_sim$list_final$A_f)
  leg_e <- leg_E
  for (i in 1:nsim){
    #initial relative concentrations e0
    leg_e[i,1:n] <- leg_E[i,1:n]/sum(leg_E[i,1:n])
    #final relative concentrations e_f
    leg_e[i,((n+1):(2*n))] <- leg_E[i,((n+1):(2*n))]/sum(leg_E[i,((n+1):(2*n))])
  }
  
  #theoretical equilibrium
  # A_base <- activities(kin_base,Keq)
  # prediction <- predict(A_base,typmut_E,correl)
  all_eq_th <- apply(all_res_sim$list_init$A0,1,predict_th,correl_fun,B_fun)#due to R transformation of matrix in numeric class when there is only one row (nsim=1)...
  #all_eq_th <- apply(leg_A[,1:n],1,predict_th,correl_fun,B_fun)
  eq_th_e <- NULL
  eq_th_r <- NULL
  for (i in 1:nsim) {
    eq_th_e <- rbind(eq_th_e,all_eq_th[[i]]$pred_e)
    eq_th_r <- rbind(eq_th_r,all_eq_th[[i]]$pred_r)
  }
  # enzpred <- eq_th$pred_e
  # cpred <- eq_th$pred_c
  lty_eq <- 2 #dashed lines
  
  #effective equilibrium
  if (correl_fun=="RegNeg"|correl_fun=="CRPos"|correl_fun=="CRNeg") {
    lty_eq <- 3 #dotted lines
    eq_eff_e <- NULL
    eq_eff_E <- NULL
    eq_eff_tau <- NULL
    for (i in 1:nsim) {
      eq_eff_by_sim <- predict_eff(leg_E[i,1:n],B_fun,leg_A[i,1:n],correl_fun)
      eq_eff_e <- rbind(eq_eff_e,as.vector(eq_eff_by_sim$pred_e))
      eq_eff_E <- rbind(eq_eff_E,as.vector(eq_eff_by_sim$pred_E))
      eq_eff_tau <- rbind(eq_eff_tau,as.vector(eq_eff_by_sim$pred_tau))
    }
  } else {
    #effective equilibrium does not exist
    eq_eff_e <- matrix(NA,nrow=nsim,ncol=(2*n))
    eq_eff_E <- matrix(NA,nrow=nsim,ncol=(2*n))
    eq_eff_tau <- matrix(NA,nrow=nsim,ncol=1)
  }
  
  
  
  
  #########################
  #######""" Graphics
  #############################
  
  #which simulations are taken
  if (length(which.sim)==0) {which.sim <- 1:nsim}
  
  #color vector for enzymes
  mycol_E <- matrix(0, ncol=2, nrow=n)
  mycol_E <- data.frame(mycol_E)
  for (j in 1:n) {
    #column 1 = enzyme nam "Enzyme j"
    mycol_E[j,1] <- paste("Enzyme ",j)
    #column 2 = enzyme color
    mycol_E[j,2] <- j+1
  }
  #color for theoretical equilibrium (identical for all simulation if A does not change)
  #if (pmutA==0) {col_th <- 1} else {col_th <- mycol_sim}
  #color for effective equilibrium: its changes if different A0 or E0
  #if (same.kin0==TRUE&same.E0==TRUE&pmutA==0) {col_eff <- 1} else {col_eff <- as.vector(mycol_E[,2])}
  
  #time axis: number of generations, but only 'npt' observations every 'pasobs' generations
  time.axis <- seq(1,ngene,by=pasobs)
  
  
  #####"
  #Every graph follow same scheme:
  #1. eventually, empty graph with time in x axis and interesting variable in y axis
  #2. for each enzyme
  #3. add connected points for enzyme j for current variable
  #4. eventually, add predicted values
  #5. add legend of color enzymes
  
  
  
  #Graphics for every simulations
  for (i in which.sim) {
    
    #Take data for simulation i
    last_res <- tabR[tabR$sim==i,]
    last_pred_e <- tabP_e[tabP_e$sim==i,]
    last_pred_c <- tabP_c[tabP_c$sim==i,]
    enz <- last_res[,1:n]/last_res[,(2*n+1)] #relative concentrations
    
    ###### Heading with parameters E0,A 0, B, correl for simulation i
    if (gr.sim.heading==TRUE) {
      if (new.window==TRUE) {dev.new()}
      #par(mar=c(5,5,5,5))
      param_sim <- c(paste("E0 =",leg_E[i,1:n]), paste("A0 =",leg_A[i,1:n]), paste("B =",B_fun), paste("Constraint :",correl_fun))
      plot(x=1:npt,y=1:npt,xlab="",ylab="",sub="Parameters",main=paste("Simulation ", i),xlim=c(0,1),ylim=c(0,1),type="n")
      legend("center",inset=c(0,0),col=1,legend=param_sim,lwd=rep(1.5,n),cex=1.2)
    }
    
    
    ###### Flux
    if (gr.J.time==TRUE) {
      if (new.window==TRUE) {dev.new()}
      #par(mar=c(5,5,5,5))
      plot(seq(1,ngene,by=pasobs),last_res[,(2*n+3)],"l",lty=1,xlab="Generation",xlim=c(0,ngene+10),ylim=c(0,max(last_res[,(2*n+3)])),ylab="Flux",...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
      
      #plot(seq(1,ngene,by=pasobs),tabR[tabR$sim==i,(2*n+3)],xlim=c(0,ngene+10),ylim=c(0,max(tabR[tabR$sim==i,(2*n+3)])),type="l",xlab="Generation",ylab="Flux",...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
      # x axis = nb of generations, y axis = flux, ylim = [0,max value of J for this simul]
    }
    
    
    
    
    ##### Relative concentrations
    if (gr.e.time==TRUE) {
      if (new.window==TRUE) {dev.new()}
      #par(mar=c(5,5,5,5))
      
      plot(seq(1,ngene,by=pasobs),rep(1,npt),type="n",xlab="Generation",ylab="Relative enzyme concentrations",ylim=c(0,1),xlim=c(0,ngene),...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
      # ngene/pasobs=npt points in x-axis, so y-axis = npt
      
      for (j in 1:n){ #for each enzyme
        #line for enzyme j
        points(seq(1,ngene,by=pasobs),enz[,j],col=j+1,type="l",lty=1,...)#,lwd=1.5)
        
        #predict concentrations for enzyme j
        if (pmutA!=0&add.eq==TRUE) { #depending on mutation of A
          points(seq(1,ngene,by=pasobs),last_pred_e[,j], col=j+1, type="l", lty=lty_eq,lwd=lwd.eq)#lwd=1.5)
        }
        if (pmutA==0&add.eq==TRUE) { 
          abline(h=eq_th_e[i,j],lty=2,col=j+1,lwd=lwd.eq)
          abline(h=eq_eff_e[i,j],lty=3,col=j+1,lwd=lwd.eq)
        }
      }
      legend("topright",inset=c(0,0),col=mycol_E[,2],legend=mycol_E[,1],lwd=rep(1.5,n),cex=1.2)
      
    }
    
    
    
    ##### Absolute concentrations
    if (gr.E.time==TRUE) {
      if (new.window==TRUE) {dev.new()}
      #par(mar=c(5,5,5,5))
      plot(seq(1,ngene,by=pasobs),rep(1,npt),type="n",xlab="Generation",ylab="Absolute enzyme concentrations",ylim=c(0,max(last_res[,1:n])),xlim=c(0,ngene),...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
      for (j in 1:n){
        points(seq(1,ngene,by=pasobs),last_res[,j],col=j+1,type="l",lty=1,...)#,lwd=1.5)
        
        if (pmutA==0&add.eq==TRUE) { 
          abline(h=eq_eff_E[i,j],lty=3,col=j+1,lwd=lwd.eq)
        }
      }
      legend("topleft",inset=c(0,0),col=mycol_E[,2],legend=mycol_E[,1],lwd=rep(1.5,n),cex=1.2)
    }
        

    
    
    
    
    ##### Etot
    if (gr.Etot.time==TRUE) {
      if (new.window==TRUE) {dev.new()}
      #par(mar=c(5,5,5,5))
      plot(seq(1,ngene,by=pasobs),last_res[,(2*n+1)],"l",lty=1,lwd=1.5,xlab="Generation",xlim=c(0,ngene+10),ylim=c(0,max(last_res[,(2*n+1)])),ylab="Total concentration",...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
    }
    

    
    ##### Kinetic parameters
    if (gr.kin.time==TRUE) {
      if (new.window==TRUE) {dev.new()}
      #par(mar=c(5,5,5,5))
      plot(seq(1,ngene,by=pasobs),rep(1,npt),type="n",xlab="Generation",ylab="Kinetic parameter",xlim=c(0,ngene),ylim=c(0,max(last_res[,(n+1):(2*n)])),...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
      for (j in 1:n){
        points(seq(1,ngene,by=pasobs),last_res[,n+j],type="l",col=j+1,lty=1,...)#,lwd=1.5)
      }
      legend("topleft",col=mycol_E[,2],legend=mycol_E[,1],lwd=rep(1.5,n),cex=1.2)  
    }
    
 
    
    ##### Activities
    if (gr.A.time==TRUE) {
      if (new.window==TRUE) {dev.new()}
      #par(mar=c(5,5,5,5))
      plot(seq(1,ngene,by=pasobs),rep(1,npt),type="n",xlab="Generation",ylab="Activity",xlim=c(0,ngene),ylim=c(0,max(last_res[,(2*n+3+1):(2*n+3+n)])),...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
      for (j in 1:n){
        points(seq(1,ngene,by=pasobs),last_res[,2*n+3+j],type="l",col=j+1,lty=1,...)#,lwd=1.5)
      }
      legend("topleft",col=mycol_E[,2],legend=mycol_E[,1],lwd=rep(1.5,n),cex=1.2)
    }
    
    
    ##### Response coefficients
    if (gr.rep.time==TRUE) {
      cont <- NULL
      for (k in 1:nrow(last_res)) {
        #computes response coefficients for each step of simulation
        rep_row <- coef_rep(last_res[k,1:n],last_res[k,(2*n+4):(3*n)],correl_fun,beta_fun)
        cont <- rbind(cont,rep_row)
      }
      if (new.window==TRUE) {dev.new()}
      #par(mar=c(5,5,5,5))
      plot(seq(1,ngene,by=pasobs),rep(1,npt),type="n",xlab="Generation",ylab="Response coefficients",ylim=c(0,1),xlim=c(0,ngene),...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
      for (j in 1:n){
        points(seq(1,ngene,by=pasobs),cont[,j],type="l",col=j+1,lty=1,...)#,lwd=1.5)
        if (pmutA!=0&add.eq==TRUE) {
          points(seq(1,ngene,by=pasobs),last_pred_c[,j], col=j+1, type="l", lty=2,lwd=lwd.eq)
        }
        if (pmutA==0&add.eq==TRUE) { 
          abline(h=eq_th_r[i,j],col=j+1,lty=2,lwd=lwd.eq)
        }
      }
      legend("topright",inset=c(0,0),col=mycol_E[,2],legend=mycol_E[,1],lwd=rep(1.5,n),cex=1.2)
      
    }
    
    
    ##### Position tau
    #if tau exists
    if (correl_fun=="RegPos"|correl_fun=="RegNeg"|correl_fun=="CRPos"|correl_fun=="CRNeg") {
      #compute bounds of tau such as 0<e<1
      bounds_tau <- range_tau(leg_E[i,1:n],B_fun)
      #compute tau values for each row
      last_tau <- apply(last_res[,1:n],1,droite_tau,leg_E[i,1:n],B_fun)
      
      
      if (gr.tau.time==TRUE) {
        # tau = f(t)
        if (new.window==TRUE) {dev.new()}
        #par(mar=c(5,5,5,5))
        plot(seq(1,ngene,by=pasobs),last_tau,"l",lty=1,xlab="Generation",xlim=c(0,ngene),ylim=c(min(bounds_tau),max(bounds_tau,1)),ylab="Position tau",...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
        #bounds of tau
        abline(h=bounds_tau,lwd=1)
        #theoretical equilibrium
        abline(h=1,lty=2,lwd=lwd.eq)
        #effective equilibrium
        if (pmutA==0&add.eq==TRUE) { 
          abline(h=eq_eff_tau[i],col=1,lty=3,lwd=lwd.eq)
        }
        
        #J = f(tau)
        if (new.window==TRUE) {dev.new()}
        #par(mar=c(5,5,5,5))
        plot(last_tau,last_res[,2*n+3],"l",lty=1,xlab="Position tau",xlim=c(min(bounds_tau),max(bounds_tau,1)),ylim=c(0,max(last_res[,2*n+3])),ylab="Flux",...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
        abline(v=bounds_tau,lwd=1)
        abline(v=1,lty=2,lwd=lwd.eq)
        if (pmutA==0&add.eq==TRUE) { 
          abline(v=eq_eff_tau[i],col=1,lty=3,lwd=lwd.eq)
        }
      }
      
    }
    
    
    
    
    #End for simulation i
  }
  
  
}



#Other graphics, colored by simulation
#################################
#################################
#################################
#' @rdname simul.evol.graph.methods
#' @usage graph.simul.others.by.sim(all_res_sim,new.window=FALSE,add.eq=TRUE,which.sim=NULL,
#' gr.Ef.E0=FALSE, gr.Af.A0=FALSE, gr.Ef.Af=FALSE, gr.J.A.E=TRUE,
#' gr.J.e=FALSE, gr.J.E=FALSE, gr.J.A=FALSE, ...)
#' 
#' @param gr.Ef.E0 Logical. Add graph of final concentrations (absolute E and relative e) depending its initial value? 
#' @param gr.Af.A0 Logical. Add graph of final activities (resp. kinetic parameters) depending its initial value? 
#' @param gr.Ef.Af Logical. Add graph of final concentrations (absolute E and relative e) depending on final activities A? 
#' @param gr.J.A.E Logical. Add 3D-graph of flux J depending on concentrations E and activities A? 
#' @param gr.J.e Logical. Add graph of flux J depending on \emph{relative} concentrations e?
#' @param gr.J.E Logical. Add graph of flux J depending on \emph{absolute} concentrations E? 
#' @param gr.J.A Logical. Add graph of of flux J depending on activities A? 
#' 
#' @return Function \code{graph.simul.others.by.sim} returns nothing.
#' @export


graph.simul.others.by.sim <- function(all_res_sim,new.window=FALSE,add.eq=TRUE,which.sim=NULL,
                                      gr.Ef.E0=FALSE, gr.Af.A0=FALSE, gr.Ef.Af=FALSE, gr.J.A.E=TRUE,
                                      gr.J.e=FALSE,gr.J.E=FALSE, gr.J.A=FALSE,...) {
  #here, variable 'i' is always used to indicate simulation number and 'j' for enzyme number
  
  ####### Take interessing parameters
  
  #input parameters
  n <- all_res_sim$param$n #number of enzymes
  nsim <- all_res_sim$param$nsim #number of simulations
  Keq_fun <- all_res_sim$param$Keq #equilibrium constants
  X_fun <- all_res_sim$param$X #numerator of flux
  beta_fun <- all_res_sim$param$beta #co-regulation coefficients
  B_fun <- all_res_sim$param$B #global co-regulation coefficients
  correl_fun <- all_res_sim$param$correl #constraints on system
  N_fun <- all_res_sim$param$N #population size
  pasobs <- all_res_sim$param$pasobs
  npt <- all_res_sim$param$npt #number of observations
  ngene <- npt*pasobs #number of generations
  pmutA <- all_res_sim$param$pmutA #proba mutation of kin
  same.E0 <- all_res_sim$param$same.E0
  is.random.E0 <- all_res_sim$param$is.random.E0
  same.kin0 <- all_res_sim$param$same.kin0
  is.random.kin0 <- all_res_sim$param$is.random.kin0
  
  
  #simulation results
  tabR <- all_res_sim$tabR
  tabP_e <- all_res_sim$tabP_e
  tabP_c <- all_res_sim$tabP_r
  
  #initial and final values of E, kin, A and relative concentrations e
  #matrix such as: each simulation in row, initial values from column 1 to n and final values from column n+1 to 2n
  leg_E <- cbind(all_res_sim$list_init$E0,all_res_sim$list_final$E_f)
  leg_kin <- cbind(all_res_sim$list_init$kin0,all_res_sim$list_final$kin_f)
  leg_A <- cbind(all_res_sim$list_init$A0,all_res_sim$list_final$A_f)
  leg_e <- leg_E
  for (i in 1:nsim){
    #initial relative concentrations e0
    leg_e[i,1:n] <- leg_E[i,1:n]/sum(leg_E[i,1:n])
    #final relative concentrations e_f
    leg_e[i,((n+1):(2*n))] <- leg_E[i,((n+1):(2*n))]/sum(leg_E[i,((n+1):(2*n))])
  }
  
  #theoretical equilibrium
  # A_base <- activities(kin_base,Keq)
  # prediction <- predict(A_base,typmut_E,correl)
  all_eq_th <- apply(all_res_sim$list_init$A0,1,predict_th,correl_fun,B_fun)#due to R transformation of matrix in numeric class when there is only one row (nsim=1)...
  #all_eq_th <- apply(leg_A[,1:n],1,predict_th,correl_fun,B_fun)
  eq_th_e <- NULL
  eq_th_r <- NULL
  for (i in 1:nsim) {
    eq_th_e <- rbind(eq_th_e,all_eq_th[[i]]$pred_e)
    eq_th_r <- rbind(eq_th_r,all_eq_th[[i]]$pred_r)
  }
  # enzpred <- eq_th$pred_e
  # cpred <- eq_th$pred_c
  lty_eq <- 2 #dashed lines
  
  #effective equilibrium
  if (correl_fun=="RegNeg"|correl_fun=="CRPos"|correl_fun=="CRNeg") {
    lty_eq <- 3 #dotted lines
    eq_eff_e <- NULL
    eq_eff_E <- NULL
    eq_eff_tau <- NULL
    for (i in 1:nsim) {
      eq_eff_by_sim <- predict_eff(leg_E[i,1:n],B_fun,leg_A[i,1:n],correl_fun)
      eq_eff_e <- rbind(eq_eff_e,as.vector(eq_eff_by_sim$pred_e))
      eq_eff_E <- rbind(eq_eff_E,as.vector(eq_eff_by_sim$pred_E))
      eq_eff_tau <- rbind(eq_eff_tau,as.vector(eq_eff_by_sim$pred_tau))
    }
  } else {
    #effective equilibrium does not exist
    eq_eff_e <- matrix(NA,nrow=nsim,ncol=(2*n))
    eq_eff_E <- matrix(NA,nrow=nsim,ncol=(2*n))
    eq_eff_tau <- matrix(NA,nrow=nsim,ncol=1)
  }
  
  
  
  
  #########################
  #######""" Graphics
  #############################
  
  #which simulations will be represented. If NULL, all
  if (length(which.sim)==0) {which.sim <- 1:nsim}
  
  #color vector for simulations
  mycol_sim <- rainbow(nsim)
  #color vector for enzymes
  mycol_E <- matrix(0, ncol=2, nrow=n)
  mycol_E <- data.frame(mycol_E)
  for (j in 1:n) {
    #column 1 = enzyme nam "Enzyme j"
    mycol_E[j,1] <- paste("Enzyme ",j)
    #column 2 = enzyme color
    mycol_E[j,2] <- j+1
  }
  #color for theoretical equilibrium (identical for all simulation if A does not change)
  #if (pmutA==0) {col_th <- 1} else {col_th <- mycol_sim}
  #color for effective equilibrium: its changes if different A0 or E0
  if (same.kin0==TRUE&same.E0==TRUE&pmutA==0) {col_eff <- 1} else {col_eff <- mycol_sim}
  
  
  
  
  if (gr.Ef.E0==TRUE) {
    
    ###### Ef in relation to E0
    if (new.window==TRUE) {dev.new()}
    #par(mar=c(5,5,5,5))
    plot(which.sim,which.sim,xlim=c(0,max(leg_E[which.sim,1:n])),ylim=c(0,max(leg_E[which.sim,(n+1):(2*n)])),type="n",xlab="Initial concentration",ylab="Final concentration",...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
    #empty graph of: x axis = E0, y axis = Ef
    for (j in 1:n) {
      #points(x=leg_E[which.sim,j],y=leg_E[which.sim,(n+j)],type="p",col=j+1)
      #text of simulation number, colored according to enzyme
      text(x=leg_E[which.sim,j],y=leg_E[which.sim,(n+j)],labels=which.sim,col=j+1)
    }
    legend("topright",inset=c(0,0),col=mycol_E[,2],legend=mycol_E[,1],lwd=rep(1.5,n),cex=1.2)
    
    
    
    ###### ef in relation to e0
    if (new.window==TRUE) {dev.new()}
    #par(mar=c(5,5,5,5))
    plot(which.sim,which.sim,xlim=c(0,1),ylim=c(0,1),type="n",xlab="Initial relative concentration",ylab="Final relative concentration",...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
    #graphe vide de : x axis = e0, y axis = ef
    for (j in 1:n) {
      #points(x=leg_e[which.sim,j],y=leg_e[which.sim,(n+j)],type="p",col=j+1)
      text(x=leg_e[which.sim,j],y=leg_e[which.sim,(n+j)],labels=which.sim,col=j+1)
    }
    legend("topright",inset=c(0,0),col=mycol_E[,2],legend=mycol_E[,1],lwd=rep(1.5,n),cex=1.2)
    
  }
  
  
  if (gr.Af.A0==TRUE) {
    
    ##### kinf en fonction de kin0
    if (new.window==TRUE) {dev.new()}
    #par(mar=c(5,5,5,5))
    plot(which.sim,which.sim,xlim=c(0,max(leg_kin[which.sim,1:n])),ylim=c(0,max(leg_kin[which.sim,(n+1):(2*n)])),type="n",xlab="Initial kinetic parameter",ylab="Final kinetic parameter",...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
    #graphe vide de : x axis = A0, y axis = Af
    for (j in 1:n) {
      #points(x=leg_kin[which.sim,j],y=leg_kin[which.sim,(n+j)],type="p",col=j+1)
      text(x=leg_kin[which.sim,j],y=leg_kin[which.sim,(n+j)],labels=which.sim,col=j+1)
    }
    legend("topright",inset=c(0,0),col=mycol_E[,2],legend=mycol_E[,1],lwd=rep(1.5,n),cex=1.2)
    
    
    ##### Af en fonction de A0
    if (new.window==TRUE) {dev.new()}
    #par(mar=c(5,5,5,5))
    plot(which.sim,which.sim,xlim=c(0,max(leg_A[which.sim,1:n])),ylim=c(0,max(leg_A[which.sim,(n+1):(2*n)])),type="n",xlab="Initial activity",ylab="Final activity",...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
    #graphe vide de : x axis = A0, y axis = Af
    for (j in 1:n) {
      #points(x=leg_A[which.sim,j],y=leg_A[which.sim,(n+j)],type="p",col=j+1)
      text(x=leg_A[which.sim,j],y=leg_A[which.sim,(n+j)],labels=which.sim,col=j+1)
    }
    legend("topright",inset=c(0,0),col=mycol_E[,2],legend=mycol_E[,1],lwd=rep(1.5,n),cex=1.2)
    
  }
  
  
  if (gr.Ef.Af==TRUE) {
    
    ##### Concentration relative finale ef en fonction de Af
    if (new.window==TRUE) {dev.new()}
    #par(mar=c(5,5,5,5))
    plot(which.sim,which.sim,xlim=c(0,max(leg_A[which.sim,(n+1):(2*n)])),ylim=c(0,1),type="n",xlab="Final activity",ylab="Final relative concentration",main="",cex.main=2,cex.lab=1.6,cex.axis=1.3)
    #graphe vide de : x axis = Af, y axis = Ef
    for (j in 1:n) {#pour chaque enzyme
      #points(x=leg_A[which.sim,n+j],y=leg_e[which.sim,(n+j)],type="p",col=j+1)
      text(x=leg_A[which.sim,n+j],y=leg_e[which.sim,(n+j)],labels=which.sim,col=j+1) #graphe de Ef en fonction de Af, numérotée pour chaque simulation, colorée selon l'enzyme
    }
    legend("topright",inset=c(0,0),col=mycol_E[,2],legend=mycol_E[,1],lwd=rep(1.5,n),cex=1.2)
    
    
    ##### Ef en fonction de Af
    if (new.window==TRUE) {dev.new()}
    #par(mar=c(5,5,5,5))
    plot(which.sim,which.sim,xlim=c(0,max(leg_A[which.sim,(n+1):(2*n)])),ylim=c(0,max(leg_E[which.sim,(n+1):(2*n)])),type="n",xlab="Final activity",ylab="Final concentration",...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
    #graphe vide de : x axis = Af, ordonne = Ef
    for (j in 1:n) {#pour chaque enzyme
      #points(x=leg_A[which.sim,n+j],y=leg_E[which.sim,(n+j)],type="p",col=j+1)
      text(x=leg_A[which.sim,n+j],y=leg_E[which.sim,(n+j)],labels=which.sim,col=j+1) #graphe de Ef en fonction de Af, numérotée pour chaque simulation, colorée selon l'enzyme
    }
    legend("topright",inset=c(0,0),col=mycol_E[,2],legend=mycol_E[,1],lwd=rep(1.5,n),cex=1.2)
    
  }
  
  
  
  if (gr.J.A.E==TRUE) {
    
    ##### J en fonction de A et E
    if (new.window==TRUE) {dev.new()}
    #par(mar=c(5,5,5,5))
    s3d<-scatterplot3d(x=1:npt,y=1:npt,z=1:npt,xlab="Concentration",ylab="Activity",zlab="Flux",xlim=c(0,max(tabR[,1:n])),ylim=c(0,max(tabR[,(n+1):(2*n)])),zlim=c(0,max(tabR[,2*n+3])),type="n")
    #empty 3D-plot of E in x, A in y and J in z
    for (j in 1:n){
      s3d$points3d(x=tabR[,j],y=tabR[,n+j],z=tabR[,2*n+3],col=j+1,cex=0.2)
    }
  }
  
  if (gr.J.e==TRUE) {
    
    ####### 2D plot J en fonction de ej, colored by sim
    for (j in 1:n){
      if (new.window==TRUE) {dev.new()}
      #par(mar=c(5,5,5,5))
      plot(x=1:npt,y=1:npt,xlab=paste("Relative concentration of enzyme ",j,sep=""),ylab="Flux",xlim=c(0,1),ylim=c(0,max(tabR[,2*n+3])),type="n",...)#)
      for (i in which.sim) {#pour chaque simulation
        points(x=tabR[tabR$sim==i,j]/tabR[tabR$sim==i,2*n+1],y=tabR[tabR$sim==i,2*n+3],col=mycol_sim[i],cex=0.2)
      }
    }
  }
  
  if (gr.J.E==TRUE) {
    
    ####### 2D plot J en fonction de ej, colored by sim
    for (j in 1:n){
      if (new.window==TRUE) {dev.new()}
      #par(mar=c(5,5,5,5))
      plot(x=1:npt,y=1:npt,xlab=paste("Concentration of enzyme ",j,sep=""),ylab="Flux",xlim=c(0,1),ylim=c(0,max(tabR[,2*n+3])),type="n",...)#)
      for (i in which.sim) {#pour chaque simulation
        points(x=tabR[tabR$sim==i,j],y=tabR[tabR$sim==i,2*n+3],col=mycol_sim[i],cex=0.2)
      }
    }
  }
  
  if (gr.J.A==TRUE) {
    ###### J fct de A par sim
    for (j in 1:n){
      if (new.window==TRUE) {dev.new()}
      #par(mar=c(5,5,5,5))
      plot(x=1:npt,y=1:npt,xlab=paste("Activity of enzyme ",j,sep=""),ylab="Flux",xlim=c(0,max(tabR[,(n+1):(2*n)])),ylim=c(0,max(tabR[,2*n+3])),type="n",...)#)
      for (i in which.sim) {#pour chaque simulation
        points(x=tabR[tabR$sim==i,n+j],y=tabR[tabR$sim==i,2*n+3],col=mycol_sim[i],cex=0.2)
      }
    }
    
  }
  
}







#Through time, Range of Neutral Variation
#################################
#################################
#################################
#' @rdname simul.evol.graph.methods
#' @usage graph.simul.by.time.RNV(all_res_sim,new.window=FALSE,add.eq=TRUE,which.sim=NULL,
#' gr.RNV.E=TRUE,gr.RNV.size=FALSE,gr.RNV.delta=FALSE,
#' gr.RNV.J=FALSE,zoom.RNV.J=NULL,gr.sim.heading=FALSE,
#' col_RNV=c("grey60","grey80"),lty_RNV=c("dashed","longdash"),lwd.eq=1.5,...)
#' 
#' @param gr.RNV.E,gr.RNV.size,gr.RNV.delta Logical. Add graph concentrations \code{E} with its RNV / RNV size and RNV size divided by Etot / apparent mutation effect \code{delta} corresponding to RNV in relation to time?
#' 
#' @param gr.RNV.J Logical. Add graph of flux depending on time and depending on concentrations with RNV and neutral zone? 
#' @param zoom.RNV.J Numeric vector of length 2, corresponding to \code{ylim} of graphics flux with RNV. If \code{NULL} (by default), \code{ylim} is zoomed around maximal flux of current simulation.
#' If it is nor \code{NULL} nor vector of length 2, there is no zoom on flux.
#' @param col_RNV,lty_RNV Vector of length 2, for color (resp. lty, see plot function) of RNV lines. First element correspond to inferior bounds and second one to superior bounds of RNV.
#'
#' @details
#' \emph{col_RNV and lty_RNV}
#' 
#' Vector of length 2, for color (resp. lty, see \code{plot} function) of RNV lines. First element correspond to inferior bounds and second one to superior bounds of RNV.
#' 
#' These parameters are only available for plot \code{gr.RNV.E} and \code{gr.RNV.J}.
#'    
#' @return Function \code{graph.simul.by.time.RNV} returns invisible list of 2 elements:
#' \itemize{
#'    \item \code{$RNV_all_sim}: List of \code{nsim} elements (which is number of simulation). Every element is the output of function \code{\link{RNV.for.simul}} for corresponding simulation.
#'    If simulation \emph{i} is not contained in \code{which.sim}, \code{$RNV_all_sim[[i]]} is \code{NULL}. 
#'    \item \code{$RNV_mean_size}: Numeric matrix of \code{n+2} columns. Row number is between \code{nsim} and \code{2*nsim}, depending on applied constraint.
#'    \code{n} first columns correspond to RNV mean size for corresponding enzyme for last half simulation; column \code{n+1} indicates simulation number and column \code{n+2} RNV number (between 1 and 2).
#' } 
#' 
#' @export

graph.simul.by.time.RNV <- function(all_res_sim,new.window=FALSE,add.eq=TRUE,which.sim=NULL,
                                    gr.RNV.E=TRUE,gr.RNV.size=FALSE,gr.RNV.delta=FALSE,
                                    gr.RNV.J=FALSE,zoom.RNV.J=NULL,gr.sim.heading=FALSE,
                                    col_RNV=c("grey60","grey80"),lty_RNV=c("dashed","longdash"),lwd.eq=1.5,...) {
  #here, variable 's' is always used to indicate simulation number and 'i' for mutant enzyme number
  
  ####### Take interessing parameters
  
  #input parameters
  n <- all_res_sim$param$n #number of enzymes
  nsim <- all_res_sim$param$nsim #number of simulations
  Keq_fun <- all_res_sim$param$Keq #equilibrium constants
  X_fun <- all_res_sim$param$X #numerator of flux
  beta_fun <- all_res_sim$param$beta #co-regulation coefficients
  B_fun <- all_res_sim$param$B #global co-regulation coefficients
  correl_fun <- all_res_sim$param$correl #constraints on system
  N_fun <- all_res_sim$param$N #population size
  Etot_fun <- all_res_sim$param$Etot0 #initial total concentration
  pasobs <- all_res_sim$param$pasobs
  npt <- all_res_sim$param$npt #number of observations
  ngene <- npt*pasobs #number of generations
  pmutA <- all_res_sim$param$pmutA #proba mutation of kin
  same.E0 <- all_res_sim$param$same.E0
  is.random.E0 <- all_res_sim$param$is.random.E0
  same.kin0 <- all_res_sim$param$same.kin0
  is.random.kin0 <- all_res_sim$param$is.random.kin0
  #to keep track of user choice
  choice.zoom.J <- zoom.RNV.J
  
  
  #simulation results
  tabR <- all_res_sim$tabR
  tabP_e <- all_res_sim$tabP_e
  tabP_c <- all_res_sim$tabP_r
  
  #initial and final values of E, kin, A and relative concentrations e
  #matrix such as: each simulation in row, initial values from column 1 to n and final values from column n+1 to 2n
  leg_E <- cbind(all_res_sim$list_init$E0,all_res_sim$list_final$E_f)
  leg_kin <- cbind(all_res_sim$list_init$kin0,all_res_sim$list_final$kin_f)
  leg_A <- cbind(all_res_sim$list_init$A0,all_res_sim$list_final$A_f)
  leg_e <- leg_E
  for (i in 1:nsim){
    #initial relative concentrations e0
    leg_e[i,1:n] <- leg_E[i,1:n]/sum(leg_E[i,1:n])
    #final relative concentrations e_f
    leg_e[i,((n+1):(2*n))] <- leg_E[i,((n+1):(2*n))]/sum(leg_E[i,((n+1):(2*n))])
  }
  
  #theoretical equilibrium
  # A_base <- activities(kin_base,Keq)
  # prediction <- predict(A_base,typmut_E,correl)
  all_eq_th <- apply(all_res_sim$list_init$A0,1,predict_th,correl_fun,B_fun)#due to R transformation of matrix in numeric class when there is only one row (nsim=1)...
  #all_eq_th <- apply(leg_A[,1:n],1,predict_th,correl_fun,B_fun)
  eq_th_e <- NULL
  eq_th_r <- NULL
  for (i in 1:nsim) {
    eq_th_e <- rbind(eq_th_e,all_eq_th[[i]]$pred_e)
    eq_th_r <- rbind(eq_th_r,all_eq_th[[i]]$pred_r)
  }
  # enzpred <- eq_th$pred_e
  # cpred <- eq_th$pred_c
  lty_eq <- 2 #dashed lines
  
  #effective equilibrium
  if (correl_fun=="RegNeg"|correl_fun=="CRPos"|correl_fun=="CRNeg") {
    lty_eq <- 3 #dotted lines
    eq_eff_e <- NULL
    eq_eff_E <- NULL
    eq_eff_tau <- NULL
    eq_eff_J <- NULL
    for (i in 1:nsim) {
      eq_eff_by_sim <- predict_eff(leg_E[i,1:n],B_fun,leg_A[i,1:n],correl_fun)
      eq_eff_e <- rbind(eq_eff_e,as.vector(eq_eff_by_sim$pred_e))
      eq_eff_E <- rbind(eq_eff_E,as.vector(eq_eff_by_sim$pred_E))
      eq_eff_tau <- rbind(eq_eff_tau,as.vector(eq_eff_by_sim$pred_tau))
      eq_eff_J <- rbind(eq_eff_J,flux(as.vector(eq_eff_by_sim$pred_E),leg_A[i,1:n],X_fun))
    }
  } else {
    #effective equilibrium does not exist
    eq_eff_e <- matrix(NA,nrow=nsim,ncol=(2*n))
    eq_eff_E <- matrix(NA,nrow=nsim,ncol=(2*n))
    eq_eff_tau <- matrix(NA,nrow=nsim,ncol=1)
    eq_eff_J <- rep(NA,nsim)
  }
  
  
  
  
  #########################
  #######""" Graphics
  #############################
  
  #which simulations are taken
  if (length(which.sim)==0) {which.sim <- 1:nsim}
  
  #color vector for enzymes
  mycol_E <- matrix(0, ncol=2, nrow=n)
  mycol_E <- data.frame(mycol_E)
  for (j in 1:n) {
    #column 1 = enzyme nam "Enzyme j"
    mycol_E[j,1] <- paste("Enzyme ",j)
    #column 2 = enzyme color
    mycol_E[j,2] <- j+1
  }
  #color for theoretical equilibrium (identical for all simulation if A does not change)
  #if (pmutA==0) {col_th <- 1} else {col_th <- mycol_sim}
  #color for effective equilibrium: its changes if different A0 or E0
  #if (same.kin0==TRUE&same.E0==TRUE&pmutA==0) {col_eff <- 1} else {col_eff <- as.vector(mycol_E[,2])}
  
  #color for RNV bounds
  #list elements = inf or sup bounds
  # df <- data.frame(col=c(0,0),lty=c(0,0))
  # mycol_RNV <- list(inf=df,sup=df)
  # #inside data.frame columns = color and lty for RNV bounds
  # # colnames(mycol_RNV$inf) <- c("col","lty")
  # # colnames(mycol_RNV$sup) <- c("col","lty")
  # #inside data.frame rows = which RNV (row 1 = nearest RNV; row 2 = farthest RNV)
  # mycol_RNV$inf$col <- c("grey60") #inf bounds in dark grey
  # mycol_RNV$sup$col <- c("grey80") #sup bounds in light grey
  # mycol_RNV$sup$lty <- c(3) #sup bounds in dotted line
  # mycol_RNV$inf$lty <- c(2,4) #near inf bounds in dashed line ; far inf bounds in dotted-dashed line
  mycol_RNV <- list(col_RNV,lty_RNV)
  
  #time axis: number of generations, but only 'npt' observations every 'pasobs' generations
  time.axis <- seq(1,ngene,by=pasobs)
  
  
  #####"
  #Every graph follow same scheme:
  #1. eventually, empty graph with time in x axis and interesting variable in y axis
  #2. for each enzyme
  #3. add connected points for enzyme j for current variable
  #4. eventually, add predicted values
  #5. add legend of color enzymes
  
  
  
  #Graphics for every simulations
  
  #save RNV for all simulations
  RNV_all_sim <- vector("list",length=nsim)
  #take mean RNV size for all simul
  RNV_mean_size <- NULL
  
  for (s in which.sim) {
    
    #Take result for simulation s
    res_sim <- tabR[tabR$sim==s,]
    #compute RNV elements
    RNV_all <- RNV.for.simul(res_sim,n,N_fun,correl_fun,beta_fun)
    #save
    RNV_all_sim[[s]] <- RNV_all
    
    nb_RNV <- RNV_all$nb_RNV
    
    RNV_dis_delta_graph <- RNV_all$RNV_delta
    RNV_dis_enz_graph <- RNV_all$RNV_enz
    RNV_dis_size <- RNV_all$RNV_size
    RNV_dis_size_divEtot <- RNV_all$RNV_size_divEtot
    RNV_proxy <- RNV_all$RNV_proxy
    RNV_proxy_divEtot <- RNV_all$RNV_proxy_divEtot
    RNV_dis_flux <- RNV_all$RNV_flux
    limits_NZ <- RNV_all$limits_NZ
    
    #save mean RNV
    #add simul number 's' and RNV number
    RNv_promean <- cbind(RNV_proxy,s,1:nb_RNV)
    RNV_mean_size <- rbind(RNV_mean_size,RNv_promean)
    
    
    # #Take concentrations and total concentration for simulation s
    tabE <- res_sim[,c(1:n,2*n+1)]
    #tabE <- cbind(tabE,tabR[tabR$sim==s,(2*n+1)]) #total concentration
    #Take activities for simulation s
    tabA <- res_sim[,(2*n+4):(3*n+3)]
    tabe_pred <- tabP_e[tabP_e$sim==s,1:n]

    # 
    # 
    # #how many RNV
    # if (correl_fun=="SC"|correl_fun=="RegPos") {
    #   nb_RNV <- 1
    # } else {
    #   nb_RNV <- 2
    # }
    # 
    # 
    # #########" Computes RNV
    # 
    # #list of n elements, one by enzyme. For each enzyme, rows = simulation and columns = delta corresponding to inf bounds and sup bounds of RNV (x2 if there is a 2nd RNV)
    # RNV_dis_delta_graph <- vector("list",length=n)
    # #idem, but for enzyme concentrations + delta
    # RNV_dis_enz_graph <- vector("list",length=n)
    # 
    # #for each result of current simulation s
    # for (r in 1:nrow(tabE)) {
    #   #extract resident concentrations for this line
    #   E_resid <- as.matrix(tabE[r,1:n])
    #   A_resid <- as.matrix(tabA[r,1:n])
    # 
    #   # compute delta at limits of neutral zone
    #   delta_RNV_all <- RNV.delta.all.enz(E_resid,A_resid,N_fun,correl_fun,beta_fun)
    #   #delta_all <- delta.RNV.all(E_resid,A_base,N,correl,beta,B,n,lim_calc) #E_fun,A_fun,N_fun,correl_fun,beta_fun,B_fun,n_fun,d_lim
    #   
    #   #for each enzyme
    #   for (i in 1:n) {
    #     #take corresponding delta at bounds of RNV
    #     RNV_dis_delta_graph[[i]] <- rbind(RNV_dis_delta_graph[[i]],delta_RNV_all[[i]])
    #     #add delta to concentrations to have corresponding "mutant" concentration at bounds of RNV
    #     RNV_dis_enz_graph[[i]] <- rbind(RNV_dis_enz_graph[[i]],E_resid[i] + delta_RNV_all[[i]])
    #     
    #   }
    # }
    # 
    # 
    # ####### RNV size
    # ##list of one or two elemnts (number of RNV), and in each element, matrix of n columns and npt rows (nb points for current simulation)
    # #size of RNV (delta_sup - delta_inf)
    # RNV_dis_size <- vector("list",length=nb_RNV)
    # #RNV size divided by resident total concentration
    # RNV_dis_size_divEtot <- vector("list",length=n)
    # 
    # ## matrix of one or two rows (number of RNV) and n columns
    # #mean RNV size
    # RNV_proxy <- matrix(0,nrow=nb_RNV,ncol=n)
    # #mean RNV size divided by resident total concentration
    # RNV_proxy_divEtot <- matrix(0,nrow=nb_RNV,ncol=n)
    # 
    # #for each RNV
    # for (nor in 1:nb_RNV) {
    #   #for each enzyme
    #   for (i in 1:n) {
    #     #take RNV size for all result 'r'
    #     RNV_size_vec_r <- NULL
    #     
    #     #for each result
    #     for (r in 1:nrow(tabE)) {
    #       
    #       #superior and inferior bounds of current RNV 'nor' and current enzyme 'i' and current result 'r'
    #       RNV_sup_bounds_delta <- RNV_dis_delta_graph[[i]][r,2*nor]
    #       RNV_inf_bounds_delta <- RNV_dis_delta_graph[[i]][r,2*nor-1]
    #       
    #       #if there is two RNV and sup bounds is not accessible
    #       if (nb_RNV==2 & is.na(RNV_sup_bounds_delta)) {
    #         # RNV size = | inf bounds (RNV 2) - inf bounds (RNV 1) |
    #         RNV_size_num_r <- abs(RNV_dis_delta_graph[[i]][r,3] - RNV_dis_delta_graph[[i]][r,1])
    #         #RNV_dis_size[[nor]] <- cbind(RNV_dis_size[[nor]], abs(RNV_dis_delta_graph[[i]][,3] - RNV_dis_delta_graph[[i]][,1]))
    #       } else {
    #         #RNV size = |sup bounds - inf bounds| for RNV 'nor'
    #         RNV_size_num_r <- abs(RNV_sup_bounds_delta - RNV_inf_bounds_delta)
    #         #RNV_dis_size[[nor]] <- cbind(RNV_dis_size[[nor]], abs(RNV_sup_bounds_delta - RNV_inf_bounds_delta))
    #       }
    #       
    #       #add this RNV size for this result 'r' with precedent
    #       RNV_size_vec_r <- c(RNV_size_vec_r,RNV_size_num_r)
    #     }
    #     
    #     #add RNV size for all rows with precedent enzyme 
    #     RNV_dis_size[[nor]] <- cbind(RNV_dis_size[[nor]], RNV_size_vec_r)
    #   }
    #   
    #   #division by Etot : tabE has same nrows than RNV_dis_size element, and R compute column by column
    #   RNV_dis_size_divEtot[[nor]] <- RNV_dis_size[[nor]]/tabE[,n+1]
    #   
    #   #mean of RNV size for half end of simul
    #   RNV_proxy[nor,] <- apply(RNV_dis_size[[nor]][round(npt/2):npt,],2,mean,na.rm=TRUE)
    #   RNV_proxy_divEtot[nor,] <- apply(RNV_dis_size_divEtot[[nor]][round(npt/2):npt,],2,mean,na.rm=TRUE)
    # }
    #print(RNV_proxy)
    
    # delta_mean <- vector("list",length=n)
    # for (i in 1:n) {
    #   bfff <- apply(RNV_dis_delta_graph[[i]][round(npt/2):npt,],2,mean,na.rm=TRUE)
    #   delta_mean[[i]] <- c(delta_mean[[i]],bfff)
    # }
    
    ###### Heading with parameters E0,A 0, B, correl for simulation s
    if (gr.sim.heading==TRUE) {
      if (new.window==TRUE) {dev.new()}
      #par(mar=c(5,5,5,5))
      param_sim <- c(paste("E0 =",leg_E[s,1:n]), paste("A0 =",leg_A[s,1:n]), paste("B =",B_fun), paste("Constraint :",correl_fun))
      plot(x=1:npt,y=1:npt,xlab="",ylab="",sub="Parameters",main=paste("Simulation ", s),xlim=c(0,1),ylim=c(0,1),type="n")
      legend("center",inset=c(0,0),col=1,legend=param_sim,lwd=rep(1.5,n),cex=1.2)
      
      
      if (new.window==TRUE) {dev.new()}
      #par(mar=c(5,5,5,5))
      plot(x=1:npt,y=1:npt,xlab="",ylab="",sub="RNV legend",main=paste("Simulation ", s),xlim=c(0,1),ylim=c(0,1),type="n")
      #legend("center",legend=c("Sup bounds","Inf bounds for small RNV","Inf bounds for big RNV"),lwd=rep(1.5,3),lty=c(3,2,4))
      legend("center",legend=c("Sup bounds","Inf bounds for small RNV","Inf bounds for big RNV"),lwd=rep(1.5,3),lty=c(mycol_RNV$sup$lty[1],mycol_RNV$inf$lty),col=c(mycol_RNV$sup$col[1],mycol_RNV$inf$col))
      
    }
    
    
    ## RNV on concentrations
    if (gr.RNV.E==TRUE) {
      if (new.window==TRUE) {dev.new()}
      #par(mar=c(5,5,5,5))
      plot(seq(1,ngene,by=pasobs),rep(1,npt),type="n",xlab="Generation",ylab="Enzyme concentrations with RNV",xlim=c(0,ngene),ylim=c(0,max(tabE[,1:n])),...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
      #for each enzyme
      for (j in 1:n){
        #for each RNV
        for (nor in 1:nb_RNV) { 
          #sup bounds of RNV in light grey
          lines(seq(1,ngene,by=pasobs),RNV_dis_enz_graph[[j]][,2*nor],col=col_RNV[2],lty=lty_RNV[2],...)#col="grey80",lty=3,...)#,lwd=1.5)
          #inf bounds of RNV in dark grey
          lines(seq(1,ngene,by=pasobs),RNV_dis_enz_graph[[j]][,2*nor-1],col=col_RNV[1],lty=lty_RNV[1],...)#col="grey60",lty=(2*(nor)),...)#,lwd=1.5)
        }
        #simul
        lines(seq(1,ngene,by=pasobs),tabE[,j],col=j+1,lty=1,...)#,lwd=1.5)
        #predict concentrations for enzyme j
        if (pmutA==0&add.eq==TRUE) { 
          if (correl_fun=="Comp") {abline(h=eq_th_e[s,j]*Etot_fun,lty=3,col=j+1,lwd=lwd.eq)}
          abline(h=eq_eff_E[s,j],lty=3,col=j+1,lwd=lwd.eq)
        }
      }
      #legend("topright",legend=c("Sup bounds","Inf bounds for small RNV","Inf bounds for big RNV"),lwd=rep(1.5,3),lty=c(3,2,4))
    }
    
    
    ## RNV size
    if (gr.RNV.size==TRUE) {
      
      if (new.window==TRUE) {dev.new()}
      #par(mar=c(5,5,5,5))
      plot(seq(1,ngene,by=pasobs),rep(1,npt),type="n",xlab="Generation",ylab="RNV size",xlim=c(0,ngene),ylim=c(0,max(unlist(RNV_dis_size),na.rm = TRUE)),...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
      for (j in 1:n) {
        for (nor in 1:nb_RNV) {
          lines(seq(1,ngene,by=pasobs),RNV_dis_size[[nor]][,j],col=j+1,lty=(2*(nor)),...)#,lwd=1.5)
        }
      }
      
      ## RNV size divided by Etot
      if (new.window==TRUE) {dev.new()}
      #par(mar=c(5,5,5,5))
      plot(seq(1,ngene,by=pasobs),rep(1,npt),type="n",xlab="Generation",ylab="RNV size divided by Etot",xlim=c(0,ngene),ylim=c(0,max(unlist(RNV_dis_size_divEtot),na.rm = TRUE)),...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
      for (j in 1:n) {
        for (nor in 1:nb_RNV) {
          lines(seq(1,ngene,by=pasobs),RNV_dis_size_divEtot[[nor]][,j],col=j+1,lty=(2*(nor)),...)#,lwd=1.5)
          abline(h=RNV_proxy_divEtot[nor,j],col=j+1,lty=(2*nor),lwd=lwd.eq)
        }
      }

    }
    
    
    ## Delta
    if (gr.RNV.delta==TRUE) {
      range_delta_graph_enz <- NULL
      
      #to adapt window size of graph
      for (j in 1:n) {
        #take absolute of extreme values of delta for each enzyme
        range_delta_graph_enz <- c(range_delta_graph_enz,abs(range(RNV_dis_delta_graph,na.rm = TRUE)))
      }
      ylim_delta <- c(-min(range_delta_graph_enz),min(range_delta_graph_enz))
      
      if (new.window==TRUE) {dev.new()}
      #par(mar=c(5,5,5,5))
      plot(seq(1,ngene,by=pasobs),rep(1,npt),type="n",xlab="Generation",ylab="Delta",main="RNV",xlim=c(0,ngene),ylim=ylim_delta,...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
      for (j in 1:n) {
        for (nor in 1:nb_RNV) {
          #sup bounds
          lines(seq(1,ngene,by=pasobs),RNV_dis_delta_graph[[j]][,2*nor],col=j+1,lty=3,...)#,lwd=1.5)
          #inf bounds
          lines(seq(1,ngene,by=pasobs),RNV_dis_delta_graph[[j]][,2*nor-1],col=j+1,lty=(2*(nor)),...)#,lwd=1.5)
        }
      }
      #legend("topright",legend=c("Sup bounds","Inf bounds for small RNV","Inf bounds for big RNV"),lwd=rep(1.5,3),lty=c(3,2,4))
    }


    
    
    
    if (gr.RNV.J==TRUE) {
      
      #maximal flux
      max.J <- max(res_sim[,2*n+3])
      
      #if no demand for zoom on flux, by default
      if (length(choice.zoom.J)==0) {
        #by default, zoom on last half of simul 
        zoom.RNV.J <- range(RNV_dis_flux[round(npt/2):npt,])
        flag.zoom.E <- "default"
      } else {
        #if no want of zoom
        if (all(choice.zoom.J==FALSE)) {
          flag.zoom.E <- FALSE
          zoom.RNV.J <- c(0,max.J)
        } else {
          #if particular demand for zoom on flux but not correctly input
          if (length(choice.zoom.J)!=2) {
            #no zoom and warning message
            flag.zoom.E <- FALSE
            zoom.RNV.J <- c(0,max.J)
            if (s==which.sim[1]) {
              #warning message only for first viewed simul
              warning("Zoom on graph flux with RNV is not accurate. Please input a vector of length 2 for 'zoom.RNV.J' of 'NULL' for default.")
            }
          } else {
            #if particular demand for zoom on flux with correct input
            zoom.RNV.J <- choice.zoom.J
            flag.zoom.E <- "choice"
          }
        }
        
      }
      
      
      
      ## J by time with RNV
      if (new.window==TRUE) {dev.new()}
      #par(mar=c(5,5,5,5))
      plot(seq(1,ngene,by=pasobs),res_sim[,(2*n+3)],"l",lty=1,xlab="Generation",xlim=c(0,ngene+10),ylim=zoom.RNV.J,ylab="Flux",...)#,cex.main=2,cex.lab=1.6,cex.axis=1.3)
      #sup bounds of RNV in light grey
      lines(seq(1,ngene,by=pasobs),RNV_dis_flux[,2],col=col_RNV[2],lty=lty_RNV[2],...)#col="grey80",lty="longdash",...)#,lwd=1.5) lty=5
      #inf bounds of RNV in dark grey
      lines(seq(1,ngene,by=pasobs),RNV_dis_flux[,1],col=col_RNV[1],lty=lty_RNV[1],...)#,lwd=1.5) lty=2
      if (pmutA==0&add.eq==TRUE) {
        abline(h=eq_eff_J[s],lty=3,col=1,lwd=lwd.eq)
      }
      
      
      ## J in relation to enzymes
      #zoom.RNV.E <- c(0,max(max(tabE[,1:n]),range(RNV_dis_enz_graph,na.rm=TRUE)[2]))
      zoom.RNV.E <- c(0,max(tabE[,1:n]))
      
      #take only E1 to En and RNV for flux (inf then sup bounds) and flux in simul (col n+3)
      RNV_enz_flux <- cbind(res_sim[,1:n],RNV_dis_flux,res_sim[,(2*n+3)])
      
      #one graph by enzyme
      for (j in 1:n){
        #max.E <- max(max(tabE[round(npt/2):npt,j]),range(RNV_dis_enz_graph[[j]][round(npt/2):npt,],na.rm=TRUE)[2])
        if (flag.zoom.E=="default") {
          #by default, zoom on last half of simul for current enzyme
          #zoom.RNV.E <- range(RNV_dis_enz_graph[[j]][round(npt/2):npt,],na.rm=TRUE)
          zoom.RNV.E <- range(tabE[round(npt/2):npt,j],na.rm=TRUE)
        }
        if (flag.zoom.E=="choice") {
          #which rows as J is between chosen min and chosen max for zoom on J for current simul s
          choice.J <- which(res_sim[,2*n+3]>=min(zoom.RNV.J)&res_sim[,2*n+3]<=max(zoom.RNV.J))
          #take E_j on these rows for current enzyme j
          select_E <- res_sim[choice.J,j]
          zoom.RNV.E <- range(select_E)
        }
        
        #ordered depending on enzyme j
        RNV_enz_flux_ordered <- RNV_enz_flux[order(res_sim[,j]),]
        
        #J in relation to enzymes, with RNV and ordered enzymes
        if (new.window==TRUE) {dev.new()}
        #par(mar=c(5,5,5,5))
        plot(x=1:npt,y=1:npt,xlab=paste("Concentration of enzyme ",j,sep=""),ylab="Flux",xlim=zoom.RNV.E,ylim=zoom.RNV.J,type="n",...)#)
        #inf bounds of RNV in dark grey
        lines(RNV_enz_flux_ordered[,j],RNV_enz_flux_ordered[,n+1],col=col_RNV[1],lty=lty_RNV[1],...)#col="grey60",lty="dashed",...)#,lwd=1.5)
        #sup bounds of RNV in light grey
        lines(RNV_enz_flux_ordered[,j],RNV_enz_flux_ordered[,n+2],col=col_RNV[2],lty=lty_RNV[2],...)#col="grey80",lty="longdash",...)#,lwd=1.5)
        #simul
        lines(RNV_enz_flux_ordered[,j],RNV_enz_flux_ordered[,n+3],col=j+1,lty=1,...)#,lwd=1.5)
        #predict concentrations for enzyme j
        if (pmutA==0&add.eq==TRUE) { 
          abline(v=eq_eff_E[s,j],lty=3,col=j+1,lwd=lwd.eq)
        }
        #max flux
        abline(h=max(res_sim[,(2*n+3)]))#,col="grey80")
      }
      
      #J in relation to enzyme, but flux fit curve
      # if (new.window==TRUE) {dev.new()}
      # #par(mar=c(5,5,5,5))
      # plot(x=1:npt,y=1:npt,xlab=paste("Concentration of enzyme ",j,sep=""),ylab="Flux",xlim=zoom.RNV.E,ylim=zoom.RNV.J,type="n")
      # #for each RNV
      # for (nor in 1:nb_RNV) { 
      #   #sup bounds of RNV in light grey
      #   lines(RNV_dis_enz_graph[[j]][,2*nor],RNV_dis_flux[,2],col="grey80",lty=3,lwd=1.5)
      #   #inf bounds of RNV in dark grey
      #   lines(RNV_dis_enz_graph[[j]][,2*nor-1],RNV_dis_flux[,1],col="grey60",lty=(2*(nor)),lwd=1.5)
      # }
      # #simul
      # lines(tabE[,j],res_sim[,(2*n+3)],col=j+1,lty=1,lwd=1.5)
      # #predict concentrations for enzyme j
      # if (pmutA==0&add.eq==TRUE) { 
      #   abline(v=eq_eff_E[s,j],lty=3,col=j+1,lwd=lwd.eq)
      # }
      
    }
    
    
    #End for simulation s
  }
  
  return(invisible(list(RNV_all_sim=RNV_all_sim,RNV_mean_size=RNV_mean_size)))
}







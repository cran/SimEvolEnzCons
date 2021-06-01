#' @name graph.simul.triangle.diagram.e
#' @aliases graph.simul.triangle.diagram.e
#' 
#' @title Triangular diagram of relative concentrations for simulations of enzyme evolution
#'
#' @description Graphics of enzyme evolution simulations obtained by function \code{\link{simul.evol.enz.multiple}}.
#' Function \code{graph.simul.triangl.diagram.e} gives a triangular diagram of relative concentrations.
#' 
#' 
#' 
#' @details 
#' \bold{WARNING! If there is more than three enzymes in simulations (\eqn{n>3}), be careful for interpretations.}
#' Indeed, triangular diagram is projection of a plane for which sum of coordinates is equal to 1. 
#' In other words, if \eqn{n>3}, represented points are not strictly relative concentrations such as \eqn{e_i = E_i / Etot},
#' but are \eqn{x_i = E_i / (E_i + E_j + E_k)} for the three represented enzymes \code{i}, \code{j} and \code{k}.
#' 
#' Triangular diagram cannot be developed if there is less than three enzymes (\eqn{n<3} or \code{length(which.enz)<3}).
#' Indeed, it is a three-dimensional graph. 
#' 
#' If there is only one simulation (length of \code{which.sim} equal to 1), color is black. Else, color are taken in palette \code{rainbow}.
#' 
#' 
#' Option \code{add.line.eq.eff=TRUE} adds a line of effective equilibrium for all possible initial concentrations.
#' This line is only available for certain constraints, which are \code{"CRPos"} and \code{"CRNeg"}.
#' Code is written only for three enzymes \code{n=3}, and requires to use same initial kinetic parameters for all simulations (set \code{same.kin0=TRUE}) and no mutation of activities (set \code{pmutA=0}).
#' Also, same initial total concentration for all enzyme is required. Easiest way is to set \code{is.random.E0=TRUE}.
#' See function \code{\link{predict_eff_allE0}}.
#' 
#' Option \code{add.curve.lines=TRUE} adds contour lines of flux dome.
#' Contour lines are only available in case of competition, which are \code{"Comp"}, \code{"CRPos"} and \code{"CRNeg"}.
#' It is only available for three enzymes \code{n=3}, and other conditions are the same than option \code{add.line.eq.eff=TRUE}.
#' \code{nbniv} and \code{niv.palette} adjust number and color of contour lines respectively.
#' See function \code{\link{flux.dome.projections}}.
#' 
#' 
#' In triangular diagram, left axis correspond to first value of \code{which.enz} for enzyme number, bottom axis to its second value and right axis to its third value.
#' If options \code{add.line.eq.eff=TRUE} or \code{add.curve.lines=TRUE} are chosen and if all tests have been successfully passed, \code{which.enz} is automatically set to \code{c(1,2,3)}. 
#' 
#' 
#' 
#' @inheritParams simul.evol.graph.methods
#' @param which.enz Numeric. Which enzymes would be represented? Default is the first three enzymes, i.e. \code{c(1,2,3)}.
#' @param add.line.eq.eff Logical. Add line of effective equilibrium for all initial concentrations? Default is \code{FALSE}. \emph{See details.}
#' @param add.curve.lines Logical. Add curve lines for flux dome? Default is \code{FALSE}. \emph{See details.}
#' @inheritParams flux.dome.projections
#' @param cex.ini,cex.th,cex.eff Numeric. Size of remarkable points (respectively initial, theoretical and effective relative concentrations).
#'
#' @import graphics
#' @importFrom grDevices rainbow
#' @importFrom ade4 triangle.plot
#'
#'  
#' @seealso 
#' See \code{\link{simul.evol.graph.methods}} for other plots of enzyme evolution simulation.
#' 
#' See \code{\link{flux.dome.projections}} to have details about contour lines.
#'
#' @return Invisible matrix \code{why_no_lines} who explain why line of effective equilibrium or curve lines not appear even the corresponding options are set to \code{TRUE}.
#' If \code{why_no_lines} is \code{NULL}, options \code{add.line.eq.eff} and \code{add.curve.lines} have been set to \code{FALSE} (default).
#'
#' @examples 
#' 
#'  \donttest{
#'  #Saved simulation
#'  data(data_sim_RegNeg)
#'  graph.simul.triangle.diagram.e(data_sim_RegNeg, new.window=TRUE, add.line.eq.eff=TRUE)
#'  
#'  
#'  #add curve lines
#'  data(data_sim_Comp)
#'  graph.simul.triangle.diagram.e(data_sim_Comp,new.window=TRUE,add.curve.lines=TRUE)
#'  
#'  #all options
#'  data(data_sim_CRNeg)
#'  graph.simul.triangle.diagram.e(data_sim_CRNeg,new.window=TRUE,add.curve.lines=TRUE,
#'  add.line.eq.eff=TRUE)
#'  }
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
#' #beta_12 = 0.32, beta_13 = beta_12 x beta_23, beta_23 = -0.43
#' t_beta <- t(beta) #because R fills matrix column by column
#' beta[lower.tri(beta)] <- 1/t_beta[lower.tri(t_beta)] #beta_ji = 1/beta_ij
#' if (n==3) {beta[lower.tri(beta)] <- 1/beta[upper.tri(beta)]} #only available if n=3
#' correl <- "RegNeg"
#' 
#' evol_sim <- simul.evol.enz.multiple(E0,kin,Keq,nsim,N,correl,beta,npt=250,
#' is.random.E0=TRUE,same.E0=FALSE)
#' graph.simul.triangle.diagram.e(evol_sim, new.window=TRUE, add.line.eq.eff=TRUE)
#' 
#' #add curve lines 
#' nsim <- 3
#' correl <- "Comp"
#' evol_sim <- simul.evol.enz.multiple(E0,kin,Keq,nsim,N,correl,beta,npt=250,
#' is.random.E0=TRUE,same.E0=FALSE)
#' graph.simul.triangle.diagram.e(evol_sim,new.window=TRUE,add.curve.lines=TRUE)
#'
#' #all options
#' correl <- "CRNeg"
#' evol_sim <- simul.evol.enz.multiple(E0,kin,Keq,nsim,N,correl,beta,npt=250,
#' is.random.E0=TRUE,same.E0=FALSE)
#' graph.simul.triangle.diagram.e(evol_sim,new.window=TRUE,add.curve.lines=TRUE,add.line.eq.eff=TRUE)
#' 
#' #several enzyme
#' n <- 5
#' E0 <- c(30,30,30,30,30)
#' kin <- c(1,10,30,100,1000)
#' Keq <- c(1,1,1,1,1)
#' correl <- "SC"
#' evol_sim <- simul.evol.enz.multiple(E0,kin,Keq,nsim,N,correl,beta,npt=250)
#' graph.simul.triangle.diagram.e(evol_sim,new.window=TRUE)
#' }
#'
#' 
NULL











# Triangular diagram
#################################
#################################
#' @rdname graph.simul.triangle.diagram.e
#' 
#' @usage graph.simul.triangle.diagram.e(all_res_sim,which.enz=c(1,2,3),which.sim=NULL,
#' new.window=FALSE,add.eq=TRUE,add.line.eq.eff=FALSE,add.curve.lines=FALSE,
#' nbniv=9,niv.palette=NULL,posi.legend="topleft",cex.ini=1,cex.th=1.2,cex.eff=0.6)
#' 
#' @inheritParams flux.dome.projections
#' 
#' 
#' @return Function \code{graph.simul.triangle.diagram.e} returns point coordinates on triangular diagram.
#' @export




graph.simul.triangle.diagram.e <- function(all_res_sim,which.enz=c(1,2,3),which.sim=NULL,new.window=FALSE,add.eq=TRUE,
                                       add.line.eq.eff=FALSE,add.curve.lines=FALSE,nbniv=9,niv.palette=NULL,posi.legend="topleft",
                                       cex.ini=1,cex.th=1.2,cex.eff=0.6) {
  #here, variable 'i' is always used to indicate simulation number and 'j' for enzyme number
  
  
  ########## Test
  if (length(which.enz)<3) {
    stop("Triangular diagram need three dimensions: 'which.enz' should be a vector of length 3.")
  }
  if (length(which.enz)>3) {
    warning("Triangular diagram need three dimensions: only three first enzymes in 'which.enz' are taken.")
  }
  
  
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
  #tabP_c <- all_res_sim$tabP_r
  
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
      #eq_eff_E <- rbind(eq_eff_E,as.vector(eq_eff_by_sim$pred_E))
      #eq_eff_tau <- rbind(eq_eff_tau,as.vector(eq_eff_by_sim$pred_tau))
    }
  } else {
    #effective equilibrium does not exist
    eq_eff_e <- matrix(NA,nrow=nsim,ncol=(2*n))
    #eq_eff_E <- matrix(NA,nrow=nsim,ncol=(2*n))
  }
  
  
  ####### Other tests
  if (max(which.enz)>n|min(which.enz)<1) {
    stop("Selected enzymes are out of bounds. 'which.enz' need to include values between 1 and n.")
  }
  
  ############ Test conditions for options add.line.eq.eff=TRUE and add.curve.lines==TRUE
  #if a condition is not respected
  cond.line.eq.eff <- FALSE
  cond.curve.lines <- FALSE
  
  #explain why condition is not respected
  why_no_lines <- NULL
  
  # Equilibrium for all E0
  if (add.line.eq.eff==TRUE|add.curve.lines==TRUE) {
    #recall use of special options
    if (add.line.eq.eff==TRUE) {why_no_lines <- rbind(why_no_lines, ("Option 'add.line.eq.eff' is set to TRUE.\n") )}
    if (add.curve.lines==TRUE) {why_no_lines <- rbind(why_no_lines, ("Option 'add.curve.lines' is set to TRUE.\n") )}
    why_no_lines <- rbind(why_no_lines, "Several tests will be passed. See details.\n")
    
    #Total number of enzymes authorized
    if (n==3) {
      
      #if activities change during simulations or between simulations
      if (pmutA==0&same.kin0==TRUE) {
        #keep A0, which equal for all simulations
        A_fun <- leg_A[1,1:n]
        
        #compute Etot0 for all simulation
        all_sum_E0 <- apply(leg_E[,1:n],1,sum)
        #keep first Etot0
        Etot_0 <- all_sum_E0[1]
        
        #test if all Etot0 are equal to same value
        if (isTRUE(all.equal(rep(Etot_0,nsim), all_sum_E0))) {
          
          #possible constraints for effective equilibrium only
          if (add.line.eq.eff==TRUE) {
            if (correl_fun=="RegNeg"|correl_fun=="CRPos"|correl_fun=="CRNeg") {
              #memory of successful tests for add.line.eq.eff
              cond.line.eq.eff <- TRUE
              which.enz <- c(1,2,3)
              why_no_lines <- rbind(why_no_lines, "Tests have been successfully passed for option 'add.line.eq.eff=TRUE'.\n")
            } else { #not good correl_fun==
              why_no_lines <- rbind(why_no_lines, "Not convenient constraint 'correl' for effective equilibrium.\n")
            }
          }
          
          #possible constraints for flux dome only
          if (add.curve.lines==TRUE) {
            if (correl_fun=="Comp"|correl_fun=="CRPos"|correl_fun=="CRNeg") {
              #memory of successful tests for add.curve.lines
              cond.curve.lines <- TRUE
              which.enz <- c(1,2,3)
              why_no_lines <- rbind(why_no_lines, "Tests have been successfully passed for option 'add.curve.lines=TRUE'.\n")
            } else { #not good correl_fun==
              why_no_lines <- rbind(why_no_lines, "Not convenient constraint 'correl' for flux dome. Need competition.\n")
            }
          }
          
        } else { #Etot_0 are different between simulation
          why_no_lines <- rbind(why_no_lines, "Initial total concentration need to be equal between simulations for these options. Simpliest way is to set 'is.random.E0' to TRUE.\n")
        }
        
      } else { #not pmutA==0 or not same.kin0==TRUE
        why_no_lines <- rbind(why_no_lines, "Activities are changed during or between simulations. Please set 'pmutA' to 0 and 'same.kin0' to TRUE for these options.\n")
      }
      
    } else { #not n==3
      why_no_lines <- rbind(why_no_lines, "Total number of enzymes 'n' is limited to three for this option.\n")
    }
    
    #explain was passiert if tests not passed
    if (add.line.eq.eff==TRUE & cond.line.eq.eff==FALSE) {
      why_no_lines <- rbind(why_no_lines, "Tests have not been passed. Line of effective equilibrium will not appear.\n")
    }
    if (add.curve.lines==TRUE & cond.curve.lines==FALSE) {
      why_no_lines <- rbind(why_no_lines, "Tests have not been passed. Curve lines will not appear.\n")
    }
    
    #end of conditions for option add.line.eq.eff=TRUE and add.curve.lines=TRUE
  }
  
  
  
  #########################
  #######""" Graphics
  #############################
  
  #which simulations will be represented. If NULL, all
  if (length(which.sim)==0) {which.sim <- 1:nsim}
  
  #color vector for simulations
  mycol_sim <- rainbow(nsim)
  #if only one simulation, color is black
  if (length(which.sim)==1) {mycol_sim <- 1}
  #color for theoretical equilibrium (identical for all simulation if A does not change or Reg constraint)
  if (pmutA==0) {
    #col_th <- rep(1,nsim)
    #col_th <- rep('grey60',nsim)
    col_th <- rep('white',nsim)
  } else {col_th <- mycol_sim}
  if (correl_fun!="SC"&correl_fun!="Comp") {col_th <- rep("white",nsim)} #1/B does not change
  #color for effective equilibrium: its changes if different A0 or E0
  #if (same.kin0==TRUE&same.E0==TRUE&pmutA==0) {col_eff <- rep(1,nsim)} else {col_eff <- mycol_sim}
  if (pmutA==0) {col_eff <- rep(1,nsim)} else {col_eff <- mycol_sim}
  
  
  ###########" Triangular diagram for 3 enzymes i, j and k on plane (ei + ej + ek = 1)
  
  #### Use of triangle.plot
  #Need a positive matrix of 3 columns
  #If data.frame, use colnames as axis label
  #Each use create a new graph over precedent
  #cpoint=0 => erase points
  
  
  #### Relative concentrations of chosen enzymes
  tab_e <- tabR[,which.enz]/tabR[,(2*n+1)]
  # Take simulations number
  tab_e <- cbind.data.frame(tab_e,tabR$sim)
  # Name columns to add correct title on graph
  colnames(tab_e) <- c(paste("e",which.enz[1],sep=""), paste("e",which.enz[2],sep=""), paste("e",which.enz[3],sep=""), "sim")
  
  
  ##### Special points
  if (new.window==TRUE) {dev.new()}
  #par(mar=c(5,5,5,5))
  
  #Equilibrium
  if (add.eq==TRUE) {
    # theoretical equilibrium
    #in case of equilibrium out of bounds (negative value not available for triangle.plot)
    if (correl_fun=="RegNeg"|correl_fun=="CRNeg") {
      #eq_th_e[i,which.enz] = 1/B in these cases
      B_select <- B_fun[which.enz]
      #position in x for plot
      B_x <- (1/B_select[2] - 1/B_select[1])/sqrt(2)
      #B_x <- (1/B[which.enz[2]] - 1/B[which.enz[1]])/sqrt(2)
      #B_x <- (1/B[2] - 1/B[1])/sqrt(2)
      #position in y for plot
      B_y <- (2 * 1/B_select[3] - 1/B_select[2] - 1/B_select[1])/sqrt(6)
      #B_y <- (2 * 1/B[which.enz[3]] - 1/B[which.enz[2]] - 1/B[which.enz[1]])/sqrt(6)
      #B_y <- (2 * 1/B[3] - 1/B[2] - 1/B[1])/sqrt(6)
      #points(B_x,B_y,col=1,pch=8,cex=1.5,lwd=1.3) #darkgrey si CRPos
      #save nsim times these values, to not copy same test several times further
      tr_th <- matrix(c(B_x,B_y),ncol=2,nrow=nsim,byrow=TRUE)
    } else {
      #in other cases
      eq_th_select <- eq_th_e[,which.enz]
      #because R transforms matrix in numeric if there is only one row, thus eq_th_select dont have the right format
      #all() because R gives two classes for matrix rather one
      if (all(class(eq_th_select)=="numeric")) {eq_th_select <- t(as.data.frame(eq_th_select))}
      tr_th <- triangle.plot(eq_th_select,show.position = FALSE,scale = FALSE,draw.line = FALSE)
      #tr_th <- triangle.plot(t(as.data.frame(1/B)),show.position = FALSE,scale = FALSE,draw.line = FALSE)
    }
    
    # effective equilibrium
    if (correl_fun=="RegNeg"|correl_fun=="CRPos"|correl_fun=="CRNeg") {
      #NA not authorized for triangle.plot
      eq_eff_select <- eq_eff_e[,which.enz]
      tr_eff <- triangle.plot(eq_eff_select,show.position = FALSE,scale = FALSE,draw.line = FALSE)
    } else {
      tr_eff <- matrix(NA,ncol=2,nrow=nsim)
    }
  }
  
  # Equilibrium for all E0
  #conditions -> see tests up
  if (cond.line.eq.eff==TRUE) {
    #compute effective equilibrium for all E0 if Etot0 and A0 is equal for all simulations
    eq_eff_allE0 <- predict_eff_allE0(B_fun,A_fun,correl_fun,Etot_0,X_fun)
    #take relative concentrations only
    e_Jmax <- eq_eff_allE0$all_eq_eff[,1:n]
    #coordinates for triangle.plot
    tr_Jmax <- triangle.plot(e_Jmax[,which.enz[1:3]],show.position = FALSE,scale = FALSE,draw.line = FALSE)
  }
  
  
  
  ##### Triangular diagram 'skeletal'
  
  #new window
  if (new.window==TRUE) {dev.new()}
  #par(mar=c(4,4,4,4))
  
  #triangular gaph and create points for first three enzymes of selected enzymes
  #selected enzymes are already taken in good order in tab_e
  tr_x <- triangle.plot(tab_e[,1:3],scale=FALSE,show.position=FALSE,draw.line = FALSE,cpoint=0) #diagramme triangulaire de e1,e2,e3 pour toutes les simulations (convertis les points en coordonnÃ©es xy)
  
  #if conditions for curve lines are completed
  if (cond.curve.lines==TRUE) {
    flux.dome.projections(A_fun,Etot_0,add.reg=FALSE,X_fun=X_fun,nbniv=nbniv,niv.palette=niv.palette,new.window=new.window,posi.legend = posi.legend)
    #because add.reg=FALSE, projection on section tau is not created
    
    #add line around the triangle plot, because curve line is added on the triangle line
    A <- c(-1/sqrt(2), -1/sqrt(6))
    B <- c(1/sqrt(2), -1/sqrt(6))
    C <- c(0, 2/sqrt(6))
    polygon(c(A[1], B[1], C[1]), c(A[2], B[2], C[2]))
  }
  
  
  #### Add simulation points
  #for each chosen simulation
  for (i in which.sim){
    #add all points for this simulation, even first one
    points(tr_x[((i-1)*npt+2):(i*npt),],col=mycol_sim[i],lwd=0.5,cex=1,pch=16)
    
    #initial point for this simulation = triangle point (pch=17) / square (pch=15) / losange with background (pch=23)
    #revert x and y, because its a vector of length 2
    points(t(tr_x[((i-1)*npt+1),]),col=1,lwd=0.5,cex=cex.ini,pch=22,bg=mycol_sim[i])
    #points(t(tr_x[((i-1)*npt+1),]),col=mycol_sim[i],lwd=0.5,cex=0.7,pch=15)
    #points(tr_x[((i-1)*npt+1),1],tr[((i-1)*npt+1),2],col=mycol_sim[i],lwd=0.5,cex=0.7,pch=17)
  }
  
  
  #### Add other points
  #Line of effective equilibrium
  if (cond.line.eq.eff==TRUE) {
    points(tr_Jmax[,1],tr_Jmax[,2],col=1,cex=0.3,pch=16)
  }
  
  if (add.eq==TRUE) {
    
    for (i in which.sim) {
      #theoretical equilibrium = star
      points(tr_th[i,1],tr_th[i,2],col=1,pch=23,cex=cex.th,lwd=1,bg=col_th[i]) #filled square
      #points(tr_th[i,1],tr_th[i,2],col=col_th[i],pch=5,cex=1.5,lwd=1) #darkgrey si CRPos
      #points(tr_th[i,2],tr_th[i,1],col=col_th[i],pch=8,cex=1.5,lwd=1.3) #darkgrey si CRPos
      
      # then effective equilibrium = cross-point
      points(tr_eff[i,1],tr_eff[i,2],col=1,cex=cex.eff,lwd=0.5,pch=21,bg=col_eff[i])
      #points(tr_eff[i,1],tr_eff[i,2],col=col_eff[i],cex=1,lwd=1.8,pch=1) #cross-point = 16
    }
    
  }
  
  return(invisible(why_no_lines))
}
  
  
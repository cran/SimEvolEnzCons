#' Plots of flux shape computing from one point
#'
#' Graphics for illustrating results of function \code{\link{flux.shape.from.one.point}}
#' 
#'  
#' 
#'
#' @details 
#' See function \code{\link{flux.shape.from.one.point}}.
#' 
#' Gives graphics of:
#' \enumerate{
#'    \item Selection coefficient (discrete expression) in relation to apparent mutation effect \eqn{\delta}, which is \eqn{E^m-E^r} (option \code{gr.s.delta})
#'    \item Flux in relation to mutant concentration \eqn{E^m} (option \code{gr.J.E})
#'    \item Selection coefficient (continuous expression) in relation to mutant concentration (option \code{gr.s_cont.E})
#'    \item Selection coefficient (discrete expression) in relation to mutant concentration (option \code{gr.s_disc.E})
#'    \item Comparison between the two expressions of selection coefficient (option \code{gr.comp.s})
#'    \item If driving variable \eqn{\tau} exists: flux/mutant concentrations/discrete selection coefficient in relation to tau (option \code{gr.tau})
#' }
#' 
#' Each enzyme is represented by one color. If there are three enzymes enzymes or less, enzyme 1 is in red, 2 in green and 3 in blue.
#' For more than three enzymes, colors are taken in palette \code{"Spectral"} from package '\code{RColorBrewer}'.
#' 
#' 
#' \emph{Discrete expression} of selection coefficient refers to \eqn{s=(J^m-J^r)/J^r}. Available whatever values of enzymes concentrations. See function \code{\link{coef_sel.discrete}}.
#' \emph{Continuous expression} of selection coefficient refers to \eqn{s_i = R_Ei^J delta_i/E_i}. Available only for small values of delta_i. See function \code{\link{coef_sel.continue}}.
#' 
#' 
#' The RNV corresponds to enzyme concentrations such as selection coefficient is between -1/N and 1/N.
#'   
#'  \bold{Zooms}
#'   
#'  If \code{NULL}, default values are applied.
#'  
#'  If vector of length 1, symmetry is applied (i.e. -zoom,zoom).
#'  
#'  If vector of length 2, zoom values are applied directly.
#'  
#'  \bold{posi.legend}
#'  
#'  Indicates coordinates of the upper left corner of the legend.
#'  Available only for graph corresponding to \code{gr.s.delta}.
#'  See more in details of function \code{\link{flux.dome.projections}}.
#' 
#' 
#' @importFrom RColorBrewer brewer.pal		
#' 
#' @usage flux.shape.from.one.point.graphics(list_fsfop, N_fun, zoom_delta=NULL,add.eq=TRUE,
#' gr.s.delta=TRUE,gr.J.E=FALSE,gr.s_cont.E=FALSE,gr.s_disc.E=FALSE,
#' gr.comp.s=FALSE,gr.tau=FALSE,posi.legend="topleft",...)
#'     
#' @param list_fsfop Output of the function \code{\link{flux.shape.from.one.point}}
#' @param N_fun Numeric. Population size that influence neutral zone
#' @param zoom_delta Zoom on axis of apparent mutation effects. \emph{See details.}
#' @inheritParams simul.evol.graph.methods
#' @param add.eq Logical. Does equilibrium appear on graphics?
#' @param gr.s.delta Logical. Add graph of selection coefficient in relation to mutation effect?
#' @param gr.J.E Logical. Add graph of flux J in relation to enzyme concentrations?
#' @param gr.s_cont.E,gr.s_disc.E Logical. Add graph of selection coefficient (respectively from continuous or discrete expression) in relation to enzyme concentrations?
#' @param gr.comp.s Logical. Add graph of comparison of selection coefficient from different expressions (continuous vs discrete)?
#' @param gr.tau Logical. Add graph depending on tau?
#' @param posi.legend Legend position for activities. See \emph{details}.
#'
#'
#' @return Nothing
#' 
#' @seealso 
#' See function \code{\link{RNV.ranking.order.factor}} to see against which parameters the RNV depend.
#'   
#' @examples
#' \donttest{
#' fsfop <- flux.shape.from.one.point(100,c(1,10,30),"SC")
#' flux.shape.from.one.point.graphics(fsfop,1000)
#' }
#'
#' @export









flux.shape.from.one.point.graphics <- function(list_fsfop, N_fun, zoom_delta=NULL,add.eq=TRUE,
                                               gr.s.delta=TRUE,gr.J.E=FALSE,gr.s_cont.E=FALSE,gr.s_disc.E=FALSE,gr.comp.s=FALSE,gr.tau=FALSE,
                                               posi.legend="topleft",...) {

  ######" Settings
  x <- list_fsfop$x
  J_all <- list_fsfop$J
  sel_all_tru <- list_fsfop$sel_disc
  sel_all  <- list_fsfop$sel_cont
  tau_all  <- list_fsfop$tau
  
  E_fun  <- list_fsfop$param$E
  A_fun <- list_fsfop$param$A
  correl_fun <- list_fsfop$param$correl
  beta_fun <- list_fsfop$param$beta
  B_fun <- list_fsfop$param$B
  E_ini_fun <- list_fsfop$param$E0
  Etot_fun <- list_fsfop$param$Etot
  n_fun <- list_fsfop$param$n
  X_fun <- list_fsfop$param$X

  J_fun <- flux(E_fun,A_fun,X_fun)



  ### Equilibrium
  eq_th <- predict_th(A_fun,correl_fun,B_fun)
  if (correl_fun=="RegNeg"|correl_fun=="CRPos"|correl_fun=="CRNeg") {
    eq_eff <- predict_eff(E_ini_fun,B_fun,A_fun,correl_fun)
  } else {
    #because eq eff does not exist in other cases
    eq_eff <- list(pred_e=rep(NA,n_fun),pred_E=rep(NA,n_fun),pred_tau=NA)
  }
  #delta such as Ei*=Eir+delta_i*
  #note that x = E + delta
  if (correl_fun=="SC"|correl_fun=="RegPos"|correl_fun=="RegNeg") {
    #because Etot depends also on delta, 
    delta_eq_th <- (eq_th$pred_e*sum(E_fun)-E_fun)/(1-eq_th$pred_e*B_fun)
  }
  if (correl_fun=="Comp"|correl_fun=="CRPos"|correl_fun=="CRNeg") {
    delta_eq_th <- eq_th$pred_e*sum(E_fun)-E_fun
  }
  delta_eq_eff <- eq_eff$pred_E-E_fun

  

  ###########################"" Graphics
  
  
  # Color for graphics
  if (n_fun<=3) {
    mypalette <- 2:(n_fun+1)
  } else {
    mypalette <- brewer.pal(n_fun,"Spectral")
  }
  
  
  lim_NZ <- 1/N_fun #bounds of neutral zone
  
  ##### Zoom for graphics
  zoom_sel_delta <- 2*lim_NZ
  zoom_sel_E <- 10*lim_NZ
  zoom_comp <- 0.1*Etot_fun/100 #zoom for comparison between formula of selection coeff
  
  lim_E <- c(0,Etot_fun)
  lim_J <- c(0,max(J_all))
  lim_tau <- c(-1,2) #limits of tau in graph

  
  #zoom around mutation effect
  if (length(zoom_delta)==0) {
    #default value
    if (correl_fun=="SC") {zoom_delta <- 0.5}
    if (correl_fun=="Comp") {zoom_delta <- 3}
    if (correl_fun=="RegPos") {zoom_delta <- 0.2}
    if (correl_fun=="RegNeg") {zoom_delta <- 3.5}
    if (correl_fun=="CRPos") {zoom_delta <- 2}
    if (correl_fun=="CRNeg") {zoom_delta <- 2.5}
    #because neutral zone depending on Etot
    zoom_delta <- zoom_delta*Etot_fun/100
  }
  if (length(zoom_delta==1)) {
    #symmetry if only one value is input
    zoom_delta <- c(-zoom_delta,zoom_delta)
  }
  
  #Selection coefficient (discrete) of i fct of Eim-Eir (valable pour tout delta)
  if (gr.s.delta==TRUE) {
    
    plot(main=correl_fun,E_fun[1],J_fun, type='n', xlab="Mutation effect", ylab="Selection coefficient", ylim=c(-zoom_sel_delta,zoom_sel_delta),xlim=zoom_delta,...) #ylim=range(sel_all[2:nrow(sel_all),])
    
    abline(v=0)
    #limits of neutral zone
    abline(h=-lim_NZ,lty=2,...)
    abline(h=lim_NZ,lty=2,...)
    abline(h=0)
    
    for (i in 1:n_fun) {
      #curve of selection coefficient depending on delta_i
      lines(x-E_fun[i],sel_all_tru[,i], col=mypalette[i],...)
      #lines(E_fun[i]+x,sel_all_tru[,i], lwd=1.5, col=mypalette[i]) #courbe lorsque Ei varie
      
      #equilibrium for Ei
      if (add.eq==TRUE) {
        abline(v=delta_eq_th[i],lty=2,col=mypalette[i],...)
        abline(v=delta_eq_eff[i],lty=3,col=mypalette[i],...)
      }
      
    }
    
    # Legend
    if (!is.null(posi.legend)) {
      # if want a legend
      legend(posi.legend[1],posi.legend[2],legend=A_fun,col=mypalette,bty="n",cex=1,lwd=1.2)
    }
    
    
    #add legend depending on RNV-ranking-order-factor
    # if (correl_fun=="SC"|correl_fun=="Comp") {
    #   legend("topleft",col=mypalette,legend=A_fun,lwd=rep(1.5,n),inset=c(0,0),title="Activities")
    # }
    # if (correl_fun=="RegPos"|correl_fun=="RegNeg") {
    #   legend("topleft",col=mypalette,legend=c(A_fun,round(B_fun,2)),lwd=rep(1.5,n_fun),inset=c(0,0),title="A and B",ncol=2)
    # }
    # if (correl_fun=="CRPos"|correl_fun=="CRNeg") {
    #   dis_Be0 <- abs(round(1/B_fun-E_ini_fun/Etot_fun,2))
    #   legend("topleft",col=mypalette,legend=c(A_fun,dis_Be0),lwd=rep(1.5,n_fun),inset=c(0,0),ncol=2,title="A and Distance between e* and e0")
    # }
    
  }
  #dev.off()



  #J fct of Ei
  if (gr.J.E==TRUE) {
    #x11()
    plot(E_fun[1],J_fun, type='n', xlab="Enzyme concentration", ylab="Flux", ylim=lim_J,xlim=lim_E,...)
    for (i in 1:n_fun) {
      #curve of J depending on Ei
      lines(x,J_all[,i], col=mypalette[i],...)
      #lines(E[i]+x,J_all[,i], lwd=1.5, col=mypalette[i]) #courbe lorsque Ei varie
      #resident point
      abline(v=E_fun[i],col=mypalette[i],...)
      if (add.eq==TRUE) {
        #      if (correl_fun=="RegPos"|correl_fun=="CRPos"|correl_fun_fun=="SC"|correl_fun=="Comp") {
        #abline(v=eq_th$pred_e*Etot_fun,lty=2,col=mypalette[i],...)
        abline(v=E_fun[i]+delta_eq_th[i],lty=2,col=mypalette[i],...)
        #}
        #if (correl_fun=="RegNeg"|correl_fun=="CRPos"|correl_fun=="CRNeg"){
        abline(v=eq_eff$pred_E[i],lty=3,col=mypalette[i],...)
        #}
        
      }
      
    }
    abline(h=J_fun)
    #legend("bottomright",col=c(mypalette,1),legend=c(A,15.9938),lwd=rep(1.5,n+1),inset=c(0,0),title="Activities from E1 to En")
  }
  
  # if(correl_fun=="SC"|correl_fun=="RegPos"|correl_fun=="RegNeg"){
  #   x11()
  #   plot(e[1],Jm, type='n', xlab="Relative concentration", ylab="Flux", lwd=1.5, ylim=c(0,max(J_all)),xlim=c(0,1))
  #   for (i in 1:n) {
  #     lines(x/Etot_all[,i],J_all[,i], lwd=1.5, col=mypalette[i]) #courbe lorsque Ei varie
  #     #lines(E[i]+x,J_all[,i], lwd=1.5, col=mypalette[i]) #courbe lorsque Ei varie
  #     abline(v=e[i],col=mypalette[i]) #Ã©quilibre pour Ei et point de dÃ©part
  #   }
  # }


  #Selection coefficient (continuous) of i fct of Ei (available for low delta)
  if (gr.s_cont.E==TRUE) {
    #x11()
    plot(E_fun[1],J_fun, type='n', xlab="Enzyme concentration", ylab="Continuous selection coefficient", ylim=c(-zoom_sel_E,zoom_sel_E),xlim=lim_E,...) #ylim=range(sel_all[2:nrow(sel_all),])
    
    abline(h=-lim_NZ,lty=2)
    abline(h=lim_NZ,lty=2)
    abline(h=0)
    
    for (i in 1:n_fun) {
      #depending on Ei
      lines(x,sel_all[,i], col=mypalette[i],...)
      #lines(E[i]+x,sel_all[,i], lwd=1.5, col=mypalette[i]) #courbe lorsque Ei varie
      #resident point
      abline(v=E_fun[i],col=mypalette[i])
      #equilibrium
      if (add.eq==TRUE) {
        abline(v=E_fun[i]+delta_eq_th[i],lty=2,col=mypalette[i],...)
        abline(v=eq_eff$pred_E[i],lty=3,col=mypalette[i],...)
      }
      
      # if (correl_fun=="RegNeg"|correl_fun=="CRPos"|correl_fun=="CRNeg"){
      #   abline(v=eq_eff$pred_E[i],lty=2,col=mypalette[i])
      # }
      # if (correl_fun=="RegPos"|correl_fun=="CRPos"|correl_fun_fun=="SC"|correl_fun=="Comp") {
      #   abline(v=eq_th$pred_e*Etot_fun)
      # }
    }
    
    #legend("topright",col=mypalette,legend=A,lwd=rep(1.5,n),inset=c(0,0),title="Activities from E1 to En")
  }
  
  #Selection coefficient (discrete) of i fct of Ei (available for all delta)
  if (gr.s_disc.E==TRUE) {
    #x11()
    plot(E_fun[1],J_fun, type='n', xlab="Enzyme concentration", ylab="Discret selection coefficient", ylim=c(-zoom_sel_E,zoom_sel_E),xlim=lim_E,...) #ylim=range(sel_all[2:nrow(sel_all),])
    
    abline(h=-lim_NZ,lty=2)
    abline(h=lim_NZ,lty=2)
    abline(h=0)
    
    for (i in 1:n_fun) {
      #depending on Ei
      lines(x,sel_all_tru[,i], col=mypalette[i],...)
      #lines(E[i]+x,sel_all_tru[,i], lwd=1.5, col=mypalette[i]) #courbe lorsque Ei varie
      #resident point
      abline(v=E_fun[i],col=mypalette[i])
      #equilibrium
      abline(v=E_fun[i]+delta_eq_th[i],lty=2,col=mypalette[i],...)
      abline(v=eq_eff$pred_E[i],lty=3,col=mypalette[i],...)
      # if (correl_fun=="RegNeg"|correl_fun=="CRPos"|correl_fun=="CRNeg"){
      #   abline(v=eq_eff$pred_E[i],lty=2,col=mypalette[i])
      # }
      # if (correl_fun=="RegPos"|correl_fun=="CRPos"|correl_fun_fun=="SC"|correl_fun=="Comp") {
      #   abline(v=eq_th$pred_e*Etot_fun)
      # }
    }

    #legend("topright",col=mypalette,legend=A,lwd=rep(1.5,n),inset=c(0,0),title="Activities from E1 to En")
  }



  #Comparaison of s=Delta(J)/J discret and s=R*delta/Ei continuous
  if (gr.comp.s==TRUE) {
    #x11()
    plot(x=sel_all_tru[1,1],y=sel_all[1,1],type='n',xlab="Discret selection coefficient", ylab="Continuous selection coefficient",xlim=c(-zoom_comp,zoom_comp),ylim=c(-zoom_comp,zoom_comp),...)
    
    abline(v=0)
    abline(h=0)
    abline(a=0,b=1)
    
    for (i in 1:n_fun) {
      lines(sel_all_tru[,i],sel_all[,i],col=mypalette[i],...)
    }
    
  }



  #Function of tau in case of regulation
  if (correl_fun!="SC"&correl_fun!="Comp"&gr.tau==TRUE) {
    #compute tau at resident point
    tau_r <- droite_tau(E_fun,E_ini_fun,B_fun)

    #Flux en fct de tau
    #x11()
    plot(tau_all[1,1],J_fun, type='n', xlab="Tau", ylab="Flux", ylim=lim_J,xlim=lim_tau,...)
    for (i in 1:n_fun) {
      lines(tau_all[,i],J_all[,i], col=mypalette[i],...) #courbe lorsque Ei varie
    }
    abline(h=J_fun)
    #resident point
    abline(v=tau_r,...)
    if (add.eq==TRUE){
      #initial point
      abline(v=0,lty=4,...)
      #theoretical equilibrium
      abline(v=1,lty=2,...)
      #effective equilibrium
      #if (correl_fun!="RegPos"){
      abline(v=eq_eff$pred_tau,lty=3,...)
    }


    #E en fct de tau
    #x11()
    plot(tau_all[1,1],E_fun[1], type='n', ylab="Enzyme concentration", xlab="Tau", xlim=lim_tau,ylim=lim_E,...)
    for (i in 1:n_fun) {
      #Ei depending on tau
      lines(tau_all[,i],x, col=mypalette[i],...)
      #lines(E[i]+x,J_all[,i], lwd=1.5, col=mypalette[i]) #courbe lorsque Ei varie
      #resident point
      abline(h=E_fun[i],col=mypalette[i],...)
    }
    #resident point
    abline(v=tau_r,...)
    if (add.eq==TRUE){
      #initial point
      abline(v=0,lty=4,...)
      #theoretical equilibrium
      abline(v=1,lty=2,...)
      #effective equilibrium
      #if (correl_fun!="RegPos"){
      abline(v=eq_eff$pred_tau,lty=3,...)
    }

    #Coef sel discret en fct de tau
    #x11()
    plot(tau_all[1,1],J_fun, type='n', xlab="Tau", ylab="Discret selection coefficient", ylim=c(-zoom_sel_E,zoom_sel_E),xlim=lim_tau,...) #ylim=range(sel_all[2:nrow(sel_all),])
    
    #neutral zone
    abline(h=-lim_NZ,lty=2)
    abline(h=lim_NZ,lty=2)
    abline(h=0)
    
    for (i in 1:n_fun) {
      lines(tau_all[,i],sel_all_tru[,i], col=mypalette[i],...) #courbe lorsque Ei varie
    }

    #resident point
    abline(v=tau_r,...)
    if (add.eq==TRUE){
      #initial point
      abline(v=0,lty=4,...)
      #theoretical equilibrium
      abline(v=1,lty=2,...)
      #effective equilibrium
      #if (correl_fun!="RegPos"){
      abline(v=eq_eff$pred_tau,lty=3,...)
    }
  }

  #dev.off()


}


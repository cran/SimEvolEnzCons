#' Graph of e and RNV in relation to A at equilibrium
#'
#' Plot of RNV size and relative concentrations in relation to activities at equilibrium
#' 
#'
#' @details 
#' \code{list_eq$A}, \code{list_eq$RNV_size[1,]} and \code{list_eq$e_eq} are binded in a data.frame of \code{n} rows (one by enzyme) and three columns (A, RNV and e),
#' then ordered according to A.
#' 
#' Title (option \code{add.title=TRUE}) is "Comparing e and RNV at equilibrium".
#' 
#' 
#' @usage RNV.graph.double.at.eq(list_eq,posi.legend="topleft",add.title=TRUE,
#' cex.lab=1,mar.lab=2.5,enz.lab=FALSE,...)
#'
#' @param list_eq Output of function \code{\link{RNV.size.at.equilibr}}:
#' \itemize{
#'    \item \code{list_eq$A}: Numeric vector of activities (length \code{n})
#'    \item \code{list_eq$RNV_size}: Numeric matrix of \code{n} column and 2 rows. Only first row will be used, corresponding to near RNV (small mutations).
#'    \item \code{list_eq$e_eq}: Numeric vector of equilibrium (length \code{n})
#'    }
#' @param posi.legend Character string. Where would you put legend? See parameter \code{x} in function \code{legend}. If \code{NULL}, legend will not appear.
#' @param add.title Logical. Add a title to the plot?
#' @inheritParams RNV.size.at.equilibr
#' @param cex.lab Numeric. Size of axis label.
#' @param mar.lab Numeric. Distance of label from axis. 
#' @param enz.lab Logical. Add enzyme name on graph?
#'
#' @return Invisible data.frame of \code{n} rows (one by enzyme) and three columns ($A, $RNV and $e).
#' 
#' 
#' @seealso 
#' See \code{\link{RNV.size.at.equilibr}} to compute RNV. 
#'
#'
#' @examples
#' list_eq <- RNV.size.at.equilibr(20,"Comp",1000,Etot_0=100,show.plot=FALSE)
#' RNV.graph.double.at.eq(list_eq)
#'
#'
#'
#' @export





# #abscisse
# n_fun <- 3
# N_fun <- 1000
# correl_fun <- "Comp"
# A_fun <- c(1,10,30)
# inv_B <- c(0.5,0.3,0.6)
# B_fun <- sum(inv_B)/inv_B
# beta_fun <- compute.beta.from.B(B_fun)
# Etot_fun <- 100
# E0_fun <- c(30,30,30)
# E0_fun <- Etot_fun*E0_fun/sum(E0_fun)




RNV.graph.double.at.eq <- function(list_eq,posi.legend="topleft",add.title=TRUE,cex.lab=1,mar.lab=2.5,enz.lab=FALSE,...) {
  
  #save par() of user
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  
  #data at equilibrium
  #list_eq <- RNV.size.at.equilibr(n_fun,correl_fun,N_fun,A_fun,B_fun,E0_fun,Etot_0,Etot_eq,A_lim,inv_B_lim,E0_lim)

  #Mise en forme
  name_enz <- NULL
  for (j in 1:list_eq$n) {name_enz[j] <- paste0("E",j)}
  #only RNV of small mutation, even there is one RNV at equilibrium, due to computing method
  compar_at_eq <- cbind.data.frame(A=list_eq$A,RNV=list_eq$RNV_size[1,],e=list_eq$e_eq,name=name_enz)
  compar_at_eq_ord <- compar_at_eq[order(compar_at_eq$A),]
  
  plot_main <- NULL
  if (add.title==TRUE) {
    plot_main <- "Comparing e and RNV at equilibrium"
  }
  
  #Graphique
  #x11()
  #enough size at right
  par(mar=c(5,5,3,5)) 
  
  
  ### first plot (in black)
  #line RNV fct A with square
  plot(compar_at_eq_ord$A, compar_at_eq_ord$RNV, pch=15, axes=F, ylim=c(0,max(compar_at_eq_ord$RNV)), xlab="", ylab="", type="o",col="black",main=plot_main,...)#,lwd=1.5)
  #add left y-axis
  axis(2, ylim=c(0,max(compar_at_eq_ord$RNV)),col="black",...)#,cex.axis=1.3)
  #axis legend
  mtext("RNV size",side=2,line=mar.lab,cex=cex.lab)
  #encadre plot
  box()
  
  ### add name of enzyme above RNV points
  if (enz.lab==TRUE) {
    text(compar_at_eq_ord$A,compar_at_eq_ord$RNV,labels=compar_at_eq_ord$name,pos=3)
  }
  
  ### second plot (in red)
  #authorize superposition of plot
  par(new=T)
  #line e fct A with points
  plot(compar_at_eq_ord$A, compar_at_eq_ord$e, pch=16,  xlab="", ylab="", ylim=c(0,1), axes=F, type="o", col="red",...)#,lwd=1.5)
  #add right y-axis
  axis(4, ylim=c(0,1), col="red",col.axis="red",...)#,cex.axis=1.3)
  mtext("Relative concentration",side=4,col="red",line=mar.lab,lwd=2.5,cex=cex.lab)
  
  ### add x-axis
  axis(1,xlim=range(compar_at_eq_ord$A),...)#,cex.axis=1.3)
  #axis(1,pretty(range(A_fun),10)) #to pretty interval
  mtext("Pseudo-activity A",side=1,col="black",line=mar.lab,lwd=2.5,cex=cex.lab)
  
  ### legend
  if (length(posi.legend)!=0) {
    legend(x=posi.legend,legend=c("e","RNV"),text.col=c("black","red"),pch=c(16,15),col=c("black","red"))
  }
  

    
  
  return(invisible(compar_at_eq))
  #return(invisible(list_eq))
}




# 
# ##Utiliser RNV.size.at.equilibr() pour faire les calculs, en ajoutant A, beta et E0 en param ?
# 
# #abscisse
# n_fun <- 3
# N_fun <- 1000
# correl_fun <- "Comp"
# A_fun <- c(1,10,30)
# inv_B <- c(0.5,0.3,0.6)
# B_fun <- sum(inv_B)/inv_B
# beta_fun <- compute.beta.from.B(B_fun)
# Etot_fun <- 100
# E0_fun <- c(30,30,30)
# E0_fun <- Etot_fun*E0_fun/sum(E0_fun)
# 
# #ordonnée gauche
# 
# #theoretical
# eq_th <- predict_th(A_fun,correl_fun,B_fun)
# #effective
# if (correl_fun=="RegNeg"|correl_fun=="CRPos"|correl_fun=="CRNeg") {
#   eq_eff <- predict_eff(E0_fun,B_fun,A_fun,correl_fun)
# } 
# #choose equilibrium as resident to compute RNV at equilibrium
# if (correl_fun=="SC"|correl_fun=="Comp"|correl_fun=="RegPos") {
#   E_eq <- eq_th$pred_e*Etot_fun
#   e_eq <- eq_th$pred_e
# } else {
#   E_eq <- eq_eff$pred_E
#   e_eq <- eq_eff$pred_e
# }
# 
# #ordonnée droite
# false_sim <- cbind(t(E_eq),t(A_fun),sum(E_eq),sum(A_fun),flux(E_eq,A_fun),t(A_fun))
# 
# RNV_all <- RNV.for.simul(false_sim,n_fun,N_fun,correl_fun,beta_fun,end.mean=FALSE)
# 
# #recup mean of RNV size
# #because mean for only one value is equal to this value, it is easier to use RNV_proxy (matrix) rather than RNV_size (list)
# RNV_size <- RNV_all$RNV_proxy
# 
# 
# #Mise en forme
# #only RNV of small mutation, even there is one RNV at equilibrium, due to computing method
# compar_at_eq <- cbind.data.frame(A=A_fun,RNV=RNV_size[1,],e=e_eq)
# compar_at_eq_ord <- compar_at_eq[order(compar_at_eq$A),]
# 
# 
# #Graphique
# x11()
# #enough size at right
# par(mar=c(4,4,3,5)) 
# 
# 
# ### first plot (in black)
# #line e fct A
# plot(compar_at_eq_ord$A, compar_at_eq_ord$e, pch=16, axes=F, ylim=c(0,1), xlab="", ylab="", type="o",col="black", main="Comparing e and RNV at equilibrium")
# #add left y-axis
# axis(2, ylim=c(0,1),col="black")
# #axis legend
# mtext("Relative concentrations",side=2,line=2.5,lwd=2.5)
# #encadre plot
# box()
# 
# 
# ### second plot (in red)
# #authorize superposition of plot
# par(new=T)
# #line RNV fct A
# plot(compar_at_eq_ord$A, compar_at_eq_ord$RNV, pch=15,  xlab="", ylab="", ylim=c(0,max(RNV_size[1,])), axes=F, type="o", col="red")
# #add right y-axis
# axis(4, ylim=c(0,max(RNV_size[1,])), col="red",col.axis="red")
# mtext("RNV size",side=4,col="red",line=2.5,lwd=2.5)
# 
# ### add x-axis
# axis(1,xlim=range(A_fun))
# #axis(1,pretty(range(A_fun),10)) #to pretty interval
# mtext("Activities",side=1,col="black",line=2.5,lwd=2.5)
# 
# ### legend
# legend(x="topleft",legend=c("e","RNV"),text.col=c("black","red"),pch=c(16,15),col=c("black","red"))
# 



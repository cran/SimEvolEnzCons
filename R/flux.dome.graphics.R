#' @name flux.dome.graphics
#' @aliases flux.dome.graph3D
#' 
#' 
#' 
#' @title Plot of the flux dome
#'
#' @description
#' Function \code{flux.dome.graph3D} gives a 3D graph of flux dome for a three-enzyme pathway.
#' 
#' Function \code{flux.dome.projections} gives a triangular diagram of flux dome, 
#' and if \code{add.reg==TRUE}, it also gives a plot of the curve defined by the intersection of the dome with the upright plan drawn from regulation line.
#'
#'
#'
#' @details
#' 
#' \bold{General}
#' 
#' \emph{\bold{Available only for three-enzyme pathway}, i.e. length of \code{A_fun} equal to 3.}
#'  
#' \bold{Flux dome exists only in case of the competition.}
#' 
#' \code{E_ini_fun} is rescaled by a cross product to have sum of \code{E_ini_fun} equal to \code{Etot_fun}.
#' 
#' Only flux dome is returned if \code{add.reg} is \code{FALSE}.
#' 
#' If \code{add.reg=TRUE}, \code{E_ini_fun} and \code{B_fun} are required. Else, default is \code{NULL}.
#' 
#' 
#' \bold{Function \code{flux.dome.graph.3D}}
#' 
#' If \code{add.reg} is \code{TRUE}, a section corresponding to the upright plan drawn from line on which the relative concentrations move when there is co-regulations \emph{(see \code{\link{droites}})} is added on graph.
#' 
#' 
#' \bold{Function \code{flux.dome.projections}}
#' 
#' Projection on plan \code{(e1,e2,e3)} is from blue for low flux to red for high flux.
#' Colors are taken from \code{Spectral} palette of the package '\code{RColorBrewer}'.
#' 
#' Because contour lines are calibrated to correspond to integer values of flux, it is possible to have a difference between input and output \code{nbniv}. 
#'  
#' If \code{add.reg} is \code{TRUE}, line on which the relative concentrations move when there is co-regulations \emph{(see \code{\link{droites}})} is added on triangular diagram.
#' Also, a new plot of the curve defined by the intersection of the dome with the upright plan drawn from this line is drawn.
#' 
#' 
#' 
#' @import rgl
#' @importFrom ade4 triangle.plot
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices dev.off
#' @importFrom grDevices dev.new
#' @importFrom graphics legend
#' @importFrom graphics points
#' @importFrom graphics abline
#' @importFrom graphics plot
#' 
#' 
#' 
#' @inheritParams droites
#' @inheritParams flux
#' @param Etot_fun Numeric value of total concentration
#' @param add.reg Logical. Add regulation line in plots? Default is \code{FALSE}. 
#'
#'
#' 
#' 
#' @seealso 
#' See function \code{\link{droites}} to compute regulation line. 
#'
#'
#'
#'
#' @examples
#' Etot <- 100
#' A <- c(1,10,30)
#' beta <- matrix(c(1,10,5,0.1,1,0.5,0.2,2,1),nrow=3)
#' B <- compute.B.from.beta(beta)
#' 
#' #or for B :
#' B <- 1/c(0.2,0.5,0.3)
#' E0<- c(30,30,30)
#' 
#' flux.dome.graph3D(A,Etot,add.reg=TRUE,B,E0)
#' flux.dome.projections(A,Etot,add.reg=TRUE,B,E0)
#'
#'
NULL




#' @rdname flux.dome.graphics
#' @usage flux.dome.graph3D(A_fun,Etot_fun,add.reg=FALSE,B_fun=NULL,
#' E_ini_fun=NULL,X_fun=1,marge_section=0.1,marge_top.dome=0.5)
#' @param marge_section Numeric. Height of section up to surface dome. Default is \code{0.1}
#' @param marge_top.dome Numeric. Space between top dome and max of zlim. Default is \code{0.5}
#' @return Function \code{flux.dome.graph3D} returns a 3D-graph of flux dome, with section in case of regulation. 
#' @export


flux.dome.graph3D <- function(A_fun,Etot_fun,add.reg=FALSE,B_fun=NULL,
                              E_ini_fun=NULL,X_fun=1,marge_section=0.1,marge_top.dome=0.5) {
  
  ######" settings
  n_poss <- 3
  n_fun <- length(A_fun)
  correl_fun <- name.correl(TRUE,add.reg,B_fun)
  # if (add.reg==TRUE) {
  #   if (any(B_fun)<0) { correl_fun <- "CRNeg" } else { correl_fun <- "CRPos" }
  # } else {
  #   correl_fun <- "Comp"
  # }
  E_ini_fun <- E_ini_fun*Etot_fun/sum(E_ini_fun)
  
  
  ###### Test
  if (n_fun!=n_poss) {
    stop("Available only for 3 enzymes.")
  }
  if (correl_fun!="Comp"&correl_fun!="CRPos"&correl_fun!="CRNeg") {
    stop("Flux dome exists only in case of competition.")
  }
  is.B.accurate(B_fun,n_fun,correl_fun)
  if (length(E_ini_fun)==0&correl_fun!="Comp") {
    stop("In case of regulation, initial point 'E_ini_fun' is needed.")
  }
  if (length(E_ini_fun)!=n_fun&correl_fun!="Comp") {
    stop("E_ini_fun and A_fun need to have the same length.")
  }
  
  
  #####" Initialization
  
  #coordinates of dome maximal point = equilibrium in case of competition only
  max_max_e <- predict_th(A_fun,"Comp")$pred_e #(A_fun^(-1/2)/sum(A_fun^(-1/2)))
  max_max_J <- flux(max_max_e*Etot_fun,A_fun,X_fun)
  
  #in case of competition and regulation
  if (correl_fun=="CRPos"|correl_fun=="CRNeg") {
    ### Remarkable points
    #initial concentrations
    e0 <- E_ini_fun/sum(E_ini_fun)
    J0 <- flux(E_ini_fun,A_fun,X_fun)
    
    #theoretical equilibrium
    eq_th <- predict_th(A_fun,correl_fun,B_fun)
    J_th <- flux(eq_th$pred_e*Etot_fun,A_fun,X_fun)
    
    #effective equilibrium
    eq_eff <- predict_eff(E_ini_fun,B_fun,A_fun,correl_fun)
    J_eff <- flux(eq_eff$pred_E,A_fun,X_fun)
     
    
    ### Line of concentrations
    #bounds of tau to have e between 0 and 1
    borne_tau <- range_tau(E_ini_fun,B_fun)
    #all tau of the line with a safety limits
    tau <- seq(min(borne_tau)+0.0001,max(borne_tau)-0.0001,by=0.01)
    #computes E
    # E_test_ens <- matrix(0,nrow=length(tau),ncol=n_fun) #calcul de E
    # for (j in 1:length(tau)) {
    #   E_test_ens[j,] <- droite_E.CR(tau[j],E_ini_fun,B_fun)
    # }
    E_test_ens <- t(sapply(tau,droite_E.CR,E_ini_fun,B_fun))
    #line coordinate e
    e_test <- E_test_ens/sum(E_ini_fun)
    #corresponding flux J
    J_test <- apply(E_test_ens,1,flux,A_fun,X_fun)
    
  }
  
  

  #######" Value for shape
  
  # graphic parameters
  sph.radius <- 0.0125*(max_max_J*0.9) #sphere radius for remarkable points
  
  # dome surface
  e_seq <- seq(0,1,by=0.01) 
  J_dome <- matrix(nrow=length(e_seq),ncol=length(e_seq))
  #every value of e1
  for (x1 in 1:length(e_seq)) {
    #every value of e2
    for (x2 in 1:length(e_seq)) { 
      #because for n=3, e1+e2+e3 = 1 and e3 >=0
      if ((1-e_seq[x1]-e_seq[x2])>=0) {
        #E1 , E2 , and E3=1-E1-E2 for n=3
        E_vec <- c(e_seq[x1]*Etot_fun,e_seq[x2]*Etot_fun,(1-e_seq[x1]-e_seq[x2])*Etot_fun)
        J_dome[x1,x2] <- flux(E_vec,A_fun,X_fun)
      } else {
        #to avoid negative value
        J_dome[x1,x2] <- 0
      }
    }
  }
  
  #section of dome by line
  if (correl_fun=="CRPos"|correl_fun=="CRNeg") {
    #for w_plan : each line = a point ; each column = coordinates (x,y,z) + relative distance
    #we need for 4 points(i.e. 4 lines) to have a quadrilatere
    #here, coordinates are (e1,e2,J)
    w_plan <- matrix(0,ncol=4,nrow=4)
    
    #computes e at bounds of tau -> borne_e have 2 lines and 3 columns, one for each coordinates (e1,e2,e3)
    borne_e <- t(sapply(range(tau), droite_e,E_ini_fun,B_fun))
    
    #coordinates e1 and e2 for each points 
    #points 1 and 3 = min tau ; points 2 and 4 = max tau
    w_plan[1:2,1:2] <- borne_e[,1:2]
    w_plan[3:4,1:2] <- borne_e[,1:2] #due to fill system of matrix in R
    
    #coordinate J -> height of the section
    #0 to the base and max on the line is effective equilibrium
    #points 1 and 2 = base ; points 3 and 4 = max height plus a margin
    w_plan[,3] <- rep(c(0,J_eff+marge_section),each=2) 
    
    #same weight to each point ?
    w_plan[,4] <- 1
    
    #to have each point in column and coordinates in line
    w_plan <- t(w_plan)
    
    #order to rely points (1=min base -> 2=max base -> 4=max J_eff -> 3=min J_eff)
    w_indices <- c(1,2,4,3)
  } else {
    w_plan <- NaN
    w_indices <- NaN
  }
  
  
  
  #######" Graphics
  
  open3d()
  
  ## Graph 3D
  #plot3d(col_niv[,1], col_niv[,2], col_niv$color, zlim=c(0,max_max_J+marge_section),  xlab="Enzyme1", ylab="Enzyme2", zlab="Flux", col="lightgrey",type='n') #dôme J=f(E1,E2)
  #null plot
  plot3d(0,0,0,zlim=c(0,max_max_J+marge_top.dome),xlim=c(0,1),ylim=c(0,1),xlab="Enzyme 1", ylab="Enzyme 2", zlab=" ")
  #dome surface
  surface3d(e_seq,e_seq,J_dome,col='lightgrey',zlim=c(0,max_max_J+marge_top.dome),xlim=c(0,1),ylim=c(0,1))
  #section in case of regulation
  shade3d(qmesh3d(w_plan,w_indices),col='oldlace')
  
  ## Remarkable points
  # max of J without regulation
  spheres3d(max_max_e[1], max_max_e[2], max_max_J, col="violet", point_antialias=TRUE, radius=sph.radius, size=18)
  if (correl_fun=="CRPos"|correl_fun=="CRNeg") {
    # line visible on dome
    #lines3d(e_test[,1],e_test[,2],J_test,lwd=6,col="darkgreen",line_antialias=TRUE)
    #theoretical equilibrium e*=1/B
    spheres3d(eq_th$pred_e[1], eq_th$pred_e[2], J_th, col="red", point_antialias=TRUE, radius=sph.radius, size=18)
    #max J on line
    spheres3d(eq_eff$pred_e[1], eq_eff$pred_e[2], J_eff, col="yellow1", point_antialias=TRUE, radius=sph.radius, size=18)
    #initial point e0
    spheres3d(e0[1], e0[2], J0, col="black", point_antialias=TRUE, radius=sph.radius, size=18)
    
  }

  ## Projection on plane E
  points3d(max_max_e[1], max_max_e[2], 0, col="violet", point_antialias=TRUE, size=10) #max of J without regulation
  if (correl_fun=="CRPos"|correl_fun=="CRNeg") {
    lines3d(e_test[,1], e_test[,2], 0, lwd=5, col="darkgreen") #projection of line on concentration plane (e1,e2)
    points3d(eq_th$pred_e[1], eq_th$pred_e[2], 0, col="red", point_antialias=TRUE, size=10) #point 1/B=e*
    points3d(eq_eff$pred_e[1], eq_eff$pred_e[2], 0, col="yellow1", point_antialias=TRUE, size=10) #max of J on line
    points3d(e0[1], e0[2], 0, col="black", point_antialias=TRUE, size=10) #initial point
  }
  
}







#############" Projection of flux dome


#' @rdname flux.dome.graphics
#' @usage flux.dome.projections(A_fun,Etot_fun,add.reg=FALSE,B_fun=NULL,
#' E_ini_fun=NULL,X_fun=1,nbniv=9,niv.palette=NULL,marge_section=0.1,
#' posi.legend="topleft",new.window=FALSE,cex.pt=1.5,...)
#' @param nbniv Numeric. Number of contour lines, between 3 and 11. Default is \code{9}
#' @param new.window Logical. Does graphics appear in a new window ? Default is \code{TRUE}
#' @param cex.pt Numeric. Size of points in plot(s) drawn by \code{flux.dome.projections}
#' @param ... Arguments to be passed to \code{plot} function, such as \code{lwd} or \code{cex.axis}
#' @param niv.palette Character vector. Color palette to be passed to contour lines.
#' @param posi.legend Contour line legend position. See \emph{details}.
#' 
#' @details 
#' \bold{Color palette for contour lines}
#' 
#' By default, color palette for contour lines is taken in palette \code{"Spectral"} of package \code{RColorBrewer}, with cool colors for low flux values and warm colors for high ones.
#' 
#' Input your own color palette in \code{niv.palette}.
#' 
#' If color numbers in \code{niv.palette} is inferior to the number of contour lines \code{nbniv}, a warning message is printed and palette is replicated to have enough colors.
#' 
#' \code{posi.legend} is the position of the legend for contour line.
#' It is a character vector of length 1 (\code{"topleft"}, etc.) or a numeric vector of length 2, which is the coordinates of the upper left corner of the legend box.
#' If \code{NULL}, the legend will not be shown.
#' The first (resp. second) element of \code{posi.legend} is passed to argument \code{x} (resp. \code{y}).
#' Note that the minimum and maximum coordinates for the return triangular plot are between -0.864 and 0.864 for x-axis, and -0.66 and 1.064 for y-axis. Use \code{locator(1)} to adjust position legend.
#' 
#' @return Function \code{flux.dome.projections} returns a triangular diagram corresponding to the projection of the flux dome on the plan of relative concentrations.
#' In the case of regulation (\code{add.reg=TRUE}), a straight line is drawn.
#' If \code{add.reg=TRUE}, function \code{flux.dome.projections} also returns a plot of the curve defined by the intersection of the dome with the upright plan drawn from regulation line.
#' 
#' @examples 
#' #Position of legend at right
#' flux.dome.projections(c(1,10,30),100,posilegend=c(0.55,0.9))
#' 
#' @export


flux.dome.projections <- function(A_fun,Etot_fun,add.reg=FALSE,B_fun=NULL,
                              E_ini_fun=NULL,X_fun=1,nbniv=9,niv.palette=NULL,marge_section=0.1,
                              posi.legend="topleft",new.window=FALSE,cex.pt=1.5,...) {
  
  ######" settings
  n_poss <- 3
  n_fun <- length(A_fun)
  correl_fun <- name.correl(TRUE,add.reg,B_fun)
  # if (add.reg==TRUE) {
  #   if (any(B_fun)<0) { correl_fun <- "CRNeg" } else { correl_fun <- "CRPos" }
  # } else {
  #   correl_fun <- "Comp"
  # }
  E_ini_fun <- E_ini_fun*Etot_fun/sum(E_ini_fun)
  
  
  ###### Test
  if (n_fun!=n_poss) {
    stop("Available only for 3 enzymes.")
  }
  if (correl_fun!="Comp"&correl_fun!="CRPos"&correl_fun!="CRNeg") {
    stop("Flux dome exists only in case of competition.")
  }
  is.B.accurate(B_fun,n_fun,correl_fun)
  if (length(E_ini_fun)==0&add.reg==TRUE) {
    stop("In case of regulation, initial point 'E_ini_fun' is needed.")
  }
  if (length(E_ini_fun)!=n_fun&add.reg==TRUE) {
    stop("E_ini_fun and A_fun need to have the same length.")
  }
  
  
  #####" Initialization
  
  #coordinates of dome maximal point = equilibrium in case of competition only
  max_max_e <- predict_th(A_fun,"Comp")$pred_e #(A_fun^(-1/2)/sum(A_fun^(-1/2)))
  max_max_J <- flux(max_max_e*Etot_fun,A_fun,X_fun)
  
  #in case of competition and regulation
  if (correl_fun=="CRPos"|correl_fun=="CRNeg") {
    ### Remarkable points
    #initial concentrations
    e0 <- E_ini_fun/sum(E_ini_fun)
    J0 <- flux(E_ini_fun,A_fun,X_fun)
    
    #theoretical equilibrium
    eq_th <- predict_th(A_fun,correl_fun,B_fun)
    J_th <- flux(eq_th$pred_e*Etot_fun,A_fun,X_fun)
    
    #effective equilibrium
    eq_eff <- predict_eff(E_ini_fun,B_fun,A_fun,correl_fun)
    J_eff <- flux(eq_eff$pred_E,A_fun,X_fun)
    
    
    ### Line of concentrations
    #bounds of tau to have e between 0 and 1
    borne_tau <- range_tau(E_ini_fun,B_fun)
    #all tau of the line with a safety limits
    tau <- seq(min(borne_tau)+0.0001,max(borne_tau)-0.0001,by=0.01)
    #computes E
    # E_test_ens <- matrix(0,nrow=length(tau),ncol=n_fun) #calcul de E
    # for (j in 1:length(tau)) {
    #   E_test_ens[j,] <- droite_E.CR(tau[j],E_ini_fun,B_fun)
    # }
    E_test_ens <- t(sapply(tau,droite_E.CR,E_ini_fun,B_fun))
    #line coordinate e
    e_test <- E_test_ens/sum(E_ini_fun)
    #corresponding flux J
    J_test <- apply(E_test_ens,1,flux,A_fun,X_fun)
    
  }
  
  #######" Value for shape

  
  #spacing between contour lines
  saut <- ceiling(max_max_J/nbniv)
  #ceiling : to have integer
  #saut <- max_max_J/nbniv #will be more exact in term of nbniv

  
  ## Values for contour lines
  #values of J for contour lines
  J_select <- seq(0.0001,max_max_J,by=saut)
  param.nbniv <- nbniv #save
  nbniv <- length(J_select) #to avoid future problems
  
  #tested values of e to search coordinates (e1,e2,e3) for a corresponding J 
  e_seq <- seq(0.01,1,by=0.01)

  #A list of coordinates: each element of the list is a vector corresponding to one contour line
  e_curve1 <- vector("list",length=nbniv)
  #for each chosen level of J
  for (k in 1:nbniv) {
    #value of J for this level/contour line
    J_niv <- J_select[k]
    #all possible coordinates for this level
    sol_niv <- NULL
    
    #for each enzyme, to scan a maximum of coordinates triplet
    for (i in 1:n_fun) {
      #i= number of set variable/enzyme : relative concentration x
      #j = number of searched variable/enzyme : relative concentration y
      #l = number of deduced variable/enzyme from the two others : relative concentration 1-x-y
      ## cf forme of equation e_y = fct(e_x,Jniv) to understand importance of enzyme placement (by transforming equation Jniv = f(e1,e2,e3))

      # j <- (i+1) %% n_fun
      # l <- (i+2) %% n_fun
      
      if (i ==1) { #if variable is enzyme 1
        j <- 2 ; l <- 3 }
      if (i==2) { j <- 3 ; l <- 1 }
      if (i==3) { j <- 1 ; l <- 2 }
      
      for (x in e_seq) {
        #x = concentration of set variable/enzyme
        a <- A_fun[j]*A_fun[l]*(X_fun*Etot_fun/J_niv - 1/(A_fun[i]*x))
        b <- A_fun[j]-A_fun[l]-a*(1-x)
        c <- A_fun[l]*(1-x)
        
        #to avoid infinite solution
        if (a!=0) {
          #solve 2nd degre polynom:  a y^2 + b y + c = 0
          sol_2dg <- solv.2dg.polynom(a,b,c)
          
          #relative concentrations are all positive
          if (all(sol_2dg>=0)&all((1-x-sol_2dg)>=0)) {
            #to keep value : in column = enzyme number ; in line = each solution for 2nd equation (0,1 or 2 solutions)
            sol_e <- matrix(0,nrow=length(sol_2dg),ncol=n_fun)
            #value of set enzyme
            sol_e[,i] <- x
            #value of searched enzyme
            sol_e[,j] <- sol_2dg
            #value of deduced enzyme, because e1+e2+e3=1
            sol_e[,l] <- 1-x-sol_2dg 
            
            #add these solutions to the precedent ones
            sol_niv <- rbind(sol_niv,sol_e)
          }
        }
      }
    }
    #put solutions to the corresponding J contour line
    e_curve1[[k]] <- sol_niv
  }
 
  #unlist e_curve1
  col_niv <- NULL
  #for each level/contour line
  for (k in 1:nbniv) {
    #take value for this level
    matniv <- e_curve1[[k]]
    #add a new column with the level number
    matniv <- cbind(matniv,k)
    #add to values of precedent level
    col_niv <- rbind(col_niv,matniv)
  }
  col_niv <- as.data.frame(col_niv)
  #add names to the columns
  colnames(col_niv) <- c("e1","e2","e3","color")
  
  
  
  #######" Graphics
  #######" Projection on plan e
  
  
  ## Graph parameters
  #color palette to contour line
  if (is.null(niv.palette)) {
    #default value for palette
    niv.palette<-brewer.pal(nbniv,"Spectral")
    #reverse colors, to have red=top and blue=down
    niv.palette<-rev(niv.palette)
  }
  #not enough colors in palette
  if (length(niv.palette)<nbniv) {
    warning("Not enough colors for contour line.")
    #replicate palette to have enough colors
    niv.palette <- rep(niv.palette,ceiling(nbniv/length(niv.palette)))
  }
  
  
   #Triangular diagram of e1,e2,e3 on section (e1+e2+e3=1)
  
  #open graphic window
  if (new.window==TRUE) {
    dev.new()
  }
  
  
    #position of maximal value of J
    tr_Jm <- triangle.plot(t(max_max_e),scale=FALSE,show.position=FALSE,draw.line = FALSE,cpoint=0)
  
    if (add.reg==TRUE) {
      #straight line
      tr_test <- triangle.plot(e_test[,1:3],scale=FALSE,show.position=FALSE,draw.line = FALSE,cpoint=0)
      # points(tr_test,col=1,type='l') #droite noire
      
      #initial point
      tr_ini <- triangle.plot(t(e0),scale=FALSE,show.position=FALSE,draw.line = FALSE,cpoint=0)
      # points(tr_ini,col=1,pch=17) #black triangle
      
      #effective equilibrium
      tr_eff <- triangle.plot(t(eq_eff$pred_e),scale=FALSE,show.position=FALSE,draw.line = FALSE,cpoint=0)
      #points(tr_eff,col=1,pch=13) #black cross-point
      
      #theoretical equilibrium
      if (correl_fun=="CRPos"|correl_fun=="RegPos") {
        #because unaccessible in case of negative regulation
        tr_th <- triangle.plot(t(eq_th$pred_e),scale=FALSE,show.position=FALSE,draw.line = FALSE,cpoint=0)
        # points(tr_th,col=1,pch=8) #black star
      }
    }
    
    #close graphic window
  #dev.off()


    if (new.window==TRUE) {
      dev.new()
    }
  #pour traits : type='p',lwd=2.5,cex=0.3
  #pour aires : type='o',lwd=5,cex=1
  
  #trianguar diagram for contour lines
  tr_niv <- triangle.plot(col_niv[,1:3],scale=FALSE,show.position=FALSE,draw.line = FALSE,cpoint=0) #diagramme triangulaire de e1,e2,e3 pour toutes les simulations (convertis les points en coordonnées xy)
  # for each contour line
  for (k in 1:nbniv){
    #color points to corresponding level
    points(tr_niv[col_niv$color==k,],col=niv.palette[k],cex=0.3,pch=20,type='o',lwd=5) #couleur de la simulation i sur les points correspondant à la simulation i sur le diagramme triangulaire
  }
  
  # Remarkable points
  #points(tr_Jm,pch=20,col='violet')
  #points(tr_Jm,pch=8,col='violet',cex=cex.pt,lwd=1.5)
  #maximal value of J -> violet point
  points(tr_Jm,pch=21,bg='violet',cex=cex.pt)
  if (add.reg==TRUE) {
    #straight line -> black line
    points(tr_test,col=1,type='l',cex=cex.pt) #droite noire
    #initial point -> black point
    points(tr_ini,col=1,pch=21,cex=cex.pt,bg="black") #point noir #black triangle pour col=1,pch=17
    #effective equilibrium -> yellow point
    points(tr_eff,col=1,pch=21,cex=cex.pt,bg="yellow2") #point jaune #black cross-point pour pch=13
    #theoretical equilibrium -> red point
    if (correl_fun=="CRPos") {
      points(tr_th,col=1,pch=21,cex=cex.pt,bg="firebrick1") #point rouge  #red star pour col=2,pch=8,lwd=1.5
    }
  }

  
  # Legend
  if (!is.null(posi.legend)) {
    # if want a legend
    legend(posi.legend[1],posi.legend[2],legend=round(J_select),lwd=rep(5,nbniv),col=niv.palette,bty="n",cex=1.2)
  }
  
  
  ######" Graphics
  ##### Projection on section J=f(tau)
  
  if (add.reg==TRUE) {
    if (new.window==TRUE) {
      dev.new()
    }
    
    #line -> black line
    plot(x=tau,y=J_test,type='l',col=1,xlab="Tau",ylab="Flux",pch=20,ylim=c(0,J_eff+marge_section),...)#,lwd=1.2
    
    ## Remarkable point
    #initial point -> black point
    points(0,J0,pch=21,cex=cex.pt,bg="black") #initial "black triangle pch=17
    #effective equilibrium -> yellow point
    abline(v=eq_eff$pred_tau,lty=2)
    points(eq_eff$pred_tau,J_eff,pch=21,cex=cex.pt,col=1,bg="yellow2") #pch=13=cross-point
    #theoretical equilibrium -> red point
    if (correl_fun=="CRPos") {
      abline(v=1)
      points(1,J_th,pch=21,cex=cex.pt,bg="firebrick1") #eq th #red star pch=8,lwd=1.5
    }
    
  }
  
  
}







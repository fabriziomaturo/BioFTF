radius_plot=function(x,n=20){

  #################################################################################
  ############################      DATA    #######################################
  #################################################################################
  library(BioFTF)
  for (i in 1:ncol(x)) colnames(x)=paste("species",c(1:ncol(x)),sep=" n.")
  for (i in 1:nrow(x)) rownames(x)=paste("community",c(1:nrow(x)),sep=" n.")
  x

  class(x)
  x=as.matrix(x)

  nsiti=nrow(x)
  nspecie=ncol(x)

  # matrix with relative abundance

  xrel=prop.table(x,1)


  #################################################################################
  ###############################     DOMAIN    ###################################
  #################################################################################

  # fix points of domain

  k=(n-1)/2
  b=seq(-1,1,1/k)
  b[b==0]=0.001    #avoid "shannon" jump close to 0
  b[1]=-0.9999  #avoid "richness" jump close to -1

  # temp matrix length(b)

  appo.beta=matrix(rep(b),length(b),nsiti)


  #################################################################################
  ##############################     PROFILES    ##################################
  #################################################################################

  # temp matrix for computing beta profile

  appo.profile=matrix(NA,length(b),nsiti)

  appo.siti=matrix(rep((xrel),each=length(b)),nsiti*length(b),nspecie)

  # compute beta profile

  all=length(b)*nsiti

  for(i in 1:all) {appo.profile[i]=(1-sum(appo.siti[i,]^(appo.beta[i]+1)))/(appo.beta[i])}
  appo.profile[1,]=round(appo.profile[1,])  #the richness must be integer


  #################################################################################
  #########################       FIRST           #################################
  #################################################################################

  # temp matrix for computing betaprofile first derivative

  appo.first=matrix(NA,length(b),nsiti)

  # compute beta profile first derivative without species with null frequencies


  appo.siti2=log(appo.siti)
  appo.siti2[appo.siti2==-Inf]=0


  for (i in 1:all) appo.first[i]=(sum(appo.siti[i,]^(appo.beta[i]+1)*(1-appo.beta[i]*appo.siti2[i,]))-1)/((appo.beta[i])^2)


  #################################################################################
  #############################     SECOND     ####################################
  #################################################################################

  # temp matrix for computing betaprofile second derivative

  appo.second=matrix(NA,length(b),nsiti)

  # compute beta profile second derivative without species with null frequencies

  for (i in 1:all) appo.second[i]=(sum((appo.siti[i,]^(appo.beta[i]+1))*(-(appo.beta[i])^2*((appo.siti2[i,])^2)-2+2*appo.beta[i]*appo.siti2[i,]))+2)/((appo.beta[i])^3)


  #################################################################################
  ##############################      RADIUS      #################################
  #################################################################################

  # temp matrix for computing beta profile radius

  appo.radius=matrix(NA,length(b),nsiti)

  # compute beta profile radius

  for (i in 1:all) appo.radius[i]=((1+appo.first[i]^2)^(3/2))/abs(appo.second[i])


  # plot beta profile radius

  for (i in 1:nsiti)
  {plot(b,appo.radius[,i],type="l",xlim=c(-1,1),ylim=c(0,max(appo.radius)),lty=i,lwd=1,col=i,xlab="Beta",ylab="Radius",main="Radius of curvature")
    par(new=TRUE)}
  legend(1,max(appo.radius), paste("Com.",c(1:nsiti)),lty=c(1:nsiti),y.intersp=0.65,ncol=2,lwd=1,col=c(1:nsiti),xjust=1,merge=TRUE)

}



#**************************************************
readRules<-function(rulefile){
  # read in the rules data file
  rules.fields.names<<-names(read.table(rulefile,colClasses="character",header=TRUE,skip=1,nrows=1,sep=','))
  rules.colClasses.tmp <- (read.table(rulefile,colClasses="character",header=TRUE,skip=1,nrows=1,sep=','))
  rules.colClasses<<-matrix(NA,ncol=2,nrow=length(rules.colClasses.tmp))
  
  for(col.val in 1:length(rules.colClasses.tmp)){
    rules.colClasses[col.val,] <<- matrix((unlist(strsplit(as.character(rules.colClasses.tmp[col.val]),split=";"))),ncol=2,byrow=T)  
    rules.colClasses[rules.colClasses[,2]==rules.colClasses[,1],2]<<-"1"
  }
  
  rules.raw<-read.table(rulefile,colClasses="character",skip=3,col.names =rules.fields.names,sep=',')
  
}

#**************************************************
getParams<-function(modelParams){
  rules<-vector("list",ncol(modelParams))
  names(rules)<-rules.fields.names
  
  for(col.val in 1:nrow(rules.colClasses))
    if(rules.colClasses[col.val,1]=="character"){
      rules[[col.val]] <- matrix(suppressWarnings(unlist(strsplit(as.character(modelParams[col.val]),split=";"))), ncol=as.numeric(rules.colClasses[col.val,2]),byrow=T)
    }else if(rules.colClasses[col.val,1]=="numeric"){
      rules[[col.val]]<- matrix(suppressWarnings(as.numeric(unlist(strsplit(as.character(modelParams[col.val]),split=";")))), ncol=as.numeric(rules.colClasses[col.val,2]),byrow=T)
    }else if(rules.colClasses[col.val,1]=="logical"){
      rules[[col.val]]<- matrix(suppressWarnings(as.logical(unlist(strsplit(as.character(modelParams[col.val]),split=";")))), ncol=as.numeric(rules.colClasses[col.val,2]),byrow=T)
    }
  return(rules)
}

#***********************************************************
#Derives Ricker projections of recruitment with random error
Ricker<-function(S,a,b,cv){			 	
  err<-rnorm(1,0,cv)#rnorm(1,-cv^2/2,cv)
  if(a>=0){
    if(b!=0) R<-S*exp(a-b*S)*exp(err)
    if(b==0)R<-S*exp(err)
  }
  if(a<0)R<-S*exp(a)*exp(err)
  return(list(R=R,err=err))
}

RickerAR1<-function(S,utminus1,a,b,cv,rho){			# where a,b,cv are in lognormal space 	
  err.1<-rnorm(1,0,cv)
  err<-utminus1*rho+err.1
  if(a>=0) {
    if(b!=0) {
      R<-S*exp(a*(1-S/b))*exp(err)*exp(-cv^2/2)
      err.next<-log(R/S)-(a*(1-S/b))+cv^2/2
    }
    if(b==0) {
      R<-S*exp(err)
      err.next<-log(R/S)-0
    }
  }
  if(a<0){
    R<-S*exp(a)*exp(err)*exp(-cv^2/2)
    err.next<-log(R/S)-0+cv^2/2
  }
  if(R=="Inf")print(c("Error in RickerAR1, err.1=",err.1))
  return(list(R=R,err.next=err.next))
}

#Derives Ricker projections of recruitment without random error (for plots)
RickerDet<-function(S,a,b=NULL){			 	
  err <- 0#rep(-cv^2/2, length(S))
  if (length(a) > 1) {
    b <- a[[2]]
    a <- a[[1]]
  }
  if(a>=0){
    if(b!=0) R<-S*exp(a-rep(b, length(S))*S)*exp(err)
    if(b==0) R<-S*exp(err)
  }
  if(a<0)R<-S*exp(a)*exp(err)
  return(R)
}
#**********************************************************


#___________________________________________________________________________________________________________
# Estimate parameters or Ricker model

inits <- function(R,S){
  lin<-lm(log(R/S)~S)
  yi.init<-lin$coef[1]   #Ricker a #y-intercept
  sl.init<-lin$coef[2]   #Ricker b # slope
  sig.init<-log(sd(resid(lin)))
  return(list(a.init=as.numeric(yi.init), b.init=-as.numeric(sl.init), sig.init=as.numeric(sig.init)))
}
rickerModel <- function(theta,R,S){
  a <- as.numeric(theta[1])
  #a <- exp(as.numeric((theta[1])))
  b <- (as.numeric(theta[2]))
  sig <- exp(as.numeric(theta[3]))
  PR=S*exp(a-b*S)
  #pRS=a-b*S
  epsilon.wna=log(R)-log(PR)	#residuals
  #epsilon.wna=log(R/S)-log(pRS)	#residuals
  #cat(epsilon.wna)
  epsilon=as.numeric(na.omit(epsilon.wna))
  nloglike=sum(dnorm(epsilon,0,sig, log=T))
  return(list(PR=PR, epsilon=epsilon, nloglike=nloglike)) #actually returns postive loglikelihood
}

fn1 <- function(theta, R,S){  -1.0*rickerModel(theta, R,S)$nloglike }

rickerSolver <- function(R,S){
  init.vals<-c(inits(R,S)[1], 1,inits(R,S)[3])#use 1 as initial value for log(b), and the slope of lm won't always be positive!
  SRfit=optim(par=init.vals,fn=fn1,R=R,S=S, method="BFGS", hessian=T)	#CH: hessian=2nd derivative, optim good for up to 40 parameters
  #SRfit$par[1]  <- exp(SRfit$par[1])  # back-transform Ricker-a, which was exponentiated in Ricker( )
  SRfit$par[1]  <- (SRfit$par[1])  
  SRfit$par[2]  <- (SRfit$par[2])  
  V=solve(SRfit$hessian) #covariance matrix
  std=sqrt(abs(diag(V)))
  X=V/(std%o%std)	#outer product of standard dev matrix (CH comment)
  return(list(SRfit=SRfit, etheta=SRfit$par, V=V, X=X))
}
#*********************************************************************

sGenOptimum <- function(S, theta, sMSY) {
  a = theta[1]
  b = theta[2]
  sig = exp(theta[3])
  prt <- S * exp(a - b * S)
  epsilon <- log(sMSY) - log(prt)
  nLogLike <- sum(dnorm(epsilon, 0, sig, log = T))
  
  return(list(prt = prt, epsilon = epsilon, nLogLike = nLogLike, S = S))
}


sGenSolver <- function (theta1, theta2, theta3, sMSY) {
  #gives the min Ricker log-likelihood
  theta <- c(theta1,theta2,theta3)
  fnSGen <- function(S, theta, sMSY) -1.0 * sGenOptimum(S, theta, sMSY)$nLogLike
  
  fit <- optimize(f = fnSGen, interval = c(0, ((theta[1] / theta[2]) * (0.5 - 0.07 * theta[1]))),
                  theta = theta, sMSY = sMSY)
  #return(list(fit = fit$minimum))
  return(fit$minimum)
}

sGenSolver_v <- Vectorize(sGenSolver)

sMSYfunc<-function(a,b,ricSig=1){
  sMSY <-(1 - gsl::lambert_W0(exp(1 - a))) / b
  return(sMSY)
}  
sMSYfunc_v <- Vectorize(sMSYfunc)


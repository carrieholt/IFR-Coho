#Code to generate KF Ricker a values for 5 IFR Coho CUs

library(dplyr) 

#GetIFRcoho<-function(){
  data<-as.tbl(read.csv("data/Stockrecruitdatausedfor2019excl2010Chilko.csv", 
                        col.name=c("stk",	"CU",	"yr",	"cyc", "rec2", "rec3",
                                   "rec4",	"rec5",	"rec",	"eff",	"ETS",	
                                   "juv",	"Smax_sR",	"Smax_tvR",
                                   "n>Smax_sR",	"n>Smax_tvR",	"S>Smax",	"S0",
                                   "0.2S0",	"n<0.2S0"), skip=0))
  
  cu<-as.tbl(read.csv("data/cuToStock.csv", col.name=c("cuName","stock", 
                                                       "stYear", "endYear"),
                      skip=0))
  ricALR <- ricBLR <- ricLRAICc <- sMaxLR <- sMSYLR <- sGenLR <- ricBKF <- 
    ricKFAICc <- sMaxKF <- ricAKFMean <- ricAKFLast <- snr <- nllKF <- varKF <- 
    sMSYKFMean <- sMSYKFLast <- sGenKFLast <- sGenKFMean <- sGenKFLast <-
    ppnHigh <- ppnLow <- ppnHRec <- nYrs <- ricSig <- ricSige <- ricSigw <- NA
  
  maxYears<-max(cu$endYear)-min(cu$stYear)+1
  y1<-min(cu$stYear)
  
  ricAKF_mat <- ricAKFVar_mat <- sMSYKF_mat <- sGenKF_mat <- 
    matrix(NA, nrow=maxYears, ncol=length(cu$cuName))
  
  for(i in 1:length(cu$cuName)){
    
    R <- S <- NA
    R<-data%>%filter(CU==as.character(cu$cuName[i]))%>%
      filter(yr>=cu$stYear[i]&yr<=cu$endYear[i])%>%select(rec)
    S<-data%>%filter(CU==as.character(cu$cuName[i]))%>%
      filter(yr>=cu$stYear[i]&yr<=cu$endYear[i])%>%select(ETS)
    
    nYrs[i] <- length(S$ETS)
    out.lm <- lm(log(R$rec/S$ETS)~S$ETS)
    ricALR[i] <- out.lm$coef[1]
    ricBLR[i] <- -out.lm$coef[2]
    nllLR <- log (sum (out.lm$residuals^2) / length (out.lm$residuals)) *
      (length (out.lm$residuals)/2) + length(out.lm$residuals)/2 #concentrated neg log-likelihood = loge(sig^2)
    ricSig[i] <- sqrt(mean(out.lm$residuals^2))
    nAIC <- length(R)
    kAIC <- 2  #number of  parameters
    ricLRAICc[i] <- 2 * nllLR + 2 * kAIC + 2 * kAIC * (kAIC + 1) / (nAIC - kAIC - 1) #equation for AICc
    sMaxLR[i] <- 1/ricBLR[i]   
    sMSYLR[i] <- (1 - gsl::lambert_W0(exp(1 - ricALR[i]))) / ricBLR[i]
    sGenLR[i] <- sGenSolver(ricALR[i], ricBLR[i], 1, sMSYLR[i])
    
    
    
    if(is.na(ricBLR[i])==FALSE) {initial <- list(ricBLR[i], log(1), log(1), ricALR[i], 1, 1, "True")}
    if(is.na(ricBLR[i])==TRUE) {initial <- list(-1/max(R$rec, na.rm=T), log(1), log(1), 1, 1, 1, "True")}
    names(initial) <- c("b", "ln.sig.e", "ln.sig.w", "mean.a", "var.a", "Ts", "EstB")
    output <- kf.rw(initial=initial,x=S$ETS,y=log(R$rec/S$ETS))
    
    if(output$b<0){ #If SMax is positive
      ricAKF <- output$smoothe.mean.a
      ricAKFVar <- output$smoothe.var.a
      ricBKF[i] <- -output$b
      ricKFAICc[i] <- output$AICc
      #ricBKF.LCL[i] <- output$bCLs[1]
      #ricBKF.UCL[i] <- output$bCLs[2]
      sMaxKF[i]<-1/ricBKF[i]
      
      ricAKFMean[i] <- mean(ricAKF)
      ricAKFLast[i] <- ricAKF[length(ricAKF)]
      ricAKF_mat[(cu$stYear[i]-y1+1):(cu$endYear[i]-y1+1),i] <- ricAKF
      ricAKFVar_mat[(cu$stYear[i]-y1+1):(cu$endYear[i]-y1+1),i] <- ricAKFVar
      snr[i] <- output$sig.w/output$sig.e
      nllKF[i] <- output$cum.neg.log.lik
      varKF[i] <- output$sig.e^2 + output$sig.w^2
      ricSige[i] <- output$sig.e
      ricSigw[i] <- output$sig.w
      sMSYKF_mat[(cu$stYear[i]-y1+1):(cu$endYear[i]-y1+1), i] <- 
        (1 - gsl::lambert_W0(exp(1 - ricAKF_mat[(cu$stYear[i]-y1+1):(cu$endYear[i]-y1+1),i]))) / 
        rep(ricBKF[i],length(ricAKF))
      sMSYKFMean[i] <- sMSYKF_mat[(cu$stYear[i]-y1+1):(cu$endYear[i]-y1+1), i] %>% mean()
      sMSYKFLast[i] <-  sMSYKF_mat[(cu$endYear[i]-y1+1-3):(cu$endYear[i]-y1+1),i] %>% mean()
      
      sGenKF_m <- as_tibble ( x = list(ricA=ricAKF_mat[(cu$stYear[i]-y1+1):(cu$endYear[i]-y1+1),i], 
                                       ricB = rep(ricBKF[i],length(ricAKF)), ricSig=1, 
                                       sMSY=sMSYKF_mat[(cu$stYear[i]-y1+1):(cu$endYear[i]-y1+1),i]))
      sGenKF_m <- sGenKF_m %>% transmute(sGen= sGenSolver_v(ricA,ricB,ricSig,sMSY))
      sGenKF_mat[(cu$stYear[i]-y1+1):(cu$endYear[i]-y1+1), i] <- sGenKF_m$sGen
      sGenKFMean[i] <- sGenKF_mat[1:length(ricAKF), i] %>% mean()
      sGenKFLast[i] <- sGenKF_mat[(cu$endYear[i]-y1+1-3):(cu$endYear[i]-y1+1),i] %>% mean()
    }#End of if(output$b<0)
    
    offset <- cu$stYear[i]-min(cu$stYear)
    ppnHigh[i] <- length (which (S$ETS > sMaxKF[i])) / length (S$ETS)
    ppnLow[i] <- length (which (S$ETS < 0.2*(ricAKF_mat[(offset+1):(offset+nYrs[i]),i])/c(ricBKF[i]))) / length (S$ETS)
    hRec <-R$rec[which(S$ETS >sMaxKF[i])] # recruitments for large spawning events
    y<-which(S$ETS >sMaxKF[i])
    rMax <- NA
    for (m in 1:length(y)){
      #rMax[m] <- exp(ricAKF_mat[y[m]+offset,i])/(ricBKF[i]*exp(1)) # exp(a)/b*e
      rMax[m] <- exp(mean(ricAKF_mat[,i],na.rm=T))/(ricBKF[i]*exp(1)) # exp(a)/b*e
    }
    ppnHRec[i] <- length (which (hRec > rMax)) / length (which (S$ETS > sMaxKF[i]))
    
    
  }# End of i stocks
  
  out<- list (nYrs = nYrs, ricALR = ricALR, sMaxLR = sMaxLR, ricSig = ricSig, ricLRAICc = ricLRAICc, 
              sMSYLR = sMSYLR, ricSige = ricSige, ricSigw = ricSigw,
              sGenLR = sGenLR, ppnHigh = ppnHigh, ppnHRec = ppnHRec, ricAKFMean = ricAKFMean,
              ricAKFLast = ricAKFLast, sMaxKF = sMaxKF, ricKFAICc = ricKFAICc, snr = snr, 
              sMSYKFMean = sMSYKFMean, sMSYKFLast = sMSYKFLast, sGenKFMean = sGenKFMean, 
              sGenKFLast = sGenKFLast, ricAKF_mat = ricAKF_mat, sMSYKF_mat = sMSYKF_mat, 
              sGenKF_mat = sGenKF_mat)
  
  table <- data.frame (nYrs = nYrs, ppnHigh = ppnHigh, ppnHRec = ppnHRec, ricALR = ricALR,
                       ricAKFMean = ricAKFMean, sMaxLR = sMaxLR, sMaxKF = sMaxKF, 
                       ricSig = ricSig, ricSige = ricSige, ricSigw = ricSigw, 
                       ricLRAICc = ricLRAICc, ricKFAICc = ricKFAICc, sMSYLR = sMSYLR,
                       sMSYKFLast = sMSYKFLast, sGenLR = sGenLR, sGenKFLast = sGenKFLast)            
  row.names(table)<-cu$stock
  return(out)
#}
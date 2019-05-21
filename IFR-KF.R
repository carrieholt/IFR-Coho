# Code to generate KF Ricker a values for 5 IFR Coho CUs
# Developed by Carrie Holt

library(dplyr) 
source("KFcode.R")
source("Funcs.R")

#GetIFRcoho<-function(){
  data <- as.data.frame(read.table("data/IFC_SR.dat", header=TRUE))
 
  strYr <- data%>%group_by(CUid)%>%summarise(Min=min(Byr))
  endYr <- data%>%group_by(CUid)%>%summarise(Max=max(Byr))
  maxYears <- max (endYr$Max - strYr$Min)+1
  
  
  # cu<-as.tbl(read.csv("data/cuToStock.csv", col.name=c("cuName","stock", 
  #                                                      "stYear", "endYear"),
                      
  ricALR <- ricBLR <- ricLRAICc <- sMaxLR <- sMSYLR <- sGenLR <- ricBKF <- 
  ricKFAICc <- sMaxKF <- ricAKFMean <- ricAKFLast <- snr <- nllKF <- varKF <- 
  sMSYKFMean <- sMSYKFLast <- sGenKFLast <- sGenKFMean <- sGenKFLast <-
  ppnHigh <- ppnLow <- ppnHRec <- nYrs <- ricSig <- ricSige <- ricSigw <- NA
  
  nCUs <- length(unique(data$CUnm))
  
  ricAKF_mat <- ricAKFVar_mat <- sMSYKF_mat <- sGenKF_mat <- 
    matrix(NA, nrow=maxYears, ncol=nCUs)
  
  for(i in 1:nCUs){
    
    R <- S <- NA
    R<-data%>%filter(CUid==i)%>%select(Rec_total)
    S<-data%>%filter(CUid==i)%>%select(Sp)
    
    nYrs[i] <- length(S$Sp)
    out.lm <- lm(log(R$Rec_total/S$Sp)~S$Sp)
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
    if(is.na(ricBLR[i])==TRUE) {initial <- list(-1/max(R$Rec_total, na.rm=T), log(1), log(1), 1, 1, 1, "True")}
    names(initial) <- c("b", "ln.sig.e", "ln.sig.w", "mean.a", "var.a", "Ts", "EstB")
    output <- kf.rw(initial=initial,x=S$Sp,y=log(R$Rec_total/S$Sp))
    
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
      ricAKF_mat[,i] <- ricAKF
      ricAKFVar_mat[,i] <- ricAKFVar
      snr[i] <- output$sig.w/output$sig.e
      nllKF[i] <- output$cum.neg.log.lik
      varKF[i] <- output$sig.e^2 + output$sig.w^2
      ricSige[i] <- output$sig.e
      ricSigw[i] <- output$sig.w
      # sMSYKF_mat[, i] <- 
      #   (1 - gsl::lambert_W0(exp(1 - ricAKF_mat[,i]))) / 
      #   rep(ricBKF[i],length(ricAKF))
      # sMSYKFMean[i] <- sMSYKF_mat[, i] %>% mean()
      # sMSYKFLast[i] <-  sMSYKF_mat[(nYrs-3):(nYrs),i] %>% mean()
      # 
      # sGenKF_m <- as_tibble ( x = list(ricA=ricAKF_mat[,i], 
      #                                  ricB = rep(ricBKF[i],length(ricAKF)), ricSig=1, 
      #                                  sMSY=sMSYKF_mat[,i]))
      # sGenKF_m <- sGenKF_m %>% transmute(sGen= sGenSolver_v(ricA,ricB,ricSig,sMSY))
      # sGenKF_mat[, i] <- sGenKF_m$sGen
      # sGenKFMean[i] <- sGenKF_mat[1:length(ricAKF), i] %>% mean()
      # sGenKFLast[i] <- sGenKF_mat[(nYrs-3):(nYrs),i] %>% mean()
    }#End of if(output$b<0)
    

    
  }# End of i stocks
  
  out<- list (nYrs = nYrs, ricALR = ricALR, sMaxLR = sMaxLR, ricSig = ricSig, ricLRAICc = ricLRAICc, 
              sMSYLR = sMSYLR, ricSige = ricSige, ricSigw = ricSigw,
              sGenLR = sGenLR, ppnHigh = ppnHigh, ppnHRec = ppnHRec, ricAKFMean = ricAKFMean,
              ricAKFLast = ricAKFLast, sMaxKF = sMaxKF, ricKFAICc = ricKFAICc, snr = snr)#, 
              # sMSYKFMean = sMSYKFMean, sMSYKFLast = sMSYKFLast, sGenKFMean = sGenKFMean, 
              # sGenKFLast = sGenKFLast, ricAKF_mat = ricAKF_mat, sMSYKF_mat = sMSYKF_mat, 
              # sGenKF_mat = sGenKF_mat)
  
  table <- data.frame (nYrs = nYrs, ricALR = ricALR,
                       ricAKFMean = ricAKFMean, sMaxLR = sMaxLR, sMaxKF = sMaxKF, 
                       ricSig = ricSig, ricSige = ricSige, ricSigw = ricSigw, 
                       ricLRAICc = ricLRAICc, ricKFAICc = ricKFAICc)#, sMSYLR = sMSYLR,
                       # sMSYKFLast = sMSYKFLast, sGenLR = sGenLR, sGenKFLast = sGenKFLast)            
  row.names(table)<-unique(data$CUnm)
  return(out)
#}
  
  png("IFC_KF.png", width=6, height=8, units="in", res=400)
  par(mfrow=c(3,2), las=1, cex=0.8, mar=c(2, 4, 3, 1) + 0.1)
  #par(mfrow=c(1,1))
  for (i in 1:5){
    plot(x=1998:2013, y=ricAKF_mat[,i], xlab="", ylab="Productivity", type="l", ylim=c(-0.7,2.5))
    mtext(unique(data$CUnm)[i], side=3, line=0.5, cex=0.8)
  }
  dev.off()
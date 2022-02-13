  conto_estimation <- function(input_pvalues,lambda=0.5){
    if (is.null(ncol(input_pvalues)))
      stop("input_pvalues should be a matrix or data frame")
    if (ncol(input_pvalues) !=2)
      stop("inpute_pvalues should have 2 column")
    input_pvalues <- matrix(as.numeric(input_pvalues),nrow=nrow(input_pvalues))
    if (sum(complete.cases(input_pvalues))<nrow(input_pvalues))
      warning("input_pvalues contains NAs to be removed from analysis")
    input_pvalues <- input_pvalues[complete.cases(input_pvalues),]
    if (!is.null(nrow(input_pvalues)) & nrow(input_pvalues)<1)
      stop("input_pvalues doesn't have valid p-values")
    pcut <- seq(0.1,0.8,0.1) 
    frac1 <- rep(0,8)
    frac2 <- rep(0,8)
    frac12<- rep(0,8)
    for (i in 1:8) {
      frac1[i] <- mean(input_pvalues[,1]>=pcut[i])/(1-pcut[i])
      frac2[i] <- mean(input_pvalues[,2]>=pcut[i])/(1-pcut[i]) 
      frac12[i]<- mean(input_pvalues[,2]>=pcut[i] & input_pvalues[,1]>=pcut[i])/(1-pcut[i])^2
    }  
    ## use the median estimates for pi00 ##
    
    alpha00 <- min(frac12[pcut==lambda],1)
    
    ## alpha1 is the proportion of nulls for first p-value 
    ## alpha2 is the proportion of nulls for second p-value 
    
    if (ks.test(input_pvalues[,1],"punif",0,1,alternative="greater")$p>0.05) alpha1 <- 1 else   alpha1 <- min(frac1[pcut==lambda],1)  
    if (ks.test(input_pvalues[,2],"punif",0,1,alternative="greater")$p>0.05) alpha2 <- 1 else   alpha2 <- min(frac2[pcut==lambda],1)
    
    
    if (alpha00==1) {
      alpha01 <- 0
      alpha10 <- 0
      alpha11 <- 0
    } else {    
      if (alpha1==1  & alpha2==1) {
        alpha01 <- 0
        alpha10 <- 0
        alpha11 <- 0
        alpha00 <- 1
      }  
      
      if (alpha1==1  & alpha2!=1) {
        alpha10 <- 0
        alpha11 <- 0
        alpha01 <- alpha1-alpha00
        alpha01 <- max(0,alpha01)
        alpha00 <- 1-alpha01
      }  
      
      if (alpha1!=1  & alpha2==1) {
        alpha01 <- 0
        alpha11 <- 0
        alpha10 <- alpha2-alpha00
        alpha10 <- max(0,alpha10)
        alpha00 <- 1-alpha10
      }  
      
      if (alpha1!=1  & alpha2!=1) {
        alpha10 <- alpha2-alpha00
        alpha10 <- max(0,alpha10)
        alpha01 <- alpha1-alpha00
        alpha01 <- max(0,alpha01)
        
        if ((1-alpha00-alpha01-alpha10)<0) {
          alpha11 <- 0
          alpha10 <- 1- alpha1
          alpha01 <- 1- alpha2
          alpha00 <- 1- alpha10 - alpha01
        }  else {
          alpha11 <-  1-alpha00-alpha01-alpha10
        }  
      }  
    }
    alpha.null <- list(alpha10=alpha10,alpha01=alpha01,alpha00=alpha00,alpha1=alpha1,alpha2=alpha2)
    
    nullprop=alpha.null

    if (is.null(ncol(input_pvalues)))
    stop("input_pvalues should be a matrix or data frame")
    if (ncol(input_pvalues) !=2)
      stop("inpute_pvalues should have 2 column")
    input_pvalues <- matrix(as.numeric(input_pvalues),nrow=nrow(input_pvalues))
    if (sum(complete.cases(input_pvalues))<nrow(input_pvalues))
      warning("input_pvalues contains NAs to be removed from analysis")
    input_pvalues <- input_pvalues[complete.cases(input_pvalues),]
    if (!is.null(nrow(input_pvalues)) & nrow(input_pvalues)<1)
      stop("input_pvalues doesn't have valid p-values")
    
    pmax <- apply(input_pvalues,1,max)
    nmed <- length(pmax)
    efdr11 <- rep(0,nmed)
    efdr10 <- rep(0,nmed)
    efdr01 <- rep(0,nmed)
    
    exact=1
    
    # library(fdrtool)
    nmed  <- nrow(input_pvalues)  
    cdf12 <- input_pvalues
    
    xx1 <- c(0,input_pvalues[order(input_pvalues[,1]),1])
    yy1 <- c(0,seq(1,nmed,by=1)/nmed)
    gfit1<- gcmlcm(xx1,yy1,type="lcm")
    xknots1 <- gfit1$x.knots[-1]
    Fknots1 <- cumsum(diff(gfit1$x.knots)*gfit1$slope.knots)
    
    xx2 <- c(0,input_pvalues[order(input_pvalues[,2]),2])
    yy2 <- c(0,seq(1,nmed,by=1)/nmed)
    gfit2<- gcmlcm(xx2,yy2,type="lcm")
    xknots2 <- gfit2$x.knots[-1]
    Fknots2 <- cumsum(diff(gfit2$x.knots)*gfit2$slope.knots)
    
    if (nullprop$alpha1!=1) Fknots1 <- (Fknots1 - nullprop$alpha1*xknots1)/(1-nullprop$alpha1) else Fknots1 <- rep(0,length(xknots1))
    if (nullprop$alpha2!=1) Fknots2 <- (Fknots2 - nullprop$alpha2*xknots2)/(1-nullprop$alpha2) else Fknots2 <- rep(0,length(xknots2))
    
    
    orderq1 <- pmax
    orderq2 <- pmax
    
    gcdf1 <- pmax
    gcdf2 <- pmax
    for (i in 1:length(xknots1)) {
      if (i==1) {
        gcdf1[orderq1<=xknots1[i]] <- (Fknots1[i]/xknots1[i])*orderq1[orderq1<=xknots1[i]] 
      } else {   
        if (sum(orderq1>xknots1[i-1] & orderq1<=xknots1[i])>0){
          temp <- orderq1[orderq1>xknots1[i-1] & orderq1<=xknots1[i]] 
          gcdf1[orderq1>xknots1[i-1] & orderq1<=xknots1[i]] <- Fknots1[i-1] + (Fknots1[i]-Fknots1[i-1])/(xknots1[i]-xknots1[i-1])*(temp-xknots1[i-1])
        }
      }
    }
    
    for (i in 1:length(xknots2)) {
      if (i==1) {
        gcdf2[orderq2<=xknots2[i]] <- (Fknots2[i]/xknots2[i])*orderq2[orderq2<=xknots2[i]] 
      } else {   
        if (sum(orderq2>xknots2[i-1] & orderq2<=xknots2[i])>0){
          temp <- orderq2[orderq2>xknots2[i-1] & orderq2<=xknots2[i]] 
          gcdf2[orderq2>xknots2[i-1] & orderq2<=xknots2[i]] <- Fknots2[i-1] + (Fknots2[i]-Fknots2[i-1])/(xknots2[i]-xknots2[i-1])*(temp-xknots2[i-1])
        } 
      }
    }
    
    
    gcdf1 <- ifelse(gcdf1>1,1,gcdf1)
    gcdf2 <- ifelse(gcdf2>1,1,gcdf2)
    
    cdf12[,1] <- gcdf1
    cdf12[,2] <- gcdf2
    
    for (i in 1:nmed) {
      fdr11 <-  (pmax[i]*cdf12[i,2]*nullprop$alpha01)/mean(pmax<=pmax[i])
      fdr12 <-  (pmax[i]*cdf12[i,1]*nullprop$alpha10)/mean(pmax<=pmax[i])          
      fdr2  <-  (pmax[i]*pmax[i]*nullprop$alpha00)/mean(pmax<=pmax[i])   
      efdr11[i] <- fdr11+fdr12+fdr2
      efdr10[i] <- fdr12+fdr2
      efdr01[i] <- fdr11+fdr2
      
    }  
    
    efdr11.order <- efdr11[order(pmax,decreasing=T)]
    efdr10.order <- efdr10[order(pmax,decreasing=T)]
    efdr01.order <- efdr01[order(pmax,decreasing=T)]
    
    for (i in 2:nmed)  {
      efdr11.order[i] <- min(efdr11.order[i],efdr11.order[i-1])
      efdr10.order[i] <- min(efdr10.order[i],efdr10.order[i-1])
      efdr01.order[i] <- min(efdr01.order[i],efdr01.order[i-1])
    }
    
    efdr11 <- efdr11.order[rank(-pmax)]
    efdr10 <- efdr10.order[rank(-pmax)]
    efdr01 <- efdr01.order[rank(-pmax)]
    
    
    all=data.frame(input_pvalues,efdr11,efdr10,efdr01)
    return(all)
  }

require(LongMemoryTS)
source("tsGLASSO.R")

###########################################
## The msVAR function heavily relies on 
## the inital code introduced in Davis et. al (2016)
## Most of the functions were taken from the code
## provided in the supplement of the Davis et. al (2016)

#------------------
# matr.sum(matr)
#------------------
matr.sum = function(matr){
    n = nrow(matr)
    m = ncol(matr)
    p = m/n
    result = diag(0,n)
    for(i in 1:p){
        result = result + matr[,(i*n-n+1):(i*n)]
    }
    return(result)
}
#-------------------------------------------------------------------------------------------------------------
# simulateVAR(coefMatr=NULL,inteVect,size,burn=100,Sigma,error.dist=c("normal","t"),df=NULL,checkCausal=TRUE)
# simulate causal VAR(p) process
#-------------------------------------------------------------------------------------------------------------
simulateVAR = function(coefMatr=NULL,intercept,size,burn=100,Sigma,
                         error.dist=c("normal","t"),df=NULL,
                         checkCausal=TRUE){
  
    error.dist  = match.arg(error.dist)
    if(!is.null(coefMatr)) {coefMatr = as.matrix(coefMatr)}
    inteVect    = intercept
    #	inteVect    = matrix(inteVect,ncol=1) Command was in actual code
    K           = nrow(Sigma)
  #  print(K)
    if(error.dist=="t"){
        if(is.null(df)|df<=2) {
          stop("error: the d.f. of t-distribution is NULL or <=2")
          }
    }
    if (min(eigen(Sigma)$val)<0){
      stop("error: the covariance matrix is not positive-definite")
      }
    if (!is.null(coefMatr)){ 
        if(checkCausal & !companionVAR(coefMatr)$isCausal){stop("error: the VAR model is non-causal")}
    }
    S.scale        = diag(NA,K) 
    if(error.dist=="t") {S.scale = Sigma*(df-2)/df}
    if (!is.null(coefMatr)) {p = ncol(coefMatr)/K} else {p=0}
    simData        = matrix(NA,nrow=K,ncol=size+burn)
    cholSigma      = t(chol(Sigma))
    if(p!=0){
        simData[,1:p]  = matrix(rnorm(K*p),nrow=K,ncol=p)
        for (i in (p+1):(size+burn)){
            if (error.dist=="normal") {error = cholSigma %*% matrix(rnorm(K),ncol=1)}
            if (error.dist=="t")      {error = matrix(rmt(n=1,rep(0,K),S.scale,df),ncol=1)}
            simData[,i] = inteVect + coefMatr %*% matrix(c(simData[,(i-1):(i-p)]),ncol=1) + error
        }}
    if(p==0){
        for (i in (p+1):(size+burn)){
            if (error.dist=="normal") {error = cholSigma %*% matrix(rnorm(K),ncol=1)}
            if (error.dist=="t")      {error = matrix(rmt(n=1,rep(0,K),S.scale,df),ncol=1)}
            simData[,i] = inteVect + error
        }}
    simData  = t(simData[,(1+burn):(size+burn)])
    result   = list(simData=simData,coefMatr=coefMatr,inteVect=inteVect,
                    K=K,p=p,size=size,checkCausal=checkCausal,
                    Sigma=Sigma,error.dist=error.dist,df=df)
    return(result)
}


#-------------------------------------------------------------------------------------------------------------
# simulateWN(inteVect,size,burn=100,Sigma,error.dist=c("normal","t"),df=NULL)
# simulate white noise process
#-------------------------------------------------------------------------------------------------------------
simulateWNTG = function(inteVect,size,burn=100,
                        Sigma,error.dist=c("normal","t"),df=NULL){
  
    error.dist  = match.arg(error.dist)
    inteVect    = matrix(inteVect,ncol=1)
    K           = nrow(Sigma)
    if(error.dist=="t"){
        if(is.null(df)|df<=2){
            stop(" error: the d.f. of t-distribution is NULL or <=2")
        }
    }
    if(min(eigen(Sigma)$val)<0){
        stop(" error: the covariance matrix is not positive-definite")
    }
    S.scale        = diag(NA,K) 
    if(error.dist=="t") {
        S.scale = Sigma*(df-2)/df
    }
    simData        = matrix(NA,nrow=K,ncol=size+burn)
    cholSigma      = t(chol(Sigma))
    for (i in 1:(size+burn)){
        if(error.dist=="normal"){
            error = cholSigma %*% matrix(rnorm(K),ncol=1)
        }
        if(error.dist=="t"){
            error = matrix(rmt(n=1,rep(0,K),S.scale,df),ncol=1)
        }
        simData[,i] = inteVect + error
    }
    simData  = t(simData[,(1+burn):(size+burn)])
    result   = list(simData=simData,inteVect=inteVect,K=K,size=size,
                    Sigma=Sigma,error.dist=error.dist,df=df)
    return(result)
}

#------------------------------------------------------------------------------------------------------------------------
# companionVAR(A)
# compute the companion matrix of the VAR model's coefficient matrix A and decide whether or not the VAR model is causal.
#------------------------------------------------------------------------------------------------------------------------
companionVAR = function(A){ 
    A         = as.matrix(A)
    K         = nrow(A)
    p         = ncol(A)/K
    companion = matrix(0, nrow=K*p, ncol=K*p)
    companion[1:K,1:(K*p)] = A
    if (p > 1) {companion[(K+1):(K*p),1:(K*p-K)] = diag(1,K*p-K)}
    eigens    = eigen(companion)$val
    mods      = Mod(eigens)
    result    = list(coefMatr=A, isCausal = max(mods) < 1, 
                     companion=companion, eigens=eigens, mods=mods)
    return(result)
}

#-------------------------------------------------
# isInteger(number)
# check whether or not number is an integer.
#------------------------------------------------
isInteger = function(number){
    return (floor(number)==number)
}

#-------------------------------------------------------------------
# sqrtInv(matr)                                       
# matr is a square positive-definite matrix. Compute matr^{-1/2}
#-------------------------------------------------------------------
sqrtInv = function(matr){
    eigenval = eigen(matr)$val
    eigenvec = eigen(matr)$vec
    if (min(eigenval)<0){
        stop("The matrix is not positive definite.")
    }
    return(eigenvec %*% diag(1/sqrt(eigenval)) %*% t(eigenvec))
}

#------------------------------------------------------------------------
# find.uniquevalue.breakpoint(vec)
# find the breakpoints of the unique values in vec (probabaly with ties.)
#------------------------------------------------------------------------
find.uniquevalue.breakpoint = function(vec){
    unique.value = unique(vec)
    temp         = rep(NA,length(unique.value))
    for (i in 1:length(unique.value)){temp[i] = sum(vec==unique.value[i])}
    return(cumsum(temp))                       
}

#----------------------------------------------------------------------------
# specCompute(dta, spans = NULL, kernel = NULL, taper = 0.1, pad = 0, 
#             fast = TRUE, demean = FALSE, detrend = TRUE, plot = FALSE, 
#             upweight = NULL, na.action = na.fail,...) 
# compute the spectral domain quantities.
#----------------------------------------------------------------------------
specComputeTG <- function(dta, spans = NULL, 
                          kernel = NULL, taper = 0.1, 
                          pad = 0,ADMM_ITER = 50, 
                          lambda,
                          rho = 10, alpha = 1.5,
                          thresh = 1e-6, 
                          halfWindowLength, rho.flex = FALSE, 
                          diag = FALSE,
                        fast = TRUE, demean = FALSE,
                        detrend = TRUE, plot = FALSE, 
                        ABSTOL = 1e-6, RELTOL = 1e-4, 
                        upweight = NULL, na.action = na.fail,...){
    if(ncol(dta)<=2){	
        stop("The series must be of at least 3 dimensions")
    }
    x       <- dta
    series  <- deparse(substitute(x))
    x       <- na.action(as.ts(x))
    xfreq   <- frequency(x)
    x       <- as.matrix(x)
    N       <- nrow(x) 
    N0      <- nrow(x)
    nser    <- ncol(x)
    if (!is.null(spans)){ 
        kernel <- {
            if (is.tskernel(spans)) 
                spans
            else kernel("modified.daniell", spans%/%2)
        }
    }
    if (!is.null(kernel) && !is.tskernel(kernel)){
        stop("must specify 'spans' or a valid kernel")
    }
    if (detrend){
        t     <- 1L:N - (N + 1)/2
        sumt2 <- N * (N^2 - 1)/12
        for (i in 1L:ncol(x)){
            x[, i] <- x[, i] - mean(x[, i]) - sum(x[, i] * t) * t/sumt2
        }
    }
    if (demean){
        x <- sweep(x, 2, colMeans(x), check.margin = FALSE)
    }
    x   <- spec.taper(x, taper)
    u2  <- (1 - (5/8) * taper * 2)
    u4  <- (1 - (93/128) * taper * 2)
    if (pad > 0) {
        x <- rbind(x, matrix(0, nrow = N * pad, ncol = ncol(x)))
        N <- nrow(x)
    }
    if (fast){
        NewN <- nextn(N)
    }
    if (!fast){
        NewN <- N
    }
    x     <- rbind(x, matrix(0, nrow = (NewN - N), ncol = ncol(x)))
    N     <- nrow(x)
    Nspec <- floor(N/2)
    freq  <- seq.int(from = xfreq/N, by = xfreq/N, length.out = Nspec)
    xfft  <- mvfft(x)

    xfft  <- mvfft(x) ## Convert to Fourier 
    pgram = Peri(x)
    kernel <- kernel("modified.daniell",halfWindowLength)
    
    if (!is.null(kernel)){
        for (i in 1L:ncol(x)){ 
            for (j in 1L:ncol(x)){ 
                pgram[i, j, ] <- kernapply(pgram[i, j, ], kernel, circular = TRUE)
            }
        }
        df        <- df.kernel(kernel)
        bandwidth <- bandwidth.kernel(kernel)
    }
    if (is.null(kernel)){
        df        <- 2
        bandwidth <- sqrt(1/12)
    }
    df        <- df/(u4/u2^2)
    df        <- df * (N0/N)
    bandwidth <- bandwidth * xfreq/N
    if(!is.null(upweight)){
        A <- matrix(0,nrow=nser,ncol=nser)
        for(i in 1L:nser){ 
            A[i,i] <- max(Re(pgram[,i,i]))
        }
        for(k in 1:Nspec){
            pgram[,,k] <- pgram[,,k] + upweight * A
        }
    }
    fxx = pgram
    rm(pgram) 
#    gmodel = ADMM_LASSO_TimeSeries(S= fxx, lambda=lambda, rho = rho, alpha = alpha, MAX_ITER = ADMM_ITER,
 #                         thr = thr)
    #### Modification here 11/19/2020
    gmodel = glasso_TimeSeries(dta = dta, lambda = lambda,
                               rho = rho, alpha = alpha,
                             MAX_ITER = ADMM_ITER, thresh = thresh, 
                             ABSTOL   = ABSTOL,
                             RELTOL   = RELTOL, rho.flex = rho.flex,
                             diag = diag, halfWindowLength = halfWindowLength)
    gpre <- gmodel$Z
  #  fr = gmodel$fr
  #  Nspec = length(fr)
    Nspec = dim(gpre)[3] ## Added 11/19/2020
    gpresmry <- matrix(NA, nrow=Nspec, ncol=nser*(nser-1)/2)
    for (k in 1:Nspec){
        for (i in 1:(nser - 1)){
            for (j in (i+1):nser){
                gpresmry[k, i + (j - 1) * (j - 2)/2]   <- gpre[i,j,k]
            }
        }
    }

    
    #-----------------------------  
    spec.output <- list(freq = freq, gpresmry = gpresmry,
                        kernel = kernel, df = df, bandwidth = bandwidth,
                        n.used = N, upweight = upweight,
                        orig.n = N0, series = series, snames = colnames(x), 
                        method = ifelse(!is.null(kernel), "Smoothed Periodogram",
                                        "Raw Periodogram"), 
                        taper = taper, pad = pad, detrend = detrend, 
                        demean = demean, u2=u2, u4=u4)
    class(spec.output) <- "spec"
    return(spec.output)
}

#------------------------------------------------------------------------------------------------------------------------------------------
# sparseVAR(dta,p,p.ub,d,sigPairs,method=c("mle"),startValue=c("uniYW","sparseOLS"),  
#           startIntercept=NULL,startA=NULL,
#           forceAutoRegression=TRUE,showStatus=NULL,
#           iteMax=150,reltol=1e-4,returnConstrMatr = FALSE,computeOneStepMSE=FALSE,computeAsyVar=FALSE)
# stage 1 screening of the 2-stage approach
#------------------------------------------------------------------------------------------------------------------------------------------
sparseVARTG   = function(dta, p, p.ub=NULL, d, sigPairs=NULL,
                         method=c("mle"), startValue = c("uniYW","sparseOLS"),
                       startIntercept=NULL, startA=NULL,
                       nonZeroAR=NULL, Sigma.given=NULL,
                       forceAutoRegression=TRUE, showStatus=NULL,
                       iteMax=iteMax, reltol=1e-4, computeOneStepMSE=FALSE,
                       computeAsyVar=FALSE){
  
    if((is.null(sigPairs)) & (p>0)) {forceAutoRegression=TRUE}						 
    method      = match.arg(method)
    startValue  = match.arg(startValue)
    dimNames    = colnames(dta)
    if(is.null(p.ub)) {p.ub=p}
    N.eff       = nrow(dta)-p.ub
    K           = ncol(dta)
    dta         = t(dta)
    result      = NULL 
    
    
    if(p==0){
      cat("p is zero")
        computeAsyVar = FALSE
        computeOneStepMSE = FALSE 
        para.num.VAR  = 0
        para.num.cov  = K*d - d*(d-1)/2 + ifelse(d<K,1,0)
        para.num      = para.num.cov
        SIGMA.est     = diag(NA,K)
        dta.eff       = NULL
        if(p.ub==0)   {dta.eff=dta}
        if(p.ub>0)    {dta.eff=dta[,-(1:p.ub)]}
        S             = 1/N.eff * dta.eff %*% t(dta.eff)
        S.evd         = eigen(S)
        sigma2.est    = NULL
        if (d==K) {sigma2.est=0}
        if (d<K)  {sigma2.est = mean(S.evd$val[(d+1):K])}
        lambda.est    = S.evd$val[1:d]-sigma2.est
        U.est         = matrix(S.evd$vec[,1:d],ncol=d)
        if (d > 1){SIGMA.est = U.est %*% diag(lambda.est) %*% 
          t(U.est) + sigma2.est * diag(1,K)}
        
        if (d==1) {SIGMA.est = lambda.est * U.est %*% 
          t(U.est) + sigma2.est * diag(1,K)} 
        
        sum.temp      = 0
        SIGMA.est.inv = solve(SIGMA.est)
        for(i in 1:N.eff){ sum.temp  = sum.temp + matrix(dta.eff[,i],nrow=1) %*% 
          SIGMA.est.inv %*% matrix(dta.eff[,i],ncol=1)}
        negLogLike    = 0.5*(K*N.eff*log(2*pi) + N.eff*log(det(SIGMA.est)) + sum.temp)
        aic           = 2*negLogLike + 2 * para.num
        bic           = 2*negLogLike + log(N.eff)*para.num
        result        = list(order=0,order.max=p.ub,d=d,N=N.eff+p.ub,N.eff=N.eff,
                             sigPairs=NULL,method=method,forceAutoRegression=NULL,
                             nonZeroAR=NULL,
                             iteMax=NULL,reltol=NULL,startValue=NULL,
                             traceNegLogLike=NULL,negLogLike=negLogLike,
                             aic=aic,bic=bic, 
                             para.num.VAR=0, para.num.cov=para.num.cov,
                             para.num.total=para.num, residual=NULL,
                             estIntercept=rep(0,K),estA=NULL, 
                             estSigma=SIGMA.est,oneStepMSE=NULL)
    }
    
    if(p>=1){
        Y           = dta[,-(1:p.ub)] 
        Z           = matrix(NA, nrow=K*p, ncol=N.eff)
        if (p==1){Z = dta[,p.ub:(p.ub+N.eff-1)]}        
        if (p > 1){
            for (i in 1:N.eff){
                Z.temp = dta[,(p.ub+i-p):(p.ub+i-1)]
                Z.temp = Z.temp[,p:1]
                Z.temp = c(Z.temp)
                Z[,i] = Z.temp 
            } 
        } 
        
        #CHANGE HERE FOR NULL sigPairs
        Rmatr        = NULL
        para.num.VAR = NULL
        if(is.null(nonZeroAR)){
            temp         = constraintMatr(K, p, sigPairs, forceAutoRegression)
            Rmatr        = temp$constrMatr
            para.num.VAR = temp$M
        }
        if(!is.null(nonZeroAR)){
            nonZeroAR     = matrix(as.numeric(nonZeroAR),nrow=nrow(nonZeroAR))
            M.temp        = sum(nonZeroAR!=0)
            Rmatr         = matrix(0,nrow=K*K*p,ncol=M.temp)
            nonZeroAR.loc = which(c(nonZeroAR)!=0)
            for(i in 1:M.temp){
                Rmatr[nonZeroAR.loc[i],i] = 1
            }
            rm(nonZeroAR.loc)
            para.num.VAR = M.temp
        }
        Rindex       = which(apply(Rmatr,1,function(x){max(x)==1}))
        para.num.cov = K*d - d*(d-1)/2 + ifelse(d<K,1,0)
        para.num     = para.num.VAR + para.num.cov
        B.old        = NULL
        B.new        = NULL
        SIGMA.new    = NULL
        SIGMA.old    = NULL
        
        if (is.null(startA)){
            if(startValue=="uniYW"){
                A0 = matrix(0,nrow=K,ncol=K*p)
                for(i in 1:K){
                    artemp = ar(dta[i,],FALSE,p)   # univariate Yule-Walker is used.
                    A0[i,i+(K*(0:(p-1)))] = artemp$ar
                }
                B.old = A0
                colnames(B.old) = NULL
                rm(A0); rm(artemp); 
            }
        }
        if (!is.null(startA)){
            B.old           = startA
            colnames(B.old) = NULL 	
        }
        
        rslt.seq    = rep(NA,iteMax+2)
        rslt.seq[1] = 10^10
        rslt.seq[2] = 10^8
        ite         = 2
        
        while((ite <= (iteMax+2)) & (abs((rslt.seq[ite]-rslt.seq[ite-1])/
                                         rslt.seq[ite-1]) > reltol)){
            ite    = ite + 1
            # update SIGMA
            S      = 1/N.eff * (Y - B.old%*%Z) %*% t((Y - B.old%*%Z))
            S.evd  = eigen(S)
            sigma2.mle = NULL
            if (d==K) {sigma2.mle=0}
            if (d<K)  {sigma2.mle = mean(S.evd$val[(d+1):K])}
            lambda.mle = S.evd$val[1:d]-sigma2.mle
            U.mle      = matrix(S.evd$vec[,1:d],ncol=d)
            if (d > 1){SIGMA.new = U.mle %*% diag(lambda.mle) %*% t(U.mle) + 
              sigma2.mle * diag(1,K)}
            if (d==1) {SIGMA.new = lambda.mle * U.mle %*% t(U.mle) + 
              sigma2.mle * diag(1,K)} 
            if(!is.null(Sigma.given)){
                SIGMA.new = Sigma.given
            }
            # updata B
            sqrtSIGMAinv = sqrtInv(SIGMA.new)
            Xnew         = kronecker(t(Z),sqrtSIGMAinv)
            #Xnew         = Xnew %*% Rmatr
            Xnew          = as.matrix(Xnew[,Rindex])
            #Ynew         = kronecker(diag(1,ncol(Z)),sqrtSIGMAinv) %*% c(Y)
            Ynew         = c(sqrtSIGMAinv%*%Y)
            gamma.til     = solve(crossprod(Xnew),crossprod(Xnew,Ynew))
            #B.new        = matrix(Rmatr%*%gamma.til,nrow=K)
            beta.til      = matrix(0,nrow=nrow(Rmatr),ncol=ncol(gamma.til))
            beta.til[Rindex,] = gamma.til
            B.new         = matrix(beta.til,nrow=K)
  #          cat("B.new =",B.new,"\n")
            # swap B
            B.old        = B.new
            # compute -log-likelihood
            A            = B.new
            Y0           = dta[,(p.ub+1):(p.ub+N.eff)]
            X            = matrix(NA,nrow=K*p,ncol=N.eff)
            for (i in 1:N.eff){
                X.temp   = dta[,(p.ub+i-1):(p.ub+i-p)]
                X[,i]    = matrix(X.temp,ncol=1) 
            } 
            rslt.seq[ite]= 0.5 * (K * N.eff * log(2*pi) + N.eff * 
                                    log(det(SIGMA.new))+ sum(diag(t(Y0 - A%*%X) %*%
                                    solve(SIGMA.new) %*% (Y0 - A%*%X))))
            if(!is.null(showStatus)){
                if((ite-2)%%showStatus==1){ print(paste("ite= ", ite-2, ", -log(likelihood)= ",
                                                        round(rslt.seq[ite],digits=3),sep=""))}      
            }
        }
        A.est               = B.new
#        cat("A.est =",A.est,"\n")
        rownames(A.est)     = dimNames
        rownames(SIGMA.new) = dimNames
        colnames(SIGMA.new) = dimNames
        residual            = t(Y - B.new%*%Z)
        residual            = rbind(matrix(NA,nrow=p.ub,ncol=K),residual)
        colnames(residual)  = dimNames
        negLogLike          = rslt.seq[ite]
        aic                 = 2*negLogLike + 2 * para.num
        bic                 = 2*negLogLike + log(N.eff) * para.num
        if (computeOneStepMSE){
            GAMMA       = Z %*% t(Z)
            middle.temp = Rmatr %*% solve(t(Rmatr)%*%kronecker(GAMMA,
                                                               solve(SIGMA.new))%*%Rmatr)%*%t(Rmatr)
            OMEGA       = diag(0,K)
            for (i in 1:N.eff){
                Zt         = matrix(Z[,i],ncol=1)
                first.temp = kronecker(t(Zt),diag(K))
                OMEGA      = OMEGA + first.temp %*% middle.temp %*% t(first.temp)
            }
            oneStepMSE = SIGMA.new + (OMEGA/N.eff)
        }
        if (!computeOneStepMSE){ 
            oneStepMSE = NULL 
        }
        if(computeAsyVar){
            X            = matrix(NA,nrow=K*p,ncol=N.eff)
            for (i in 1:N.eff){
                X.temp   = dta[,(p.ub+i-1):(p.ub+i-p)] - matrix(rep(matrix(
                  apply(dta[,-(1:p.ub)],1,mean),ncol=1),p),ncol=p)
                X[,i]    = matrix(X.temp,ncol=1) 
            } 
            asyVarA         = diag(Rmatr%*%solve(t(Rmatr)%*%kronecker(X%*%t(X),
                                          solve(SIGMA.new))%*%Rmatr)%*%t(Rmatr)) 
            tRatioA         = rep(0,K^2*p)
            tRatioA[asyVarA!=0] = c(A.est)[asyVarA!=0]/sqrt(asyVarA[asyVarA!=0])
            asyVarA         = matrix(asyVarA,nrow=K)
            tRatioA         = matrix(tRatioA,nrow=K)
        }
        if(!computeAsyVar){
            asyVarA         = NULL
            tRatioA         = NULL
        }
        result= list(order=p,order.max=p.ub,d=d,N=N.eff+p.ub,
                     N.eff=N.eff,sigPairs=sigPairs,method=method,
                     forceAutoRegression=forceAutoRegression,
                     nonZeroAR=nonZeroAR,
                     iteMax=iteMax, reltol=reltol, startValue=startValue, 
                     traceNegLogLike=rslt.seq[3:ite], negLogLike=negLogLike, 
                     aic=aic, bic=bic, 
                     para.num.VAR=para.num.VAR, para.num.cov=para.num.cov, 
                     para.num.total=para.num, residual=residual,
                     estA=A.est, estSigma=SIGMA.new, oneStepMSE=oneStepMSE,
                     asyVarA=asyVarA,tRatioA=tRatioA)
        
    }
    return(result)
}

#--------------------------------------------------------------------------------------------------
# findSigPairs(dta, permuteNumber, permuteBy=c("temporal","all"), showPermuteNumber=FALSE, 
#              rankPairsBy=c("testStat","pValue"), forceAutoRegression=TRUE,
#              kernel=NULL, halfWindowLength=NULL, spans=NULL, taper=0.1, pad=0, fast=TRUE, demean=TRUE, detrend=FALSE, plot=FALSE,
#              upweight=NULL, na.action=na.fail)
# find the significant pairs based on the test statistics or the p-values.
#--------------------------------------------------------------------------------------------------
findSigPairsTG = function(dta,permuteNumber=NULL,forceAutoRegression=TRUE,
                          ADMM_ITER = 50, lambda = 0.1, rho = 10 , 
                          alpha = 1.5, thr = 1e-5,
                        kernel=NULL, halfWindowLength=NULL, spans = NULL, 
                        taper=0.1, pad=0, fast=TRUE,demean=TRUE,
                        detrend=FALSE,plot=FALSE,
                        upweight=NULL,na.action =na.fail,
                        rho.flex = FALSE, ...){  #permuteBy=c("all","temporal"),showPermuteNumber=FALSE,rankPairsBy=c("testStat","pValue")

      if(is.null(kernel)) {kernel = kernel("modified.daniell",halfWindowLength)}  ## RKim: revised
    if(is.null(halfWindowLength)) {halfWindowLength = kernel$m}
    spec.temp    = specComputeTG(dta=dta, kernel = kernel, ADMM_ITER = ADMM_ITER, 
                                 lambda = lambda, rho = rho, alpha = alpha,
                                 thr=thr, rho.flex = rho.flex, spans = spans, 
                                 taper = taper, pad = pad, fast = fast, 
                                 demean = demean, detrend= detrend, plot = plot,
                                 upweight = upweight, na.action = na.action, 
                                 halfWindowLength = halfWindowLength, ...) ## halfWindowlength added 11/19/2020
    specfreq     = spec.temp$freq
    nser <- ncol(dta)
    gpresmry = spec.temp$gpresmry                                  ## RKim: added from here
    gpre.max           = matrix(NA, nrow=nser*(nser-1)/2,ncol=3)
    colnames(gpre.max) = c("from","to","supPre")
    rowpos               = 0
    for (i in 1:(nser-1)){
        for (j in (i+1):nser){
            rowpos                 = rowpos + 1
            k                      = i+(j-1)*(j-2)/2
            gpre.max[rowpos,1:2] = c(i,j)
            gpre.max[rowpos,3]   = max(abs(gpresmry[,k]))
        }
    }
     gpre.maxsort <- gpre.max[which(gpre.max[,3]>0),]  #gpre.max[order(gpre.max[,3], decreasing=TRUE),]  ## RKim: to here

    result = list(dta=dta,forceAutoRegression=forceAutoRegression,
                  specfreq = specfreq, gpre.max = gpre.maxsort,   ## revised!!!!!
                  kernel=kernel,halfWindowLength=halfWindowLength,
                  spans=spans,taper=taper,pad=pad,fast=fast,demean=demean,
                  detrend=detrend,upweight=upweight,specComputeResult=spec.temp)
    return (result)
}

#---------------------------------------------------------------------------------------
# constraintMatr(K,p,sigPairs,forceAutoRegression)                                 
# construct the matrix of linear constraints on the autoregression coefficient matrices. 
# see of Lutkepohl (1991)
#---------------------------------------------------------------------------------------
constraintMatr = function(K,p,sigPairs,forceAutoRegression){
    result = NULL
    if(is.null(sigPairs)){
        eachA                    = diag(1,K)
        eachP                    = K
        M                        = eachP * p
        Rconstr                  = matrix(0, nrow=K^2*p,ncol=M)
        rowLoc                   = which(c(eachA)==1)
        colLoc                   = 1:eachP
        for (k in 1:p){
            rowLoc.temp = rowLoc + K^2*(k-1)
            colLoc.temp = colLoc + eachP*(k-1)
            for (i in 1:eachP){
                Rconstr[rowLoc.temp[i], colLoc.temp[i]] = 1
            }
        }
        result = list(constrMatr=Rconstr, M=M)
    }
    if(!is.null(sigPairs)){
        sigPairs = as.matrix(sigPairs)
 #       print(sigPairs)
        eachA    = diag(0,K)
        for (k in 1:nrow(sigPairs)){
            i = sigPairs[k,1]
            j = sigPairs[k,2]
            eachA[i,j] = eachA[j,i] = 1
        }
        if(forceAutoRegression){
            eachP = nrow(sigPairs)*2 + K
            diag(eachA) = 1
        }
        if(!forceAutoRegression){
            eachP = nrow(sigPairs)*2-sum(sigPairs[,1]==sigPairs[,2])
        }
        M                        = eachP * p
        Rconstr                  = matrix(0, nrow=K^2*p,ncol=M)
        rowLoc                   = which(c(eachA)==1)
        colLoc                   = 1:eachP
        for (k in 1:p){
            rowLoc.temp = rowLoc + K^2*(k-1) 
            colLoc.temp = colLoc + eachP*(k-1)
            for (i in 1:eachP){
                Rconstr[rowLoc.temp[i], colLoc.temp[i]] = 1
            }
        }

        result = list(constrMatr=Rconstr, M=M)
    }
#    cat("Rconstr =",Rconstr,"\n")
    return(result)
}


#-----------------------------------------------------------------------------------------------------------------------------------------------
# stage 2 refining of the 2-stage approach
#-----------------------------------------------------------------------------------------------------------------------------------------------

sparseVARTG.2ndStage = function(dta, p.ub=NULL, d, 
                                tRatioMatr, arCoefMatr,
                                startValue = c("uniYW","sparseOLS"), 
                                selectBy=c("bic","aic"), 
                                Sigma.given=NULL, fdr.q = 0.1,
                                showStatus=NULL, iteMax=200, stepMax=NULL,
                                reltol=1e-4, computeAsyVar=FALSE){
    
    startValue  = match.arg(startValue)
    dimNames    = colnames(dta)
    tRatioMatr  = as.matrix(tRatioMatr)
    p           = ncol(tRatioMatr)/nrow(tRatioMatr)
    if(is.null(p.ub)) {p.ub=p}
    K           = ncol(dta)
    N.eff       = nrow(dta)-p.ub
    dta         = t(dta)
    result      = NULL
    Y           = dta[,-(1:p.ub)] 
    Z           = matrix(NA, nrow=K*p, ncol=N.eff)
    if (p==1){Z = dta[,p.ub:(p.ub+N.eff-1)]}        
    if (p > 1){
        for (i in 1:N.eff){
            Z.temp = dta[,(p.ub+i-p):(p.ub+i-1)]
            Z.temp = Z.temp[,p:1]
            Z.temp = c(Z.temp)
            Z[,i] = Z.temp 
        } 
    } 
    ## Implementation of FDR
    
    odr <- order(abs(c(tRatioMatr)), decreasing=T)
    absT <- abs(tRatioMatr[odr])
    pval <- 1-pnorm(absT[absT!=0])
    Nbh <- length(pval)
    i <- 1:Nbh
    q1<- fdr.q
    odrtrim <- odr[pval < i/Nbh*q1]
    
    tR <- numeric(length(c(tRatioMatr)))
    tR[odrtrim] <- c(tRatioMatr)[odrtrim]
    tRA <- matrix(tR, nrow=nrow(tRatioMatr))
    
    arCoefMatr = as.matrix(arCoefMatr)
    AR <- numeric(length(c(arCoefMatr)))
    AR[odrtrim] <- c(arCoefMatr)[odrtrim]
    A.est <- matrix(AR, nrow=nrow(arCoefMatr))
    residual            = t(Y - A.est%*%Z)
    residual            = rbind(matrix(NA,nrow=p.ub,ncol=K),residual)
    
    SIGMA.new    = NULL
    B.old <- A.est
    S      = 1/N.eff * (Y - B.old%*%Z) %*% t((Y - B.old%*%Z))
    S.evd  = eigen(S)
    sigma2.mle = NULL
    if (d==K) {sigma2.mle=0}
    if (d<K)  {sigma2.mle = mean(S.evd$val[(d+1):K])}
    lambda.mle = S.evd$val[1:d]-sigma2.mle
    U.mle      = matrix(S.evd$vec[,1:d],ncol=d)
    if (d > 1){SIGMA.new = U.mle %*% diag(lambda.mle) %*% t(U.mle) + sigma2.mle * diag(1,K)}
    if (d==1) {SIGMA.new = lambda.mle * U.mle %*% t(U.mle) + sigma2.mle * diag(1,K)} 
    if(!is.null(Sigma.given)){
      SIGMA.new = Sigma.given
    }
    
    
    rownames(A.est)     = dimNames
    rownames(SIGMA.new) = dimNames
    colnames(SIGMA.new) = dimNames
    colnames(residual)  = dimNames
    p.stage2            = ceiling(max(which(apply(A.est!=0,2,any)))/K)	
    result= list(order=p.stage2,d=d,
                 N=N.eff+p.ub,N.eff=N.eff,
                 residual=residual,
                 estA=A.est,estSigma=SIGMA.new)
    return(result)
}
    
#     step.total          = sum(tRatioMatr!=0) 
#     order.tRatio        = order(c(abs(tRatioMatr)),decreasing=TRUE)
#     if(!is.null(stepMax)){
#         step.total = ifelse(step.total<stepMax,step.total,stepMax)
#     }
#     negLogLike.seq      = rep(NA,step.total)
#     aic.seq             = rep(NA,step.total)
#     bic.seq             = rep(NA,step.total)
#     stat.min            = 0
#     startA              = NULL
#     
#     for(ii in 1:step.total){
#         index.temp      = sort(order.tRatio[1:ii])
#         Rmatr           = matrix(0,nrow=K^2*p,ncol=length(index.temp))
#         for(jj in 1:length(index.temp)){
#             Rmatr[index.temp[jj],jj] = 1
#         }
#         Rindex          = which(apply(Rmatr,1,function(x){max(x)==1}))
#         para.num.VAR.temp = ncol(Rmatr)
#         para.num.cov.temp = K*d - d*(d-1)/2 + ifelse(d<K,1,0)
#         para.num.temp     = para.num.VAR.temp + para.num.cov.temp
#         B.old        = NULL
#         B.new        = NULL
#         SIGMA.new    = NULL
#         SIGMA.old    = NULL
#         
#         if(ii==1){
#             if(startValue=="uniYW"){
#                 A0 = matrix(0,nrow=K,ncol=K*p)
#                 for(i in 1:K){
#                     artemp = ar(dta[i,],FALSE,p)   # univariate Yule-Walker is used.
#                     A0[i,i+(K*(0:(p-1)))] = artemp$ar
#                 }
#                 B.old = A0
#                 colnames(B.old) = NULL
#                 rm(A0); rm(artemp); 
#             }
#         }
#         if(ii > 1){
#             B.old = startA
#             colnames(B.old) = NULL
#         }	
#         
#         rslt.seq    = rep(NA,iteMax+2)
#         rslt.seq[1] = 10^10
#         rslt.seq[2] = 10^8
#         ite         = 2
#         
#         while((ite <= (iteMax+2)) & (abs((rslt.seq[ite]-rslt.seq[ite-1])/
#                                          rslt.seq[ite-1]) > reltol)){
#             ite    = ite + 1
#             # update SIGMA
#             S      = 1/N.eff * (Y - B.old%*%Z) %*% t((Y - B.old%*%Z))
#             S.evd  = eigen(S)
#             sigma2.mle = NULL
#             if (d==K) {sigma2.mle=0}
#             if (d<K)  {sigma2.mle = mean(S.evd$val[(d+1):K])}
#             lambda.mle = S.evd$val[1:d]-sigma2.mle
#             U.mle      = matrix(S.evd$vec[,1:d],ncol=d)
#             if (d > 1){SIGMA.new = U.mle %*% diag(lambda.mle) %*% 
#               t(U.mle) + sigma2.mle * diag(1,K)}
#             if (d==1) {SIGMA.new = lambda.mle * U.mle %*% t(U.mle) +
#               sigma2.mle * diag(1,K)} 
#             if(!is.null(Sigma.given)){
#                 SIGMA.new = Sigma.given
#             }
#             # updata B
#             sqrtSIGMAinv = sqrtInv(SIGMA.new)
#             Xnew         = kronecker(t(Z),sqrtSIGMAinv)
#             #Xnew         = Xnew %*% Rmatr
#             Xnew          = as.matrix(Xnew[,Rindex])
#             #Ynew         = kronecker(diag(1,ncol(Z)),sqrtSIGMAinv) %*% c(Y)
#             Ynew         = c(sqrtSIGMAinv%*%Y)
#             #B.new        = matrix(Rmatr%*%gamma.til,nrow=K)
#             gamma.til     = solve(crossprod(Xnew),crossprod(Xnew,Ynew))
#             beta.til      = matrix(0,nrow=nrow(Rmatr),ncol=ncol(gamma.til))
#             beta.til[Rindex,] = gamma.til
#             B.new         = matrix(beta.til,nrow=K)
#             # swap B
#             B.old        = B.new
#             # compute -log-likelihood
#             A            = B.new
#             Y0           = dta[,(p.ub+1):(p.ub+N.eff)]
#             X            = matrix(NA,nrow=K*p,ncol=N.eff)
#             for (i in 1:N.eff){
#                 X.temp   = dta[,(p.ub+i-1):(p.ub+i-p)]
#                 X[,i]    = matrix(X.temp,ncol=1) 
#             } 
#             rslt.seq[ite]= 0.5 * (K * N.eff * log(2*pi) + N.eff * 
#                                     log(det(SIGMA.new))+ sum(diag(t(Y0 - A%*%X) %*%
#                                     solve(SIGMA.new) %*% (Y0 - A%*%X))))
#             if(showStatus){
#                 print(paste(" ARcoef ",ii,": ite= ", ite-2, ", -log(likelihood)= ", 
#                             round(rslt.seq[ite],digits=3),sep=""))      
#             }
#         }
#         negLogLike.temp     = rslt.seq[ite]
#         aic.temp            = 2*negLogLike.temp + 2 * para.num.temp
#         bic.temp            = 2*negLogLike.temp + log(N.eff) * para.num.temp
#         negLogLike.seq[ii]  = negLogLike.temp
#         aic.seq[ii]         = aic.temp
#         bic.seq[ii]         = bic.temp
#         stat.curr           = ifelse(selectBy=="aic",aic.temp,bic.temp)
#         startA              = B.new
#         if((ii==1) | ((ii>=2)&(stat.curr<stat.min))){
#             stat.min            = stat.curr
#             negLogLike          = negLogLike.temp
#             aic                 = aic.temp
#             bic                 = bic.temp
#             para.num.VAR        = para.num.VAR.temp 
#             para.num.cov        = para.num.cov.temp 
#             para.num.total      = para.num.temp
#             A.est               = B.new
#             residual            = t(Y - B.new%*%Z)
#             residual            = rbind(matrix(NA,nrow=p.ub,ncol=K),residual)
#             if(computeAsyVar){
#                 temp            = diag(Rmatr%*%solve(t(Rmatr)%*%kronecker(Z%*%t(Z),
#                                       solve(SIGMA.new))%*%Rmatr)%*%t(Rmatr))
#                 asyVarA         = temp
#                 tRatioA         = rep(0,K^2*p)
#                 tRatioA[asyVarA!=0] = c(A.est)[asyVarA!=0]/sqrt(asyVarA[asyVarA!=0])
#                 asyVarA         = matrix(asyVarA,nrow=K)
#                 tRatioA         = matrix(tRatioA,nrow=K)
#             }
#             if(!computeAsyVar){
#                 asyVarA         = NULL
#                 tRatioA         = NULL
#             }
#         }		
#     }
#     rownames(A.est)     = dimNames
#     rownames(SIGMA.new) = dimNames
#     colnames(SIGMA.new) = dimNames
#     colnames(residual)  = dimNames
#     aic.seq             = aic.seq[length(aic.seq):1]
#     bic.seq             = bic.seq[length(bic.seq):1]
#     p.stage2            = ceiling(max(which(apply(A.est!=0,2,any)))/K)	
#     result= list(order=p.stage2,d=d,N=N.eff+p.ub,
#                  N.eff=N.eff,method=method,
#                  startValue=startValue, info=selectBy, 
#                  negLogLike.seq=negLogLike.seq,aic.seq=aic.seq,
#                  bic.seq=bic.seq,negLogLike=negLogLike, 
#                  aic=aic, bic=bic,
#                  para.num.VAR=para.num.VAR,
#                  para.num.cov=para.num.cov, 
#                  para.num.total=para.num.temp, 
#                  residual=residual,
#                  estA=A.est,estSigma=SIGMA.new,
#                  asyVarA=asyVarA,tRatioA=tRatioA)
#     return(result)
# }

#-------------------------------------------------------------------
# reduceRankCovMatr(dta,info=c("bic","aic"),d=NULL,d.seq=NULL)
# reduced-rank covariance estimation
#-------------------------------------------------------------------
reduceRankCovMatr = function(dta,info=c("bic","aic"),d=NULL,d.seq=NULL){ 
    info           = match.arg(info)
    N              = nrow(dta)
    K              = ncol(dta)
    if(is.null(d.seq)){
        d.seq = seq(from=0,to=K-1,by=1)
    }
    d.sel             = NULL
    aic.seq           = matrix(NA,nrow=1,ncol=length(d.seq))
    colnames(aic.seq) = d.seq
    bic.seq           = aic.seq
    negLogLike.seq    = aic.seq
    result            = NULL 
    inter.est         = apply(dta,2,mean)
    dta.demean        = dta - matrix(1,nrow=N,ncol=1) %*% matrix(inter.est,nrow=1)
    S                 = 1/N * t(dta.demean) %*% dta.demean
    S.evd             = eigen(S)
    if(is.null(d)){
        for(i in 1:length(d.seq)){
            d.temp          = d.seq[i]
            numPara.cov     = K*d.temp - d.temp*(d.temp-1)/2 + ifelse(d.temp<K,1,0)
            sigma2.est      = ifelse(d.temp<K,mean(S.evd$val[(d.temp+1):K]),0)
            SIGMA.est = NULL
            if(d.temp==0){
                SIGMA.est = sigma2.est*diag(1,K)
            }
            if(d.temp>0){
                lambda.est  = c(S.evd$val[1:d.temp]-sigma2.est)
                U.est       = matrix(S.evd$vec[,1:d.temp],ncol=d.temp)
                SIGMA.est   = U.est%*%diag(lambda.est,d.temp)%*%t(U.est)+sigma2.est*diag(1,K)
            }
            SIGMA.est.inv   = solve(SIGMA.est)
            sum.temp        = 0
            for(ii in 1:N){ 
                sum.temp  = sum.temp + matrix(dta.demean[ii,],nrow=1) %*% 
                  SIGMA.est.inv %*% matrix(dta.demean[ii,],ncol=1)
            }
            negLogLike        = 0.5*(K*N*log(2*pi) + N*log(det(SIGMA.est)) + sum.temp)
            negLogLike.seq[i] = negLogLike
            aic.seq[i]        = 2*negLogLike + 2*numPara.cov
            bic.seq[i]        = 2*negLogLike + log(N)*numPara.cov
        }
        # selection result of d.
        info.seq      = aic.seq
        if(info=="bic"){
            info.seq  = bic.seq
        }
        d.sel = d.seq[order(info.seq)[1]]
    }
    if(!is.null(d)){
        d.sel          = d
        negLogLike.seq = NULL
        aic.seq        = NULL
        bic.seq        = NULL
    }
    numPara.cov   = K*d.sel-d.sel*(d.sel-1)/2+ifelse(d.sel<K,1,0)
    # fit using the selected reduced-rank d.sel
    lambda.est=NULL; U.est= NULL; SIGMA.est = NULL
    sigma2.est    = ifelse(d.sel<K,mean(S.evd$val[(d.sel+1):K]),0)
    if(d.sel==0){
        SIGMA.est = sigma2.est * diag(1,K)
    }
    if(d.sel>0){
        lambda.est    = c(S.evd$val[1:d.sel]-sigma2.est)
        U.est         = matrix(S.evd$vec[,1:d.sel],ncol=d.sel)
        SIGMA.est     = U.est%*%diag(lambda.est,d.sel)%*%t(U.est)+sigma2.est*diag(1,K)
    }
    SIGMA.est.inv = solve(SIGMA.est)
    sum.temp      = 0
    for(ii in 1:N){ 
        sum.temp  = sum.temp + matrix(dta.demean[ii,],nrow=1) %*% 
          SIGMA.est.inv %*% matrix(dta.demean[ii,],ncol=1)
    }
    negLogLike  = 0.5*(K*N*log(2*pi) + N*log(det(SIGMA.est)) + sum.temp)
    aic         = 2*negLogLike + 2*numPara.cov
    bic         = 2*negLogLike + log(N)*numPara.cov
    
    mean.est    = inter.est
    factor.est  = NULL
    residual    = dta.demean
    if(d.sel>0){
        factor.est  = dta.demean %*% U.est
        residual    = dta.demean - factor.est %*% t(U.est)
    }
    result      = list(data=dta,info=info,negLogLike=negLogLike,
                       aic=aic,bic=bic,d=d.sel,d.seq=d.seq,
                       negLogLike.seq=negLogLike.seq,
                       aic.seq=aic.seq, bic.seq=bic.seq, 
                       numPara.cov=numPara.cov,estIntercept=mean.est,
                       estSigma=SIGMA.est,estU=U.est,estFactor=factor.est,
                       estLambda=lambda.est,estsigma2=sigma2.est,residual=residual)
    return(result)
}
#-------------------------------------------------------------------
# computeNegLogLike_reduceRankCovMatr(dta,Sigma)
# compute negLogLike of the reduced-rank covariance model
#-------------------------------------------------------------------
computeNegLogLike_reduceRankCovMatr = function(dta,Sigma){
    N              = nrow(dta)
    K              = ncol(dta)
    inter.est      = apply(dta,2,mean)
    dta.demean     = dta - matrix(1,nrow=N,ncol=1) %*% matrix(inter.est,nrow=1)
    Sigma.inv      = solve(Sigma)
    sum.temp       = 0
    for(ii in 1:N){ 
        sum.temp  = sum.temp + matrix(dta.demean[ii,],nrow=1) %*% 
          Sigma.inv %*% matrix(dta.demean[ii,],ncol=1)
    }
    negLogLike        = 0.5*(K*N*log(2*pi) + N*log(det(Sigma)) + sum.temp)
    negLogLike        = as.numeric(negLogLike)
    result      = list(data=dta,negLogLike=negLogLike,Sigma=Sigma)
    return(result)
}

#-------------------------------------------------------------------
# specCompute.fromAR(coefMatr,Sigma,names=NULL,freq.seq.length=NULL)
# compute coh and parcoh from AR coefficients and Sigma_{Z}
#-------------------------------------------------------------------
specCompute.fromAR = function(coefMatr,Sigma,names=NULL,freq.seq.length=NULL){
    if(is.null(freq.seq.length)){
        freq.seq.length = 100
    }
    if(is.null(names)){
        names = 1:nrow(coefMatr)
    }
    coefMatr = (-1)*coefMatr
    K        = nrow(coefMatr)
    p        = ncol(coefMatr)/nrow(coefMatr)
    A.array  = array(NA,dim=c(K,K,p))
    for(i in 1:p){
        A.array[,,i] = coefMatr[,(i*K-K+1):(i*K)]
    }
    Omega          = solve(Sigma)
    freq.seq       = seq(0,pi,length=freq.seq.length)
    inv.spec.array = array(NA,dim=c(K,K,length(freq.seq)))
    spec.array     = array(NA,dim=c(K,K,length(freq.seq)))
    coh            = matrix(NA,nrow=length(freq.seq),
                            ncol=K*(K-1)/2,
                            dimnames=list(NULL,1:(K*(K-1)/2)))
    parcoh         = matrix(NA,nrow=length(freq.seq),
                            ncol=K*(K-1)/2,
                            dimnames=list(NULL,1:(K*(K-1)/2)))
    for(ii in 1:length(freq.seq)){
        lambda = freq.seq[ii]
        re.left = diag(1,K)
        im.left = diag(0,K)
        for(k in 1:p){
            re.left = re.left + t(A.array[,,k])*cos(k*lambda)
            im.left = im.left + t(A.array[,,k])*sin(k*lambda)
        }
        re.right = t(re.left)
        im.right = (-1)*t(im.left)
        left     = matrix(complex(real=re.left,imaginary=im.left),nrow=K)
        right    = matrix(complex(real=re.right,imaginary=im.right),nrow=K)
        inv.spec.array[,,ii] = left %*% Omega %*% right
        spec.array[,,ii]     = solve(inv.spec.array[,,ii])
        norm.inv.spec = diag(1/sqrt(Re(diag(inv.spec.array[,,ii])))) %*% 
          inv.spec.array[,,ii] %*% diag(1/sqrt(Re(diag(inv.spec.array[,,ii]))))
        norm.spec     = diag(1/sqrt(Re(diag(spec.array[,,ii]))))     %*% 
          spec.array[,,ii]     %*% diag(1/sqrt(Re(diag(spec.array[,,ii]))))
      for (i in 1:(K-1)){
           for (j in (i+1):K){
              parcoh[ii, i+(j-1)*(j-2)/2]       = Mod(norm.inv.spec[i,j])^2
             colnames(parcoh)[i+(j-1)*(j-2)/2] = paste(names[i],"vs",names[j],sep=" ")
            coh[ii, i+(j-1)*(j-2)/2]          = Mod(norm.spec[i,j])^2
           colnames(coh)[i+(j-1)*(j-2)/2]    = paste(names[i],"vs",names[j],sep=" ")
      }
  }
    }
    result = list(freq.seq=freq.seq,parcoh=parcoh,
                  coh=coh,inv.spec=inv.spec.array,
                  spec=spec.array)
}

#-----------------------------------------------------------------------------------------------------------------------------------------
# sparseVAR.2stage = function(dta,p.seq=NULL,m.seq=NULL,p.ub=NULL,d=NULL,allPairs=NULL,
#                             startValueMethod = c("uniYW","sparseOLS"),standardize=FALSE, nonZeroAR=NULL, Sigma.given=NULL,
#						      forceAutoRegression=TRUE,stage1.showStatus=FALSE,stage2.showStatus=FALSE, 
#                             iteMax=200,reltol=1e-5,stage1.info=c("bic","aic"),          
#                             stage2.info=c("bic","aic"),halfWindowLength=NULL,stage1.permuteBy="all", stage1.useStat=c("parcoh","coh"),
#                             stage1.rankPairsBy=c("testStat","pValue"),stage1.permuteNumber=NULL,stage1.showPermuteNumber=FALSE)
# 2-stage approach to fitting the sparse VAR model.
#------------------------------------------------------------------------------------------------------------------------------------------
msVAR = function(dta,
                 p.seq=NULL,
                 p.ub=NULL,
                 d=NULL,
                 allPairs=NULL,
                 startValueMethod = c("uniYW","sparseOLS"),
                 standardize=FALSE, 
                 nonZeroAR=NULL,
                 Sigma.given=NULL, 
                 forceAutoRegression=TRUE,
                 stage1.showStatus=FALSE,
                 stage2.showStatus=FALSE, 
                 iteMax=200,
                 reltol=1e-5, 
                 stage1.info=c("bic","aic"), 
                 halfWindowLength=NULL ,
                 ADMM_ITER=30,
                 lambda= NULL,
                 rho=100,
                 alpha=1, 
                 rho.flex = FALSE,
                 fdr.q = 0.1,
                 thr=0.005,...){ 
    K          = ncol(dta)
    mean.vec   = matrix(apply(dta,2,mean),nrow=1)
    dta.demean = as.matrix(dta - matrix(1,nrow=nrow(dta),ncol=1) %*%
                             matrix(apply(dta,2,mean),nrow=1))
    if(is.null(p.seq)){
        p.seq = 0:(floor(nrow(dta)/K)-1)
    }
    if(is.null(d)){	
        d = K
    }
    if(is.null(p.ub)){
        p.ub           = max(p.seq)
    }
    if (is.null(lambda)) {lambda = 2*sqrt(log(p.ub)/nrow(dta))}
    aic.matr           = matrix(NA,nrow=length(p.seq)) ### matrix(NA,nrow=length(p.seq),ncol=length(m.seq))  AD
    rownames(aic.matr) = p.seq
 #   colnames(aic.matr) = m.seq
    bic.matr           = aic.matr
    negLogLike.matr    = aic.matr
    D                  = diag(1,K)
    if(standardize){
        for(i in 1:K){
            D[i,i] = 1/sd(dta.demean[,i]) 
        }
    }
    D.inv             = solve(D)
    dta.std           = dta.demean %*% D
    sigPairsResult    = NULL
    p.stage1          = NULL
    sVAR.stage1.result= NULL
    if(is.null(nonZeroAR)){
        if(is.null(allPairs)){
            if(is.null(halfWindowLength)){
                halfWindowLength = floor(sqrt(nrow(dta)))
            }
            sigPairsResult    = findSigPairsTG(dta=dta.std,
                                               kernel=kernel("modified.daniell",
                                              halfWindowLength), 
                                              rho.flex = rho.flex,  ## RKim: revised
                                           fast=TRUE,
                                           demean=FALSE,
                                           detrend=FALSE,
                                           ADMM_ITER=ADMM_ITER,
                                           lambda= lambda,
                                           rho=rho,alpha=alpha,thr=thr
                                             ) #  permuteBy=stage1.permuteBy,rankPairsBy=stage1.rankPairsBy,permuteNumber=stage1.permuteNumber,showPermuteNumber=stage1.showPermuteNumber
            if (dim(sigPairsResult$gpre.max)[2] <= 1)
            {
              allPairs = matrix(c(1,1), 1, 2)
            }else{
            allPairs = matrix(sigPairsResult$gpre.max[,c(1,2)],ncol=2)		     ## RKim: added	
            }
        }
        sigPairs =  matrix(allPairs,ncol=2)              ## AD
        for (i in 1:length(p.seq)){
            p                  = p.seq[i]
            startA             = NULL

                sVAR = sparseVARTG(dta=dta.std,p=p,p.ub=p.ub,d=d,
                                 sigPairs=sigPairs,Sigma.given=Sigma.given,
                                 startA=startA,
                                 startValue=startValueMethod,
                                 forceAutoRegression=forceAutoRegression,
                                 iteMax=iteMax,reltol=reltol,
                                 computeOneStepMSE=FALSE,computeAsyVar=FALSE)
                negLogLike.matr[i,] = sVAR$negLogLike  #negLogLike.matr[i,j]  AD
                aic.matr[i,]  = sVAR$aic   #aic.matr[i,j] AD
                bic.matr[i,]  = sVAR$bic  #bic.matr[i,j] AD
                startA         = sVAR$estA
 #           }
        }
        temp = NULL
        if(stage1.info=="aic"){
            temp = matrix(which(aic.matr==min(aic.matr),arr.ind=TRUE),ncol=2)
        }
        if(stage1.info=="bic"){
            temp = matrix(which(bic.matr==min(bic.matr),arr.ind=TRUE),ncol=2)
        }
        if(nrow(temp)>1){
            temp = temp[1,]
        }
        p.stage1    = p.seq[temp[1]]
        sVAR.stage1.result = sparseVARTG(dta=dta.std,p=p.stage1,p.ub=p.ub,d=d,
                                       sigPairs=sigPairs,Sigma.given=Sigma.given,
                                       startValue=startValueMethod,
                                       forceAutoRegression=forceAutoRegression,
                                       iteMax=iteMax,reltol=reltol,computeOneStepMSE=FALSE,
                                       computeAsyVar=TRUE)
        sVAR.stage1.result$negLogLike.matr = negLogLike.matr
        sVAR.stage1.result$aic.matr        = aic.matr
        sVAR.stage1.result$bic.matr        = bic.matr
        sVAR.stage1.result$info            = stage1.info
        sVAR.stage1.result$D               = D
    }
    
    if(!is.null(nonZeroAR)){
        p.stage1    = floor(ncol(nonZeroAR)/nrow(nonZeroAR))
        sVAR.stage1.result = sparseVARTG(dta=dta.std,p=p.stage1,
                                         p.ub=p.ub,d=d,
                                       nonZeroAR=nonZeroAR,
                                       Sigma.given=Sigma.given,
                                       startValue=startValueMethod,
                                       forceAutoRegression=forceAutoRegression,
                                       iteMax=iteMax,reltol=reltol,
                                       computeOneStepMSE=FALSE,computeAsyVar=TRUE)
        sVAR.stage1.result$negLogLike.matr = NULL
        sVAR.stage1.result$aic.matr        = NULL
        sVAR.stage1.result$bic.matr        = NULL
        sVAR.stage1.result$info            = NULL
        sVAR.stage1.result$D               = D
    }
    
    sVAR.stage2.result = NULL
    if(p.stage1>0){
        if(stage2.showStatus){
            print(" stage 2 refining... ")
        }
        sVAR.stage2.result = sparseVARTG.2ndStage(dta=dta.std, p.ub=p.ub,
                                                  d=d, tRatioMatr=sVAR.stage1.result$tRatioA,
                                                  Sigma.given=Sigma.given, fdr.q = fdr.q,
                                                  arCoefMatr=sVAR.stage1.result$estA,
                                                  startValue = startValueMethod, 
                                                  showStatus=stage2.showStatus,
                                                  iteMax=iteMax,reltol=reltol,
                                                  computeAsyVar=TRUE)
        
        
        sVAR.stage1.result$residual      = sVAR.stage1.result$residual %*% D.inv
        sVAR.stage1.result$estIntercept  = mean.vec
        temp1                            = sVAR.stage1.result$estA
        for(i in 1:p.stage1){
            temp1[,(i*K-K+1):(i*K)] = D.inv %*% temp1[,(i*K-K+1):(i*K)] %*% D
        }
        sVAR.stage1.result$estA          = temp1; rm(temp1)
        sVAR.stage1.result$estSigma      = D.inv %*% sVAR.stage1.result$estSigma %*% D.inv
        
        sVAR.stage2.result$residual      = sVAR.stage2.result$residual %*% D.inv
        sVAR.stage2.result$estIntercept  = mean.vec
        temp2                            = sVAR.stage2.result$estA
        for(i in 1:(ncol(temp2)/K)){
            temp2[,(i*K-K+1):(i*K)]  = D.inv %*% temp2[,(i*K-K+1):(i*K)] %*% D
        }
        sVAR.stage2.result$estA          = temp2; rm(temp2)
        sVAR.stage2.result$estSigma      = D.inv %*% sVAR.stage2.result$estSigma %*% D.inv
        sVAR.stage2.result$D             = D
    }
    result = list(spectral=sigPairsResult,
                  stage1=sVAR.stage1.result,
                  stage2=sVAR.stage2.result)
    return(result)
}

#--------------------------------------------------------------------------------
# fullVAR =(dta,p.seq=NULL,p.ub=NULL,d=NULL,standardize=FALSE,Sigma.given=NULL,
#		    trace=FALSE,info=c("bic","aic"),computeAsyVar=TRUE)
#--------------------------------------------------------------------------------
fullVAR = function(dta,p.seq=NULL,p.ub=NULL,d=NULL,
                   standardize=FALSE,Sigma.given=NULL,
                   trace=FALSE,info=c("bic","aic"),
                   computeAsyVar=TRUE){
    K          = ncol(dta)
    mean.vec   = matrix(apply(dta,2,mean),nrow=1)
    dta.demean = as.matrix(dta - matrix(1,nrow=nrow(dta),ncol=1) %*%
                             matrix(apply(dta,2,mean),nrow=1))
    if(is.null(p.seq)){
        p.seq = 0:(floor(nrow(dta)/K)-1)
    }
    if(is.null(d)){	
        d = K
    }
    if(is.null(p.ub)){
        p.ub           = max(p.seq)
    }
    aic.seq            = matrix(NA,nrow=1,ncol=length(p.seq))
    colnames(aic.seq)  = p.seq
    bic.seq            = aic.seq
    negLogLike.seq     = aic.seq
    D                  = diag(1,K)
    N.eff              = nrow(dta) - p.ub
    if(standardize){
        for(i in 1:K){
            D[i,i] = 1/sd(dta.demean[,i]) 
        }
    }
    D.inv             = solve(D)
    dta.std           = dta.demean %*% D
    mean.vec          = mean.vec   %*% D
    result            = NULL
    info.min          = NULL
    for(i in 1:length(p.seq)){
        p = p.seq[i]
        if(trace){
            cat(" p=",p,"\n")
        }
        result.temp = fullVAR.mle(dta=dta.std,p=p,p.ub=p.ub,d=d,
                                  computeAsyVar=computeAsyVar,
                                  Sigma.given=Sigma.given)
        aic.seq[i]        = result.temp$aic
        bic.seq[i]        = result.temp$bic
        negLogLike.seq[i] = result.temp$negLogLike
        if(i==1){
            info.min  = ifelse(info=="aic",
                               result.temp$aic,result.temp$bic)
            result    = result.temp
        } else {
            info.curr = ifelse(info=="aic",
                               result.temp$aic,result.temp$bic)
            if (info.curr < info.min){
                info.min  = info.curr
                result  = result.temp
            }
        }
    }
    result$estIntercept    = mean.vec
    result$D               = D
    result$aic.seq         = aic.seq
    result$bic.seq         = bic.seq
    result$negLogLike.seq  = negLogLike.seq
    result$info            = info
    result$residual        = result$residual %*% D.inv
    result$estSigma        = D.inv %*% result$estSigma %*% D.inv
    if(result$order > 0){
        temp                 = result$estA
        for(i in 1:(ncol(temp)/K)){
            temp[,(i*K-K+1):(i*K)] = D.inv %*% temp[,(i*K-K+1):(i*K)] %*% D
        }
        result$estA          = temp
    }
    return(result)	
}


#------------------------------------------------------
# fullVAR.mle(dta,p,p.ub,d=NULL,computeAsyVar=TRUE,Sigma.given=NULL)
#------------------------------------------------------
fullVAR.mle = function(dta,p,p.ub,d=NULL,computeAsyVar=TRUE,Sigma.given=NULL){
    K             = ncol(dta)
    sigPairs.full = NULL
    for(i in 1:(K-1)){
        for(j in (i+1):K){
            sigPairs.full = rbind(sigPairs.full,c(i,j))
        }
    }
    sigPairs.full = as.matrix(sigPairs.full)
    result.trash  = sparseVAR(dta=dta,p=p,p.ub=p.ub,d=d,
                              sigPairs=sigPairs.full,
                              Sigma.given=Sigma.given,
                              computeAsyVar=computeAsyVar)
    result               = list()
    result$order         = result.trash$order
    result$d             = result.trash$d
    result$negLogLike    = result.trash$negLogLike
    result$aic           = result.trash$aic
    result$bic           = result.trash$bic
    result$para.num.VAR  = result.trash$para.num.VAR
    result$para.num.cov  = result.trash$para.num.cov
    result$para.num.total= result.trash$para.num.total
    result$residual      = result.trash$residual
    result$estA          = result.trash$estA
    result$estSigma      = result.trash$estSigma
    result$asyVarA       = result.trash$asyVarA
    result$tRatioA       = result.trash$tRatioA
    return(result)
}

#---------------------------------------------------------
# fullVAR.mle.old(dta,p,p.ub,d,computeAsyVar,Sigma.given)
#---------------------------------------------------------
fullVAR.mle.old = function(dta,p,p.ub,d,computeAsyVar,Sigma.given){
    dta      = as.matrix(dta)
    N.eff    = nrow(dta)-p.ub
    K        = ncol(dta)
    dta      = t(dta)
    Y        = dta[,-(1:p.ub)] 
    if(p>0){
        Z        = matrix(NA, nrow=K*p, ncol=N.eff)
        if (p==1){
            Z = dta[,p.ub:(p.ub+N.eff-1)]
        }        
        if (p > 1){
            for (i in 1:N.eff){
                Z.temp = dta[,(p.ub+i-p):(p.ub+i-1)]
                Z.temp = Z.temp[,p:1]
                Z.temp = matrix(Z.temp,ncol=1)
                Z[,i] = Z.temp 
            } 
        }
    }		
    if(!is.null(Sigma.given)){
        d = K-1-sum(diff(eigen(Sigma.given)$val) < 1e-10)
    }
    para.num.VAR = K^2*p
    para.num.cov = K*d - d*(d-1)/2 + ifelse(d<K,1,0)
    para.num     = para.num.VAR + para.num.cov
    A.est        = NULL
    Sigma.est    = NULL
    
    # compute B
    if(p > 0){
        A.est = matrix(solve(crossprod(t(Z)),crossprod(t(Z),t(Y))),nrow=K)
    } else {
        A.est = NULL
    }
    # compute SIGMA.
    if(is.null(Sigma.given)){
        if(p>0){
            S     = 1/N.eff * (Y-A.est%*%Z) %*% t(Y-A.est%*%Z)
        } else {
            S     = 1/N.eff * Y %*% t(Y)
        }
        S.evd      = eigen(S)
        sigma2.est = ifelse(d==K,0,mean(S.evd$val[(d+1):K]))
        U.est      = matrix(S.evd$vec[,1:d],ncol=d)
        if(d==0){
            Sigma.est = sigma2.est*diag(1,K)
        }
        if(d>0){
            lambda.est  = c(S.evd$val[1:d]-sigma2.est)
            U.est       = matrix(S.evd$vec[,1:d],ncol=d)
            Sigma.est   = U.est%*%diag(lambda.est,d)%*%t(U.est)+sigma2.est*diag(1,K)
        }
    } else {
        Sigma.est = Sigma.given
    }	
    # compute -log-likelihood
    Y0           = dta[,(p.ub+1):(p.ub+N.eff)]
    if(p>0){
        X            = matrix(NA,nrow=K*p,ncol=N.eff)
        for (i in 1:N.eff){
            X.temp   = dta[,(p.ub+i-1):(p.ub+i-p)]
            X[,i]    = matrix(X.temp,ncol=1) 
        }
        negLogLike          = 0.5 * (K * N.eff * log(2*pi) + 
                                       N.eff * log(det(Sigma.est)) + 
                                       sum(diag(t(Y0 - A.est%*%X) %*%
                                                  solve(Sigma.est) %*% 
                                                  (Y0 - A.est%*%X))))
        residual            = t(Y - A.est%*%Z)
        residual            = rbind(matrix(NA,nrow=p.ub,ncol=K),residual)
    } else {
        negLogLike          = 0.5 * (K * N.eff * log(2*pi) + N.eff * 
                                       log(det(Sigma.est)) + 
                                       sum(diag(t(Y0) %*% 
                                                  solve(Sigma.est) %*% Y0)))
        residual            = t(Y)
        residual            = rbind(matrix(NA,nrow=p.ub,ncol=K),residual)
    }
    aic                     = 2*negLogLike + 2 * para.num
    bic                     = 2*negLogLike + log(N.eff) * para.num
    # compute asymptotic variance
    if(computeAsyVar & p>0){
        X            = matrix(NA,nrow=K*p,ncol=N.eff)
        for (i in 1:N.eff){
            X.temp   = dta[,(p.ub+i-1):(p.ub+i-p)] -
              matrix(rep(matrix(apply(dta[,-(1:p.ub)],1,mean),
                                ncol=1),p),ncol=p)
            X[,i]    = matrix(X.temp,ncol=1) 
        } 
        asyVarA         = diag(kronecker(solve(X%*%t(X)),Sigma.est) )
        tRatioA         = rep(0,K^2*p)
        tRatioA[asyVarA!=0] = c(A.est)[asyVarA!=0]/sqrt(asyVarA[asyVarA!=0])
        asyVarA         = matrix(asyVarA,nrow=K)
        tRatioA         = matrix(tRatioA,nrow=K)
    } else {
        asyVarA         = NULL
        tRatioA         = NULL
    }
    result= list(order=p,d=d,negLogLike=negLogLike, aic=aic, bic=bic, 
                 para.num.VAR=para.num.VAR, para.num.cov=para.num.cov,
                 para.num.total=para.num, residual=residual,
                 estA=A.est, estSigma=Sigma.est,
                 asyVarA=asyVarA,tRatioA=tRatioA)
    return(result)
} 

#--------------------------------------------------------------
# forecastVAR(data.train,data.test,estIntercept,estA,h.max)							
#--------------------------------------------------------------
forecastVAR = function(data.train,data.test=NULL,
                       estIntercept=NULL,estA,h.max){
  
    K              = ncol(data.train)
    if(!is.null(data.test)){
        data.test  = matrix(data.test,ncol=K)
    }
    if(is.null(estIntercept)){
        estIntercept = matrix(0,nrow=K,ncol=1)
    } else {
        estIntercept = matrix(estIntercept,nrow=K,ncol=1)
    }
    p              = floor(ncol(estA)/K)
    forecast       = matrix(NA,nrow=p+h.max,ncol=K)
    forecast[1:p,] = matrix(data.train[(nrow(data.train)-p+1):nrow(data.train),])
    for(i in (p+1):(p+h.max)){
        temp = matrix(0,nrow=K,ncol=1)
        for(j in 1:p){
            temp = temp + estA[,(j*K-K+1):(j*K)] %*% matrix(forecast[i-j,],ncol=1)
        }
        temp = temp + matrix(estIntercept,ncol=1)
        forecast[i,] = matrix(temp,nrow=1); rm(temp)
    }
    forecast           = matrix(forecast[-(1:p),],ncol=K)
    rownames(forecast) = 1:h.max
    forecast.error = NULL
    if(!is.null(data.test)){
        forecast.error           = matrix(data.test[1:h.max,]-forecast,ncol=K)
        rownames(forecast.error) = 1:h.max		
    }
    result = list(forecast=forecast,h.max=h.max,
                  forecast.error=forecast.error)
    return(result)
}						

#---------------------
# dnormal(x,mu,Sigma)
#---------------------
dnormal = function(x,mu,Sigma,return.log=FALSE){
    if(length(x)==1){
        if(return.log){
            return(log(dnorm(x=x,mean=mu,sd=sqrt(Sigma))))
        }
        if(!return.log){
            return(dnorm(x=x,mean=mu,sd=sqrt(Sigma)))
        }		
    }
    if(length(x)>1){
        x  = matrix(x,ncol=1)
        mu = matrix(mu,ncol=1)
        Sigma.inv = solve(Sigma)
        K  = nrow(Sigma)
        temp = (-0.5)*(K*log(2*pi)+log(det(Sigma))+
                         t(x-mu)%*%Sigma.inv%*%(x-mu))
        if(return.log){
            return(as.numeric(temp))
        }
        if(!return.log){
            return(as.numeric(exp(temp)))
        }
    }
}							

#--------------------------------
# score.pred.old(dta,estA,estSigma)
#--------------------------------
score.pred.old = function(dta,estA,estSigma){
    K = ncol(dta)
    p = ncol(estA)/K
    log.score = rep(NA,nrow(dta)-p)
    quad.score = rep(NA,nrow(dta)-p)
    for(i in (p+1):nrow(dta)){
        mu.temp         = matrix(0,nrow=K,ncol=1)
        for(k in 1:p){
            mu.temp = mu.temp + estA[,(k*K-K+1):(k*K)]%*%matrix(dta[i-k,],ncol=1)
        }
        log.score[i-p]   = -dnormal(x=dta[i,],mu=mu.temp,
                                    Sigma=estSigma,return.log=TRUE)
        quad.score[i-p]  = -2*dnormal(x=dta[i,],mu=mu.temp,Sigma=estSigma) + 
          1/(sqrt(det(estSigma))*(2*sqrt(pi))^K)
    }
    result = list(log.score=log.score,
                  quad.score=quad.score)
}		

#--------------------------------------
# score.pred(dta,h,estA,estSigma,)
#--------------------------------------
score.pred = function(dta,h,estA,estSigma){
    K  = ncol(dta)
    #dta = dta - matrix(1,nrow=nrow(dta),ncol=1)%*%apply(dta,2,mean)
    A1 = estA[,1:K]
    A2 = estA[,(K+1):(2*K)]
    tmp1 = NULL; tmp2 = NULL; tmpSigma = NULL;
    if(h==1){
        tmp1 = A1; tmp2 = A2; 
        tmpSigma = estSigma
    }
    if(h==2){
        tmp1 = A1%*%A1+A2; tmp2 = A1%*%A2; 
        tmpSigma = A1%*%estSigma%*%t(A1)+estSigma
    }
    if(h==3){
        tmp1 = A1%*%A1%*%A1+A1%*%A2+A2%*%A1; 
        tmp2 = (A1%*%A1+A2)%*%A2; 
        tmpSigma = (A1%*%A1+A2)%*%estSigma%*%t(A1%*%A1+A2)+
          A1%*%estSigma%*%t(A1)+estSigma
    }
    if(h==4){
        tmp1 = (A1%*%A1+A2)%*%(A1%*%A1+A2)+A1%*%A2%*%A1; 
        tmp2= (A1%*%A1+A2)%*%A1%*%A2+A1%*%A2%*%A2; 
        tmpSigma = (A1%*%A1%*%A1+A1%*%A2+A2%*%A1)%*%
          estSigma%*%t(A1%*%A1%*%A1+A1%*%A2+A2%*%A1)+(A1%*%A1+A2)%*%
          estSigma%*%t(A1%*%A1+A2)+A1%*%
          estSigma%*%t(A1)+estSigma
    }
    log.score = rep(NA,nrow(dta))
    for(i in (h+2):nrow(dta)){
        tmpMu = tmp1%*%matrix(dta[i-h,],nrow=K,ncol=1)+
          tmp2%*%matrix(dta[i-h-1,],nrow=K,ncol=1)
        log.score[i] = -dnormal(x=dta[i,],mu=tmpMu,
                                Sigma=tmpSigma,return.log=TRUE)
    }
    return(log.score)
}

#------------------------------------------------------------------------------------------------------------------------
# impulseResponseVAR(A,Sigma=NULL)
# compute the impulse response coefficients of a causal VAR model
#------------------------------------------------------------------------------------------------------------------------
impulseResponseVAR = function(A, Sigma=NULL, maxLag=10){
    A         = as.matrix(A)
    K         = nrow(A)
    p         = ncol(A)/K
    if(is.null(Sigma)){ 
        Sigma=diag(1,K)
    }
    if(min(eigen(Sigma)$values)<0){
        stop(" Error: the noise covariance matrix Sigma is not positive definite")
    }
    Sigma.chol             = t(chol(Sigma))
    companion              = matrix(0, nrow=K*p, ncol=K*p)
    companion[1:K,1:(K*p)] = A
    if (p>1){
        companion[(K+1):(K*p),1:(K*p-K)] = diag(1,K*p-K)
    }
    if(max(Mod(eigen(companion)$values))>=1){
        warning(" Warning: the VAR model is not causal")
    }
    result      = array(NA, dim=c(K,K,maxLag+1), 
                        dimnames=list(NULL,NULL, 0:maxLag))
    matrixPower = diag(1,K*p) 
    result[,,1] = matrixPower[1:K,1:K] %*% Sigma.chol
    for(i in 1:maxLag){
        matrixPower   = matrixPower%*%companion
        result[,,i+1] = matrixPower[1:K,1:K]%*%Sigma.chol
    }
    return(result)
}

#########################################################################
#####################Estimates bias, variance and 
#####################This function used for Replicating TABLE1

####@True.A - True coefficient matrix A
####@exp_estA - Expected Estimated coefficient matrix A
####@var_estA - Estimated Variance
####@max.p    - number of lag

bias_var_mse<-function(True.A, exp_estA, var_estA, max.p)
{
  ### This function estimates bias, variance,
  ### and MSE of coefficient matrix A #####
  A_augmen<-array(0, dim = c(nrow(True.A), nrow(True.A) * (max.p)))
  A_augmen[,1 : ncol(True.A)] = True.A
  bias<-array(0, dim = c(1, dim(exp_estA)[2]))
  var<-array(0,dim = c(1, dim(exp_estA)[2]))
  mse<-array(0,dim = c(1, dim(exp_estA)[2]))
  bias.temp <- 0
  var.temp<-0
  for (p in 1:max.p)
  {
    bias.tempi<-norm(exp_estA[,max(1,(p-1)*nrow(A_augmen)):
                                (p*nrow(A_augmen))]-
                       A_augmen[,max(1,(p-1)*nrow(A_augmen)):
                                  (p*nrow(A_augmen))], type = c( "F"))^2
    bias.temp = bias.temp + bias.tempi
    var.tempi<-sum(var_estA[,max(1,(p-1)*nrow(A_augmen)):(p*nrow(A_augmen))])
    var.temp = var.temp + var.tempi
    mse.temp=bias.temp+var.temp
  }
  bias=bias.temp
  var=var.temp
  mse=mse.temp
  return(list(bias_sq = bias, var = var, mse = mse))
}


######################################################################

## @estAdj - estimated Adjacency matrix
## @trueAdj - true Adjacency matrix

compareA <- function (estAdj, trueAdj) 
{
  #### This function compares esitmated
  ### adjacency matrix of A with the true Adj matrix 
  ### of coefficient matrix A
  ml <- estAdj
  mt <- trueAdj
  p <- dim(ml)[2]
  mt[mt != 0] <- rep(1, sum(mt != 0))
  ml[ml != 0] <- rep(1, sum(ml != 0))
  diffm <- ml - mt
  nmbTrueGaps <- (sum(mt == 0) - p)/2
  fpr <- if (nmbTrueGaps == 0) 
    1
  else (sum(diffm > 0)/2)/nmbTrueGaps
  diffm2 <- mt - ml
  nmbTrueEdges <- (sum(mt == 1)/2)
  tpr <- if (nmbTrueEdges == 0) 
    0
  else 1 - (sum(diffm2 > 0)/2)/nmbTrueEdges
  trueEstEdges <- (nmbTrueEdges - sum(diffm2 > 0)/2)
  tdr <- if (sum(ml == 1) == 0) {
    if (trueEstEdges == 0) 
      1
    else 0
  }
  else trueEstEdges/(sum(ml == 1)/2)
  c(tpr = tpr, fpr = fpr, tdr = tdr)
}

### This function plots the coefficient matrix
plot_mat <- function (Mat, main = NULL) 
{
  tmppar <- par(pty = "s")
  image(sign(t(apply(Mat, 2, rev))), axes = FALSE, col = c("gray50", 
                                                           "white", "black"), main = main)
  par(tmppar)
}


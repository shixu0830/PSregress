library(survival)
library(splines)
library(ATE)
library(tmle)
library(mvtnorm)
library(mgcv)
library(foreach)
library(Matrix)
library(glmnet)
expit = function(x){exp(x)/(1+exp(x))}
logit = function(p){log(p/(1-p))}

######################################################################
#####         key functions for methods under comparison         #####
##### Stratification; IPTW; Xadj; Unadj; PSadj; TMLE; DR.Xadj.PS; Matching
######################################################################
### step-func of PS: regress on indicators of strata
Stratification = function(dat,PS.omit,n.strata=5){
  my.formula = get_formula(dat,formula.head="z~",omit=PS.omit)
  ps = glm(as.formula(my.formula), data=dat,family=binomial(link = "logit"))$fitted #prob scale
  #strata = cut(ps,include.lowest=T,breaks=quantile(ps,seq(0,1,1/n.strata),labels=letters[1:n.strata]))
  strata = try( ## get PS strata indicators
    cut(ps,include.lowest=T,breaks=quantile(ps,seq(0,1,1/n.strata),names=F))
    ,silent=TRUE)
  if(inherits(strata, "try-error")){
    rslt = c(ATE=NA,se=NA,p1=NA,p0=NA,mybias=NA,b=NA)
  }else{
    dat.strata = cbind(dat,strata=strata)
    n = nrow(dat)
    dat1=data.frame(z=rep(1,n),strata=strata)
    dat0=data.frame(z=rep(0,n),strata=strata)
    m = try(glm(y~factor(z)+factor(strata),family="binomial", data = dat,
                control=glm.control(maxit = 25)), silent=TRUE)
    rslt = read.results(m,dat1,dat0,dat,ps)
  }
  
  return(rslt)
}
## IPTW: logistic regression with weighted observations
IPTW = function(dat,trim,PS.omit){
  my.formula = get_formula(dat,formula.head="z~",omit=PS.omit)
  ps = glm(as.formula(my.formula),data=dat,family=binomial(link = "logit"))$fitted #prob scale
  if(trim==T){ ## trim PS to exclude extreme val
    range = quantile(ps,probs=c(0.0025,0.9975))# so at most n*0.005=5 pts will have their ps changed
    Lbound = range[1]; Ubound = range[2]
    ps[which(ps<Lbound)]=Lbound
    ps[which(ps>Ubound)]=Ubound
  }
  my.weight = dat$z/ps + (1-dat$z)/(1-ps)
  m = try(glm(y~factor(z),family="binomial", data=dat,
              weights = my.weight,
              control=glm.control(maxit = 25)), silent=TRUE)
  if(inherits(m, "try-error")){
    biptw_z0=se_biptw_z0=p1=p0=NA
  }else{
    m_z0 = summary(m)$coefficients
    biptw_z0 = m_z0[2,1]; se_biptw_z0 = m_z0[2,2]
  }
  
  p = dat$y*my.weight
  p1 = mean(p[dat$z==1]); p0 = mean(p[dat$z==0])
  return(c(ATE=biptw_z0,se=se_biptw_z0,
           p1=p1,p0=p0
  ))
}
## g-computation: adjust for X then standardize
Xadj = function(dat,OutcomeR.omit){
  n = nrow(dat)
  my.formula = get_formula(dat,formula.head="z~",omit=F)
  ps = glm(as.formula(my.formula), data=dat,family=binomial(link = "logit"))$fitted #prob scale
  
  my.formula = get_formula(dat,formula.head="y~factor(z)+",omit=OutcomeR.omit)
  m = try(glm(as.formula(my.formula),
              family="binomial", data = dat,
              control=glm.control(maxit = 25)), silent=TRUE)
  dat1=data.frame(z=rep(1,n),dat[,grep("x",names(dat))])
  dat0=data.frame(z=rep(0,n),dat[,grep("x",names(dat))])
  rslt = read.results(m,dat1,dat0,dat,ps)
  return(rslt)
}
## naive: no confounding adjustment
Unadj = function(dat,OutcomeR.omit){
  n = nrow(dat)
  my.formula = get_formula(dat,formula.head="z~",omit=PS.omit)
  ps = glm(as.formula(my.formula), data=dat,family=binomial(link = "logit"))$fitted #prob scale
  
  my.formula = "y~factor(z)"
  m = try(glm(as.formula(my.formula),
              family="binomial", data = dat,
              control=glm.control(maxit = 25)), silent=TRUE)
  dat1=data.frame(z=rep(1,n),dat[,grep("x",names(dat))])
  dat0=data.frame(z=rep(0,n),dat[,grep("x",names(dat))])
  rslt = read.results(m,dat1,dat0,dat,ps)
  return(rslt)
}
### PS-regress: data-adaptive smoothing parameters; w/ and w/o treatment heterogeneity
PSadj = function(dat,PS.omit=F,trim){
  n = nrow(dat)
  my.formula = get_formula(dat,formula.head="z~",omit=PS.omit)
  ps = glm(as.formula(my.formula), data=dat,family=binomial(link = "logit"))$fitted #prob scale
  
  bsps2 = bs(ps,df=3)
  ind1 = target.smooth(dat,ps,interaction.spec=F)
  bsps3 = bs(ps,df=ind1$mydf,degree=ind1$mydegree)
  ind2 = target.smooth(dat,ps,interaction.spec=T)
  bsps4 = bs(ps,df=ind2$mydf,degree=ind2$mydegree)
  
  z = dat$z; y=dat$y
  #### consider four variates of PS-regress method
  ### (0) linear adjustment: biased
  dat.0 = data.frame(y,z,ps)
  ### (1) cubic spline: not nonparametric
  dat.1 = data.frame(y,model.matrix(~z+bsps2)[,-1]); formula1 = as.formula(~z+bsps2)
  ### (2) Bspline with data-adaptive smoothing param; homogeneous treatment effect
  dat.2 = data.frame(y,model.matrix(~z+bsps3)[,-1]); formula2 = as.formula(~z+bsps3)
  ### (3) Bspline with data-adaptive smoothing param; heterogeneous treatment effect
  dat.3 = data.frame(y,model.matrix(~z*bsps4)[,-1]); formula3 = as.formula(~z*bsps4)
  
  m.0 = try(glm(as.formula(paste0(names(dat.0)[1],"~",paste(names(dat.0)[2:ncol(dat.0)],collapse="+"))),
                family="binomial", data = dat.0,
                control=glm.control(maxit = 25)), silent=TRUE)
  dat1.0=data.frame(y,z=rep(1,n),ps)
  dat0.0=data.frame(y,z=rep(0,n),ps)
  rslt.0 = read.results(m.0,dat1.0,dat0.0,dat,ps)
  names(rslt.0)=paste0(names(rslt.0),".0")
  
  m.1 = try(glm(as.formula(paste0(names(dat.1)[1],"~",paste(names(dat.1)[2:ncol(dat.1)],collapse="+"))),
                family="binomial", data = dat.1,
                control=glm.control(maxit = 25)), silent=TRUE)
  z=rep(1,n);dat1.1=data.frame(y,model.matrix(formula1)[,-1])
  z=rep(0,n);dat0.1=data.frame(y,model.matrix(formula1)[,-1])
  rslt.1 = read.results(m.1,dat1.1,dat0.1,dat,ps)
  names(rslt.1)=paste0(names(rslt.1),".1")
  
  m.2 = try(glm(as.formula(paste0(names(dat.2)[1],"~",paste(names(dat.2)[2:ncol(dat.2)],collapse="+"))),
                family="binomial", data = dat.2,
                control=glm.control(maxit = 25)), silent=TRUE)
  z=rep(1,n);dat1.2=data.frame(y,model.matrix(formula2)[,-1])
  z=rep(0,n);dat0.2=data.frame(y,model.matrix(formula2)[,-1])
  rslt.2 = read.results(m.2,dat1.2,dat0.2,dat,ps)
  names(rslt.2)=paste0(names(rslt.2),".2")
  
  m.3 = try(glm(as.formula(paste0(names(dat.3)[1],"~",paste(names(dat.3)[2:ncol(dat.3)],collapse="+"))),
                family="binomial", data = dat.3,
                control=glm.control(maxit = 25)), silent=TRUE)
  z=rep(1,n);dat1.3=data.frame(y,model.matrix(formula3)[,-1])
  z=rep(0,n);dat0.3=data.frame(y,model.matrix(formula3)[,-1])
  rslt.3 = read.results(m.3,dat1.3,dat0.3,dat,ps)#m=m.3;dat1=dat1.3;dat0=dat0.3
  names(rslt.3)=paste0(names(rslt.3),".3")
  
  
  return(c(rslt.0,rslt.1,rslt.2,rslt.3,df1=unlist(ind1),df2=unlist(ind2)))
  
}
## TMLE using user-specified parametric model for Q (E[T|X,Z])
TMLE = function(w.dat,f.dat,PS.omit,OutcomeR.omit,add.superlearner){ 
  if(length(grep("x",names(w.dat)))!=length(grep("x",names(f.dat)))){
    myW = cbind(w.dat[,grep("x",names(w.dat))],f.dat[,grep("x",names(f.dat))])
  }else{
    myW = f.dat[,grep("x",names(f.dat))]
  }
  
  if(add.superlearner==F){
    g.formula = get_formula(w.dat,formula.head="A~",omit=PS.omit)
    Q.formula = get_formula(f.dat,formula.head="Y~factor(A)+",omit=OutcomeR.omit)
    m = tmle(Y=f.dat$y,A=f.dat$z,W=myW,
             family = "binomial",
             Qform = as.formula(Q.formula),
             gform = as.formula(g.formula)
    )
  }else{
    library(SuperLearner)
    ### selected a subset of the library to speed up
    mylib = list(
      c("SL.glm","screen.glmnet"),
      c("SL.glmnet","screen.glmnet")
    )
    m=tmle(Y=f.dat$y,A=f.dat$z,W=myW,
           family="binomial", 
           Q.SL.library = mylib,
           cvQinit=T,V=2,
           g.SL.library = mylib)
  }
  
  btmle_z0 = m$estimates$OR$log.psi
  se_btmle_z0 = sqrt(m$estimates$OR$var.log.psi)
  
  return(c(ATE=btmle_z0,se=se_btmle_z0,
           p1=mean(m$Qstar[,2]),p0=mean(m$Qstar[,1]),
           mybias=mean(m$estimates$IC$IC.logOR)
  ))
}
### doubly robust
DR.Xadj.PS = function(w.dat,f.dat,PS.omit,OutcomeR.omit,trim){
  n = nrow(w.dat)
  ## get the weight (PS): can omit confounder
  my.formula = get_formula(w.dat,formula.head="z~",omit=PS.omit)
  my.ps = glm(as.formula(my.formula), data=w.dat,family=binomial(link = "logit"))$fitted #prob scale
  if(trim==T){ ## trim PS to exclude extreme val
    range = quantile(my.ps,probs=c(0.0025,0.9975))# so at most n*0.0005=5 pts will have their ps changed
    Lbound = range[1]; Ubound = range[2]
    my.ps[which(my.ps<Lbound)]=Lbound
    my.ps[which(my.ps>Ubound)]=Ubound
  }
  my.weight = w.dat$z/my.ps + (1-w.dat$z)/(1-my.ps)
  my.formula = get_formula(f.dat,formula.head="y~factor(z)+",omit=OutcomeR.omit)
  m = try(glm(as.formula(my.formula),
              family="binomial", data = f.dat,
              control=glm.control(maxit = 25)), silent=TRUE)
  dat1=data.frame(z=rep(1,n),f.dat[,grep("x",names(f.dat))])
  dat0=data.frame(z=rep(0,n),f.dat[,grep("x",names(f.dat))])
  
  if(inherits(m, "try-error")){
    ATE.drx_ps=ATT.drx_ps=NA 
  }else{
    predict1 = predict(m,newdata=dat1,type="response")
    predict0 = predict(m,newdata=dat0,type="response")
    p1.drx_ps = mean((f.dat$z)/(my.ps)*(f.dat$y-predict1)+predict1)
    p0.drx_ps = mean((1-f.dat$z)/(1-my.ps)*(f.dat$y-predict0)+predict0)
    # p1.drx_ps = mean(f.dat$y*f.dat$z/my.ps)-mean(predict1*(f.dat$z-my.ps)/my.ps) ##alternatively
    # p0.drx_ps = mean(f.dat$y*(1-f.dat$z)/(1-my.ps)-predict0*(f.dat$z-my.ps)/(1-my.ps)) ##alternatively
    ATE.drx_ps = log((p1.drx_ps/(1-p1.drx_ps))/(p0.drx_ps/(1-p0.drx_ps)))
  }
  return(c(ATE=ATE.drx_ps,p1=p1.drx_ps,p0=p0.drx_ps))
}
### PS matching: not considered in the paper because it
### targets ATT, thus not comparable to other methods unless null effect, i.e. ATT=ATE=0
Matching = function(dat,PS.omit,my.replacement=T,my.M=1){
  my.formula = get_formula(dat,formula.head="z~",omit=PS.omit)
  m0 = glm(as.formula(my.formula), data=dat,family=binomial(link = "logit")) #prob scale
  ps = logit(m0$fitted)
  rr = Match(Y=NULL, Tr=dat$z, X=predict(m0), M=my.M, replace=my.replacement,caliper=0.2,ties=F)
  if(my.replacement == F){ind.control = unique(rr$index.control)}else{ind.control = rr$index.control}
  ind = c(unique(rr$index.treated),ind.control)
  id = c(1:length(unique(rr$index.treated)),rep(1:length(unique(rr$index.treated)),each=my.M))
  dat.match = cbind(dat[ind,],id=id,ps=ps[ind])
  m = try(glm(y~factor(z),family="binomial", data = dat.match,
              control=glm.control(maxit = 25)), silent=TRUE)
  if(inherits(m, "try-error")){
    b_z1 = NA; se_b_z1 = NA
  }else{
    m_z1 = summary(m)$coefficients
    b_z1 = m_z1[2,1]; se_b_z1 = m_z1[2,2]
  }
  return(c(ATE=b_z1,se=se_b_z1,
           p1=mean(dat.match$y[dat.match$z==1]),
           p0=mean(dat.match$y[dat.match$z==0])
  ))
}

#############################################
#####              others               #####  
#############################################
## target.smooth() and get.CV.mesure(): To get the smoothing parameter via cross-validation
get.CV.mesure = function(mydf,mydegree,dat.train,ps.train,dat.test,ps.test,interaction.spec){
  bsps4 = bs(ps.train,df=mydf,degree=mydegree); z=dat.train$z; y=dat.train$y
  if(interaction.spec==F){
    dat.3 = data.frame(y,model.matrix(~z+bsps4)[,-1])
  }else{
    dat.3 = data.frame(y,model.matrix(~z*bsps4)[,-1])
  }
  m.3 = try(glm(as.formula(paste0(names(dat.3)[1],"~",paste(names(dat.3)[2:ncol(dat.3)],collapse="+"))),
                family="binomial", data = dat.3,
                control=glm.control(maxit = 25)), silent=TRUE)
  bsps4 = bs(ps.test,df=mydf,degree=mydegree); z = dat.test$z; y=dat.test$y
  if(interaction.spec==F){
    z=rep(1,length(z));dat1.3=data.frame(y,model.matrix(~z+bsps4)[,-1])
    z=rep(0,length(z));dat0.3=data.frame(y,model.matrix(~z+bsps4)[,-1])
  }else{
    z=rep(1,length(z));dat1.3=data.frame(y,model.matrix(~z*bsps4)[,-1])
    z=rep(0,length(z));dat0.3=data.frame(y,model.matrix(~z*bsps4)[,-1])
  }
  
  rslt.3 = read.results(m.3,dat1.3,dat0.3,dat=dat.test,ps=ps.test)
  error = rslt.3["mybias"]^2+(rslt.3["se"]^2)
  # error = rslt.3["mybias"]^2+((rslt.3["se"])*sqrt(nrow(dat1.3)))^2
  return(error)
}
target.smooth = function(dat,ps,interaction.spec){
  n = nrow(dat)
  M = ceiling(n^(1/5))
  MSE = matrix(NA,nrow=M,ncol=M)
  ind.sample=sample(1:n,size=ceiling(n/2))
  dat.1 = dat[ind.sample,]; dat.2 = dat[-ind.sample,]
  ps.1 = ps[ind.sample]; ps.2 = ps[-ind.sample]
  ## Two-fold CV
  for(mydf in 3:M){
    for(mydegree in 3:mydf){
      ### fold one
      error1 = get.CV.mesure(
        mydf=mydf,mydegree=mydegree,
        dat.train=dat.1,ps.train=ps.1,
        dat.test=dat.2,ps.test=ps.2,
        interaction.spec
      )
      #### fold two
      error2 = get.CV.mesure(
        mydf=mydf,mydegree=mydegree,
        dat.train=dat.2,ps.train=ps.2,
        dat.test=dat.1,ps.test=ps.1,
        interaction.spec)
      MSE[mydegree,mydf] = error1+error2
    }
  }
  ind = which(MSE == min(MSE,na.rm=T), arr.ind = TRUE)
  return(list(mydegree=ind[1], mydf=ind[2]))
}
## To get ATE and var from fitted model "m"
read.results = function(m,dat1,dat0,dat,ps){
  if(inherits(m, "try-error")){
    b_zpbs0=se_b_zpbs0=ATE.ps=ATT.ps=p1.ps=p0.ps=NA
  }else{
    b_zpbs0 = summary(m)$coefficients[2,1]
    ## find ATE.ps and ATT.ps [standardization]
    predict1 = predict(m,newdata=dat1,type="response")
    predict0 = predict(m,newdata=dat0,type="response")
    
    ## compute p1, p0
    p1 = mean(predict1)
    p0 = mean(predict0)
    ATE = logit(p1)-logit(p0) 
  }
  
  ## compute variance
  # OR
  IF1 = predict1-p1 + dat$z/ps*(dat$y-predict1)
  IF0 = predict0-p0 + (1-dat$z)/(1-ps)*(dat$y-predict0)
  #IF = (1-p0) / ( p0*(1-p1)^2 )*IF1  -  (p1) / ( (1-p1)*(p0)^(2) )*IF0 ## IF for OR
  IF =  IF1/(p1*(1-p1))  -  IF0/(p0*(1-p0))  ## IF for logOR
  mysd = sqrt(  sum(IF^2) / (nrow(dat)^2)  )
  mybias = mean(IF)
  
  #   # RD var
  #   IF1 = predict1-p1 + dat$z/ps*(dat$y-predict1)
  #   IF0 = predict0-p0 + (1-dat$z)/(1-ps)*(dat$y-predict0)
  #   IF =  IF1 - IF0
  #   mysd.trimps = sqrt(   sum(IF^2) / (nrow(dat)^2)   )
  #   mybias.trimps = mean(IF)
  
  rslt = c(ATE=ATE,se=mysd,p1=p1,p0=p0,mybias=mybias,b=b_zpbs0)
  return(rslt)
}
get_formula = function(dat,formula.head=NULL,omit){
  if(omit==T){
    ind.omit = which(names(dat)=="x9")#grep("x9",names(dat))
    x.names = names(dat)[setdiff(grep("x",names(dat)),ind.omit)]
  }else{
    x.names = names(dat)[grep("x",names(dat))]
  }
  my.formula = paste0(formula.head,paste(x.names,collapse="+"))
  return(my.formula)
}


#############################################
#####  code for realistic simulation    #####  
#############################################
ordgendata <- function(n, sigma, quants.norm){
  retval = mvtnorm::rmvnorm(n = n, sigma = sigma)
  for (i in 1:ncol(sigma)) {
    retval[, i] = cut(x = retval[, i], breaks = c(-1/0, quants.norm[[i]]), 
                      right = FALSE)
  }
  retval - 1  
}
generate.cov.ord <- function(n, P.ord, Quant.norm, Corr.norm, alpha_xz)
{
  Z <- ordgendata(n, sigma=Corr.norm, quants.norm=Quant.norm)
  n.B <- length(which(unlist(lapply(P.ord,FUN=function(x){length(x)})==2)))
  n.A <- ncol(Z) - n.B
  
  if (n.B > 0 ) {
    B <- Z[, 1:n.B, drop=FALSE]
    colnames(B) <- paste0('bin', 1:n.B)
  }
  
  if(n.A>0)
  {
    A <- as.matrix(Z[,(n.B+1):ncol(Z)])  # rita added: as.matrix
    A.indicator <- NULL
    catlabs <- NULL
    for(i in 1:n.A)
    {
      dummy <- NULL
      levels <- sort(unique(A[,i]))
      for(level in levels[-1]) {
        catlabs <- c(catlabs, paste0('cat',i,'_','level',level))
        dummy <- cbind(dummy, A[,i]==level)
      }
      A.indicator <- cbind(A.indicator, dummy)
    }
    colnames(A.indicator) <- catlabs
  }
  
  ### Confounders Z
  Z.model.data <- as.matrix(cbind(B, A.indicator))
  
  return(list(B=B, A.indicator=A.indicator))
}
gen_XZ = function(n,p,X2.p,alpha0,alpha_xz,
                  P.ord, Quant.norm, Corr.norm,
                  highd=F,p.highd=100){
  #### generate X
  AB = generate.cov.ord(n, P.ord, Quant.norm, Corr.norm, alpha_xz)
  X <- cbind(AB$B, AB$A.indicator)
  if(highd==T){
    X = cbind(X,rmvnorm(n=nrow(X),mean=rep(0,p.highd)))
    alpha_xz = c(alpha_xz,rep(0,p.highd))
  }else{
    p.highd=0
  }
  ### alphaX: linear predictor in PS model
  alphaX = X%*%matrix(alpha_xz,nrow=length(alpha_xz))
  
  ## generate Z:
  ## E[Z|X]=expit(alpha0+alphaX)
  Z = rbinom(n=n,1,p = expit(alpha0+alphaX))
  XZ = list(Z=Z,X=X,alpha0=alpha0,alpha_xz=alpha_xz,highd=highd,p.highd=p.highd)
  return(XZ)
}

gen_Y = function(n,p,X2.p,gamma,gamma0,gamma_xy,X,Z,interaction
                 ,alpha0,alpha_xz,highd=F,p.highd=100){
  #### E[Y|Z,X]=expit(gamma0+betaX*Z+gammaX)
  if(highd==T){
    gamma_xy = c(gamma_xy,rep(coef.highd/p.highd,p.highd))
  }
  ## gammaX
  gammaX = X%*%matrix(gamma_xy[-1],nrow=length(gamma_xy[-1]))
  
  ### betaX is treatment effect
  if(interaction==T){
    betaX = X%*%matrix(rev(gamma_xy[-1]),nrow=length(gamma_xy[-1]))
  }else{
    betaX = gamma_xy[1]
  }
  
  myp0 = expit(gamma0+gammaX)
  myp1 = expit(gamma0+betaX+gammaX)
  Y0 = rbinom(n=n,size=1,p=myp0)
  Y1 = rbinom(n=n,size=1,p=myp1)
  p1 = mean(Y1)
  p0 = mean(Y0)
  ATE = logit(p1)-logit(p0)
  ATT = logit(mean(Y1[Z==1]))-logit(mean(Y0[Z==1]))
  ATC = logit(mean(Y1[Z==0]))-logit(mean(Y0[Z==0]))
  Y = (1-Z)*Y0+Z*Y1
  
  Yp=list(Y=Y,p0=p0,p1=p1,ATE=ATE,ATT=ATT,ATC=ATC)
  return(Yp)
}

gen_dat = function(XZ,Yp,p.highd){
  ## data sets
  x = XZ$X; z = XZ$Z
  y = Yp$Y;
  n = length(z); p = ncol(x)
  colnames(x) = paste("x",1:ncol(x),sep="")
  dat = data.frame(y=y,z=z,x)
  return(list(dat=dat))
}


#####################################################
#####  code to find the intercept that leads    #####  
#####  to the desired event rate for Y and Z    #####  
#####################################################
intercept_Z = function(rate.z=exposure.rate,alpha0,alpha_xz,myres,
                       P.ord, Quant.norm, Corr.norm){
  ## note: no highd or p.highd becuase assume coef of noise in XZ is zero
  ## find intercept for Z, fix rate.z = E[Z]
  N = 1000000
  count=0
  
  #### generate X
  AB = generate.cov.ord(N, P.ord, Quant.norm, Corr.norm, alpha_xz)
  X <- cbind(AB$B, AB$A.indicator)
  ### alphaX: linear predictor in PS model
  alphaX = X%*%matrix(alpha_xz,nrow=length(alpha_xz))
  
  x.l=-50#logit((l+u)/2)
  Z = rbinom(n=N,1,p = expit(x.l+alphaX))
  f.l=rate.z-mean(Z)
  
  x.r=50#logit((l+u)/2)
  Z = rbinom(n=N,1,p = expit(x.r+alphaX))
  f.r=rate.z-mean(Z)
  
  while(abs(x.l-x.r)>myres & count<1000){
    x.m=(x.l+x.r)/2
    Z = rbinom(n=N,1,p = expit(x.m+alphaX))
    f.m=rate.z-mean(Z)
    
    if (f.m == 0) {
      return(x.m)
    }
    else if (f.l * f.m < 0) {
      x.r <- x.m
      f.r <- f.m
    }
    else {
      x.l <- x.m
      f.l <- f.m
    }
    
    count = count+1
    # print(c(x.l,x.r,f.l,f.r))
  }
  
  return(x.m)
}

intercept_Y = function(rate.y=control.rate,gamma0,gamma_xy,myres
                       ,alpha0,alpha_xz,P.ord, Quant.norm, Corr.norm,
                       highd,p.highd){
  ## intercept_Y() func doesnt need to know "interaction" cuz holding Z=0 grp only
  
  N = 1000000
  ## now find gamma0 to fix rate.y=E[Y|Z=0]
  count=0
  
  ## generate X and Z
  AB = generate.cov.ord(N, P.ord, Quant.norm, Corr.norm, alpha_xz)
  X <- cbind(AB$B, AB$A.indicator)
  if(highd==T){
    X = cbind(X,rmvnorm(n=nrow(X),mean=rep(0,p.highd)))
    alpha_xz = c(alpha_xz,rep(0,p.highd))
    gamma_xy = c(gamma_xy,rep(coef.highd/p.highd,p.highd))
  }
  ### alphaX
  alphaX = X%*%matrix(alpha_xz,nrow=length(alpha_xz))
  Z = rbinom(n=N,1,p = expit(alpha0+alphaX))
  
  ## generate Y
  gammaX = X%*%matrix(gamma_xy[-1],nrow=length(gamma_xy[-1]))
  
  x.l=-50#logit((l+u)/2)
  Y = rbinom(n=N,size=1,p=expit(x.l+0+gammaX))
  f.l=rate.y-mean(Y[Z==0])
  
  x.r=50#logit((l+u)/2)
  Y = rbinom(n=N,size=1,p=expit(x.r+0+gammaX))
  f.r=rate.y-mean(Y[Z==0])
  
  while(abs(x.l-x.r)>myres & count<1000){
    x.m=(x.l+x.r)/2
    Y = rbinom(n=N,size=1,p=expit(x.m+0+gammaX))
    f.m=rate.y-mean(Y[Z==0])
    
    if (f.m == 0) {
      return(x.m)
    }
    else if (f.l * f.m < 0) {
      x.r <- x.m
      f.r <- f.m
    }
    else {
      x.l <- x.m
      f.l <- f.m
    }
    
    count = count+1
    # print(c(x.l,x.r,f.l,f.r))
  }
  return(x.m)
}

find_intercept = function(
  alpha_xz,alpha0,Common.P,control.rate,Corr.norm,Corr.ord,
  exposure.rate,gamma_xy,gamma0,norm.spec,P,P.ord,
  Quant.norm,Quants.norm,highd,p.highd){
  ## find_intercept() func doesnt need to know "interaction" cuz holding Z=0 grp only
  set.seed(1)
  assoc.param = c(1,1.5)
  
  newalpha0 = newgamma0 = rep(NA,length(assoc.param))
  for(my.time in 1:length(assoc.param)){
    alpha_xz = assoc.param[my.time]*alpha_xz_fix
    gamma_xy = gamma_xy_fix
    
    newalpha0[my.time] = intercept_Z(rate.z=exposure.rate,alpha0,alpha_xz,myres=1e-5,
                                     P.ord, Quant.norm, Corr.norm)
    ## generate Y #use same seed for any gamma,alpha,since X,Z are diff, Y will be diff
    newgamma0[my.time] = intercept_Y(rate.y=control.rate,gamma0,gamma_xy,myres=1e-4,
                                     alpha0=newalpha0[my.time],alpha_xz,P.ord, 
                                     Quant.norm, Corr.norm,highd,p.highd)  
  }
  intercepts = rbind(newalpha0,newgamma0)
  return(intercepts)
}


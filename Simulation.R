#############################################
######        Simulation code          ######
######    Contact: shixu@umich.edu     ######
#############################################
##### Notation: Y=outcome; Z=treatment; X=confounders
##### E[Z|X]=expit(alpha0+alphaX), alphaX=alpha_xz*X
##### E[Y|Z,X]=expit(gamma0+betaX*Z+gammaX), gammaX=gamma_xy*X
rm(list=ls())
n.rep= 100
#########################################
#####  simulation setting
#########################################
my.filepath = "***"
load(paste0(my.filepath,"Summ_Stat_GH.RData")) 
source(paste0(my.filepath,"functions_bin_NEW.R"))
args = commandArgs(TRUE)
seed = as.numeric(args[[1]])
rareZ65 = as.numeric(args[[2]])==1   ## rareZ=1: 0.69 ## rareZ=0: 0.2
PS.wrong=as.numeric(args[[3]])==1    ## whether propensity score model is misspecified
OR.wrong=as.numeric(args[[4]])==1    ## whether outcome model is misspecified
null = as.numeric(args[[5]])==1      ## null ATE or non-null
interaction = as.numeric(args[[6]])==1  ## whether the data is generated under heterogeneous treatment effect as : f(X)*Z (instead of betaX*Z)
### for example: myseed=1;rareZ65=F;PS.wrong=F;OR.wrong=F;null=F;interaction=F;my.time=1;PS.omit=PS.wrong;OutcomeR.omit=OR.wrong
highd=F;p.highd=NULL;coef.highd=NULL ## low-d scenario in the paper
n = 10000
if(rareZ65==T){ ## two scenarios: (1) majority treated (2) majority control 
  exposure.rate #0.651088: the observed rate in the real data
}else{
  exposure.rate = 0.2
}
if(null==T){
  gamma_xy_fix[1]=0 ##such that betaX=0 and myp0=myp1, i.e. E[Y(1)]=E[Y(0)]
}else{
  gamma_xy_fix[1]=log(3)
}
control.rate = control.rate*(100000/n) ## scale by sample size to keep No. events fixed 
###########################################################
##### for all scenarios (orig and strg confounding) 
##### pre-calculat the intercept to control 
##### the event rate of Y and Z (control.rate, exposure.rate)
###########################################################
intercepts = find_intercept(alpha_xz_fix,alpha0,Common.P,control.rate,Corr.norm,Corr.ord,
                            exposure.rate,gamma_xy,gamma0,norm.spec,P,P.ord,
                            Quant.norm,Quants.norm,highd,p.highd)
###############################################################
#####  run simulation under orig and strg confounding
###############################################################
assoc.param = c(1,1.5) ## strength of confounding
rslt=NULL ##2*73 matrix 
### first row is original strength of confounding eff
### second row is stronger confounding eff
for(my.time in 1:length(assoc.param)){##my.time=1 and 2: orig and strg confounding
  alpha_xz = assoc.param[my.time]*alpha_xz_fix
  gamma_xy = gamma_xy_fix
  alpha0 = intercepts[1,my.time]
  gamma0 = intercepts[2,my.time]
  #########################################
  #####  generate data
  #########################################
  ## generate confounder X and treatment Z
  set.seed(myseed)
  XZ = gen_XZ(n,p,X2.p,alpha0,alpha_xz,
              P.ord, Quant.norm, Corr.norm,highd,p.highd)
  
  ## generate outcome Y 
  Yp = gen_Y(n,p,X2.p,gamma,gamma0,gamma_xy,X=XZ$X,Z=XZ$Z,
             interaction=interaction,XZ$alpha0,XZ$alpha_xz,
             highd=XZ$highd,p.highd=XZ$p.highd)
  # exp(Yp$ATE) ## the true OR
  n.ctrl = sum(Yp$Y[XZ$Z==0]); n.trt = sum(Yp$Y[XZ$Z==1]); 
  
  ## generate full data set
  alldat = gen_dat(XZ,Yp,p.highd=XZ$p.highd)
  dat = alldat$dat
  ## check if similar to control.rate & exposure.rate
  (generated.control.rate=mean(Yp$Y[XZ$Z==0]));control.rate
  (generated.exposure.rate=mean(XZ$Z));exposure.rate
  
  #########################################
  #####  run PS methods
  #########################################
  ## whether PS/OR model is misspecified by omitting covariates
  PS.omit=PS.wrong;OutcomeR.omit=OR.wrong 
  ### (1) naive: no confounding adjustment
  unadj = Unadj(dat,OutcomeR.omit)
  names(unadj)=paste0("Unadj.",names(unadj))
  ### (2) IPTW with PS trimmed
  ipw = IPTW(dat,trim=T,PS.omit)
  names(ipw)=paste0("IPTW.",names(ipw))
  ### (3) g-computation
  xadj = Xadj(dat,OutcomeR.omit)
  names(xadj)=paste0("Xadj.",names(xadj))
  ### (4) PS-regress: data-adaptive smoothing parameters
  ### w/ and w/o treatment heterogeneity; PS not trimmed
  psadj = PSadj(dat,PS.omit,trim=F)
  names(psadj)=paste0("psadj.",names(psadj))
  ### (5) step-func of PS: regress on indicators of strata
  str5= Stratification(dat,PS.omit,n.strata=5)
  names(str5)=paste0("str5.",names(str5))
  ### (6) doubly robust with PS trimmed
  dr.ps = c(DR.Xadj.PS(w.dat=dat,f.dat=dat,PS.omit,OutcomeR.omit,trim=T))
  names(dr.ps)=paste0("DR.",names(dr.ps))
  ### (7) TMLE w/o super learner
  tmle0 = c(TMLE(w.dat=dat,f.dat=dat,PS.omit,OutcomeR.omit,add.superlearner=F))
  names(tmle0)=paste0("TMLEnoSL.",names(tmle0))
  ### (8) TMLE w super learner
  tmle1 = c(TMLE(w.dat=dat,f.dat=dat,PS.omit,OutcomeR.omit,add.superlearner=T))
  names(tmle1)=paste0("TMLEwSL.",names(tmle1))
  
  
  rslt.i = c(param=assoc.param[my.time],
             unadj,ipw,xadj,psadj,str5,
             dr.ps,tmle0,tmle1,
             ATE=Yp$ATE,ATT=Yp$ATT,ATC=Yp$ATC,
             p1=mean(Yp$p1),p0=mean(Yp$p0),
             generated.control.rate=generated.control.rate,
             generated.exposure.rate=generated.exposure.rate,
             n.ctrl=n.ctrl,n.trt=n.trt)
  abs(round(rslt.i[grep("ATE",names(rslt.i))]-Yp$ATE,5))
  
  rslt = rbind(rslt,rslt.i)
}







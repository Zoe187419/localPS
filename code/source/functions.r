library(ggplot2)

###########################
### local ps regression ###
## 
## Y is outcome
## X is treatment
## ps is propensity score
## h is bandwidth

localPSreg <- function(Y,ps,X,h=.1){
  V=beta=numeric(length(Y))
  CI = matrix(0,length(Y),2)
  # loop over observations and do a weighted regression around that obs
  # weight matrix -- 
  K <- matrix(0 , length(Y) , length(Y))
  for (i in 1:length(Y)){
    # bandwidth -- user specified for now
    #h = diff(range(ps))/15
    # weights
    K[i,] = dnorm( (ps[i]-ps)/h )/dnorm(0) 
    # standardize weights within X (useful when collapsing below)
    K[i,] = K[i,] / (X*sum(X*K[i,]) + (1-X)*sum((1-X)*K[i,]))
    # fit the model
    dta = data.frame(Y=Y,X=X,K=K[i,],ps=ps)
    D <- svydesign(id = ~1, weights = ~K, data=dta)
    M = svyglm(Y~X,D)
    # save output
    beta[i] = M$coef[2]
    V[i] = vcov(M)[2,2]
    CI[i,] = beta[i] + c(-1.96,1.96)*sqrt(V[i])
  }
  # collapse weights for estimation of ATT
  w <- colSums(K[X==1,])
  dta = data.frame(Y=Y,X=X,K=w)
  D <- svydesign(id = ~1, weights = ~w, data=dta)
  M = svyglm(Y~X,D)
  
  out = data.frame(Y=Y,X=X,ps=ps,beta=beta,LB=CI[,1],UB=CI[,2],V=V)
  return(list(out=out,ATT=summary(M)$coef[2,]))
}



##############################
### local ps/ps regression ###
## 
## Y is outcome
## X is treatment
## ps is propensity score
## prog is the prognostic score
## h is bandwidth
## estimand is either ATE or ATT

localPSPSreg <- function(Y,ps,prog,X,h=.1 , estimand="ATT"){
  V=beta=numeric(length(Y))
  CI = matrix(0,length(Y),2)
  K <- matrix(0 , length(Y) , length(Y))
  # loop over observations and do a weighted regression around that obs
  for (i in 1:length(Y)){
    # bandwidth -- user specified for now
    #h = diff(range(ps))/15
    # weights
    K[i,] = dnorm( sqrt( (ps[i]-ps)^2 + (prog[i]-prog)^2 )/h )/dnorm(0)
    # standardize weights within X (useful when collapsing below)
    K[i,] = K[i,] / (X*sum(X*K[i,]) + (1-X)*sum((1-X)*K[i,]))
    
    # fit the model
    dta = data.frame(Y=Y,X=X,K=K[i,])
    D <- svydesign(id = ~1, weights = ~K, data=dta)
    M = svyglm(Y~X,D)
    # save output
    beta[i] = M$coef[2]
    V[i] = vcov(M)[2,2]
    CI[i,] = beta[i] + c(-1.96,1.96)*sqrt(V[i])
  }
  # collapse weights for estimation of ATT
  w <- colSums(K[X==1,])
  dta = data.frame(Y=Y,X=X,K=w)
  D <- svydesign(id = ~1, weights = ~w, data=dta)
  M = svyglm(Y~X,D)
  
  out = data.frame(Y=Y,X=X,ps=ps,prog=prog,beta=beta,var=V,LB=CI[,1],UB=CI[,2])
  return(list(out=out,ATT=summary(M)$coef[2,]))
}

###############################
###  heat map for localPSPSreg ###
##
## fit is "out" from localPSPSreg
## degree is degree of loess smooth
## span is span of loess smooth

heatmap.localPSPSreg <- function(fit,degree=1,span=0.1,...){
  ##############################
  ## pull the fit information ##
  # smooth using loess
  # and expand to a regular grid
  # use loess to put heatmap data on a regular grid
  temp=data.frame(ps=fit$ps,prog=fit$prog,beta=fit$beta)
  fit.loess = loess(beta ~ ps*prog, data = temp, degree = degree, span = span , normalize=F)
  # create regular grid and predict at grid points
  temp = expand.grid(list(ps = seq(min(fit$ps),max(fit$ps), length.out=100), prog = seq(min(fit$prog), max(fit$prog), length.out=100)))
  temp[,"TreatmentEffect"] = as.numeric(predict(fit.loess, newdata = temp))
  
  ########################################
  # set up range and colors for heatmap ##
  a = max(abs(temp$TreatmentEffect),na.rm=T)
  lims = c(-a,a)
  midpoint = 0
  mid = "white"
  low = "red"
  high = "blue"
  # this plot is to grab the y-axis labels AND the legend
  p2 = ggplot() + geom_tile(data=temp, aes_string(x="ps", y="prog", fill = "TreatmentEffect")) +
    xlab("Propensity score") + ylab("Prognostic score") +
    scale_fill_gradient2(limits = lims, midpoint=midpoint , low=low ,mid=mid, high = high) +
    scale_x_continuous(expand = c(0,0)) +
    scale_y_continuous(expand = c(0,0)) +
    theme(plot.margin=unit(c(1,0,0,5),"lines"),axis.ticks=element_blank()) +
    geom_point(aes(x=ps,y=prog),size=.5,data=fit[fit$X==1,])
  p2
}

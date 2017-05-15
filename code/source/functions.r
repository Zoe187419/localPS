library(ggplot2)

##############################
### local ps/ps regression ###
## 
## Y is outcome
## X is treatment
## ps is propensity score
## prog is the prognostic score
## h.ps is bandwidth for propensity score
## h.prog is bandwidth for prognostic score
## estimand is either ATE or ATT

localPSreg <- function(Y,X,ps=NULL,prog=NULL,h.ps=.1, h.prog=1 , estimand="ATT"){
  if(is.null(ps) & is.null(prog)){
    stop("At least one of ps or prog must be specified.")
  }
  if(is.null(ps)){
    h.ps = NULL 
  }
  if(is.null(prog)){
    h.prog = NULL 
  }
  
  # set up
  ps = cbind(ps=ps , prog=prog)
  h = c(h.ps , h.prog)
  
  V=beta=numeric(length(Y))
  CI = matrix(0,length(Y),2)
  K <- matrix(0 , length(Y) , length(Y))
  # loop over observations and do a weighted regression around that obs
  for (i in 1:length(Y)){
    # bandwidth -- user specified for now
    #h = diff(range(ps))/15
    # weights
    K[i,] = exp( -0.5 * colSums( ((t(ps)-ps[i,])/h)^2) )
    
    # poorly specified bandwidth may cause K=0 for all treated or control
    w.trt = sum(X*K[i,])
    w.con = sum((1-X)*K[i,])
    if (w.trt==0 | w.con==0){
      warning(paste0("Observation number ",i," is too far from opposite treatment group. Treatment effect not calculated at this observation. Consider increasing the bandwidth to avoid this warning."))
      
      # we cant calculate a treatment effect here with the given bandwidth, 
      # set K=0 everywhere and output is NA
      K[i,] = 0
      beta[i] = NA
      V[i] = NA
      CI[i,] = NA
    }else{
      # standardize weights within X (useful when collapsing below)
      K[i,] = K[i,] / (X*w.trt + (1-X)*w.con)
      
      # fit the model
      dta = data.frame(Y=Y,X=X,K=K[i,])
      D <- svydesign(id = ~1, weights = ~K, data=dta)
      M = suppressWarnings(svyglm(Y~X,D))
      # save output
      beta[i] = M$coef[2]
      V[i] = vcov(M)[2,2]
      CI[i,] = beta[i] + c(-1.96,1.96)*sqrt(V[i])
    }
  }
  # collapse weights for estimation of ATT
  w <- colSums(K[X==1,])
  dta = data.frame(Y=Y,X=X,K=w)
  D <- svydesign(id = ~1, weights = ~w, data=dta)
  
  #---------------------------------------
  # edited by cindy
  #M = suppressWarnings(svyglm(Y~X,D))
  fitted.model <- suppressWarnings(svyglm(Y~X,D))
	Yhat <- predict(fitted.model, type="response") 
	RMSE <- sqrt(mean((Yhat-Y)^2))
	ATT.estimate <- coef(fitted.model)["X"]
	ATT.se <- SE(fitted.model)["X"]
	ATT.L95 <-confint(fitted.model)["X","2.5 %"]
	ATT.U95 <-confint(fitted.model)["X","97.5 %"]
#---------------------------------------
  
  out = data.frame(Y=Y,X=X,ps,beta=beta,var=V,LB=CI[,1],UB=CI[,2])
  return(list(out=out, ATT.estimate=ATT.estimate, ATT.se=ATT.se, ATT.L95=ATT.L95, ATT.U95=ATT.U95, RMSE=RMSE)) #edited by cindy to return ATT
  #return(list(out=out,ATT=summary(M)$coef[2,]))
}

###############################
###  heat map for localPSreg ###
##
## fit is "out" from localPSPSreg
## degree is degree of loess smooth
## span is span of loess smooth

heatmap.localPSreg <- function(fit,degree=1,span=0.1,...){
  if (!all( c("ps","prog")%in%colnames(fit))){
    stop("localPSreg must be fit with two scores.")
  }
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

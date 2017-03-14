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
  for (i in 1:length(Y)){
    # bandwidth -- user specified for now
    #h = diff(range(ps))/15
    # weights
    K = dnorm( (ps[i]-ps)/h )/dnorm(0) 
    # fit the model
    dta = data.frame(Y=Y,X=X,K=K,ps=ps)
    D <- svydesign(id = ~1, weights = ~K, data=dta)
    M = svyglm(Y~X,D)
    # save output
    beta[i] = M$coef[2]
    V[i] = vcov(M)[2,2]
    CI[i,] = beta[i] + c(-1.96,1.96)*sqrt(V[i])
  }
  return(list(beta=beta,CI=CI,V=V))
}

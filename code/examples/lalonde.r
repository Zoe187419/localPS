rm(list=ls())
source("code/source/functions.r")
library(twang)
data(lalonde)

# fit the PS model
ps.lalonde <- ps(treat ~ age + educ + black + hispan + nodegree +
                   married + re74 + re75,
                 data = lalonde, stop.method=c("ks.max"),
                 estimand = "ATE",verbose=FALSE)

# save the ps
lalonde$ps = ps.lalonde$ps[,1]

# fit LoWePS -- propensity score only!
out <- localPSreg(Y=lalonde$re78,X=lalonde$treat,ps=lalonde$ps,h.ps=0.1)

# check that the ATT estimate works
temp = out$out
setDT(temp)

# just the mean of beta's at treated propensity scores
temp[X==1,mean(beta)]

# compare to calculation using weights
out$ATT

# generate plot!
library(data.table)
setDT(lalonde)
lalonde[, beta:= out$out[,"beta"]]
setkey(lalonde , ps)

pdf("figures/lalonde_example.pdf",8,6)
lalonde[ , plot(x=ps,y=beta , type='l')]
dev.off()


# find prognostic score
lalonde$prog = predict( lm(re78 ~ age + educ + black + hispan + nodegree +
                             married + re74 + re75 , data=lalonde[treat==0,]) , newdata=lalonde)

# scale to 0-1
lalonde[ , prog := (prog-min(prog))/(max(prog)-min(prog))]

# fit LoWePSPS -- propensity score and prognostic score
out <- localPSreg(Y=lalonde$re78,X=lalonde$treat,ps=lalonde$ps,prog=lalonde$prog,h.ps=0.1,h.prog=0.1)

# check that the ATT estimate works
temp = out$out
setDT(temp)

# just the mean of beta's at treated propensity scores
temp[X==1,mean(beta)]

# compare to calculation using weights
out$ATT

# generate heatmap!
heatmap.localPSreg(out$out)


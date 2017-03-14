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

# fit LoWePS-QR
out <- localPSreg(Y=lalonde$re78,X=lalonde$treat,ps=lalonde$ps,h=0.1)

# generate plots!
library(data.table)
setDT(lalonde)
lalonde[, beta:= out$beta]
setkey(lalonde , ps)

pdf("figures/lalonde_example.pdf",8,6)
lalonde[ , plot(x=ps,y=beta , type='l')]
dev.off()

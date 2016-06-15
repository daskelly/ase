#!/usr/bin/env Rscript
library(rjags)
library(runjags)
rjags::load.module("mix")       # for dbetabin in JAGS

source('../simulate_data.R')
dat <- simulate.data()

mat <- data.frame(Y=dat$Y, N=dat$N)
mat$p <- apply(mat, 1, function(row) binom.test(row[1], row[2])$p.value)
# set cluster1 for the SNP with least ASE and cluster2 the most, 
# according to the binomial exact test. This is as per instructions in
# http://doingbayesiandataanalysis.blogspot.com/2012/06/mixture-of\
# -normal-distributions.html, for implementing a mixture model in JAGS
isNull <- rep(NA, dat$n.snps)
isNull[which.min(mat$p)] <- FALSE
isNull[which.max(mat$p)] <- TRUE
# also, for any SNPs with p=0 or p=1, set cluster to 0
# probNull[mat$p == 0 | mat$p == 1] <- FALSE
datlist <- list(nsnps=dat$n.snps, Y=dat$Y, N=dat$N)

### Run the model and examine results ###
mod <- runjags::run.jags("biased_snps.bug", data=datlist, n.chains=4,
    monitor=c("isNull", "param1", "param2"))

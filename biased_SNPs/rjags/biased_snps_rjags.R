#!/usr/bin/env Rscript
library(rjags)
library(runjags)
rjags::load.module("mix")       # for dbetabin in JAGS

### Run the model and examine results ###
mod <- runjags::run.jags("biased_snps.bug", data=datlist, n.chains=4,
    monitor=c("isNull", "param1", "param2"))

# Model for identifying biased SNPs for subsequent removal

## Overview

The code in this directory implements a model for identifying
biased SNPs described in section 1.3.1 of the supplementary 
material in 
[doi:10.1101/gr.119784.110](https://dx.doi.org/10.1101/gr.119784.110). 
There are two separate implementations of the model,
both of which are run using R and should produce very similar results:

1. The original implementation used in the paper, 
available in the directory [orig](orig). The code in this directory
was not originally released with the paper but was used to 
filter out putatively biased SNPs.
2. An implementation of the model using [JAGS](http://mcmc-jags.sourceforge.net/). JAGS is a general framework
for simulation from Bayesian hierarchical models using MCMC. Here we
call JAGS from R using the runjags and rjags packages.

## Dependencies

* [orig](orig): the optparse R package
* [JAGS](JAGS): the rjags and runjags R packages

## Which implementation should you choose? 

The choice is up to you.
The first implementation has been more thoroughly tested but is 
relatively slow and the code is specific to this paper so it may
be difficult to troubleshoot.
The second implementation should be significantly faster and much
easier to troubleshoot, but should be considered a preview release 
in that it may still contain bugs.
Over the longer term, development will focus on the JAGS implementation
and the original code will be archived.

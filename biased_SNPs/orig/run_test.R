# Run a test of the model using simulated data
source('../simulate_data.R')

dat <- simulate.data()
n.snps <- length(dat$Y)
# columns 1 and 2 are currently ignored.
df <- data.frame(gene=rep("gene", n.snps), SNP=1:n.snps, dat$Y, dat$N-dat$Y)
filename <- "test_SNP_counts.txt"
write.table(df, file=filename, col.names=F, row.names=F, sep="\t", quote=F)

n.iter <- 100000
thin <- 100
source('2componentNullModel.bySNP.R')
# the above should take less than a minute

sapply(result$final, function(item) 
    quantile(item$sims, c(.025, .25, .75, .975)))


# plot a hist of posterior prob(ASE)

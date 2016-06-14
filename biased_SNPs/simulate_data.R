# Simulate some data that we can use to test whether the model
# is giving us reasonable results.
# Note that we are simulating data according to the exact 
# statistical model we implemented, so this should be an easy task.
simulate.data <- function(seed=1234) {
    set.seed(seed)
    n.snps <- 10000
    pi0 <- 0.95
    alpha <- rlnorm(1, 4.3, 1.8)
    delta <- runif(1)
    epsilon <- runif(1)
    avg.coverage <- 50
    N <- rpois(n.snps, lambda=avg.coverage)
    # generate beta-binomially distributed values:
    isNull <- runif(n.snps) < pi0
    Y <- isNull*rbinom(n.snps, N, rbeta(n.snps, alpha, alpha)) + 
        (1 - isNull)*rbinom(n.snps, N, rbeta(n.snps, delta, epsilon))
    list(Y=Y, N=N, isNull=isNull, alpha=alpha, delta=delta, epsilon=epsilon)
}

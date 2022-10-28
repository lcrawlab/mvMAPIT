library('mvMAPIT')

ind <- 1e2
nsnp <- 100
H2 <- 0.6
rho <- 0.5
maf <- 0.05 + 0.45*runif(nsnp)
X   <- (runif(ind*nsnp) < maf) + (runif(ind*nsnp) < maf)
X   <- matrix(as.double(X),ind,nsnp,byrow = TRUE)
s <- sample.int(10000, 1)
sim <- simulate_traits(X,
                           n_causal = 80,
                           H2 = H2,
                           rho = rho,
                           logLevel = 'DEBUG',
                           seed = s)
Xmean<-apply(X, 2, mean); Xsd<-apply(X, 2, sd); X<-t((t(X)-Xmean)/Xsd)

mapit <- mvmapit(t(X), t(sim$phenotype), cores = 32, logLevel = 'DEBUG')

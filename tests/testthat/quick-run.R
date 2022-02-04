library('mvMAPIT')

ind <- 1e2
nsnp <- 100
maf <- 0.05 + 0.45*runif(nsnp)
X   <- (runif(ind*nsnp) < maf) + (runif(ind*nsnp) < maf)
X   <- matrix(as.double(X),ind,nsnp,byrow = TRUE)
Y <- simulate_phenotypes(X, H2 = 1.0, rho = 0.0, logLevel = 'DEBUG')
Xmean<-apply(X, 2, mean); Xsd<-apply(X, 2, sd); X<-t((t(X)-Xmean)/Xsd)

mapit <- MvMAPIT(t(X), t(Y), cores = 32, logLevel = 'DEBUG')

tns <- c(9, 9, 9, 0, 0)
p.true <- c(0.1, 0.2, 0.33, 0.5, 0.7)
K <- length(p.true)
tys <- rbinom(K, tns, p.true);tys

gen.prior <- function(K, J=10000){
    pss.prior <- t(apply(matrix(runif(K*J), ncol=K), 1, sort))
    pss.prior
}

pss.prior <- gen.prior(K, J=1e5)

kpidx.fn <- function(pss.prior, tys, tns){
    K <- length(tys)
    Num <- dim(pss.prior)[1]
    dist.fn <- function(i){
         ps.gen <- pss.prior[i, ]
         tys.gen <- rbinom(K, tns, ps.gen)
         data.dist <- sum((tys - tys.gen)**2)
         data.dist <= 1
    }
    kp.idx <- sapply(1:Num, dist.fn)
    kp.idx
}



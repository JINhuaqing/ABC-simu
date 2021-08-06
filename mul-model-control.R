library(arrApply)

ascending.fn <- function(vec){
    diff.vec <- diff(vec)
    ndose <- length(vec)
    if (sum(diff.vec >=0)==ndose){
        rv <- vec
    }else{
        rv <- vec
        for (i in 2:ndose){
            if ( rv[i] < rv[i-1])
                rv[i] <- rv[i-1]
        }
    }
    return(rv)
}

prob.k <- function(pss, dlts, nums){
    J <- dim(pss)[1]
    numsMat <- matrix(rep(nums, J), nrow=J, byrow=T)
    dltsMat <- matrix(rep(dlts, J), nrow=J, byrow=T)
    pps <- pss ** dltsMat * (1-pss)**(numsMat-dltsMat)
    ws <- arrApply(pps, 2, "prod")
    w <- mean(ws)
    w
}

gen.u.rand <- function(k, K=5, phi=0.3, delta=0.1){
    cps <- c(runif(1, min=phi-delta, max=phi+delta))
    #cps <- c(phi)
    if (k > 1){
        cps <- c(cps, runif(k-1, min=0, max=phi-delta))
    }
    if (k < K){
        cps <- c(cps, runif(K-k, min=phi+delta, 1))
    }
    sort(cps)
}

gen.u.rand.over <- function(k, K=5, phi=0.3, delta=0.1){
    cps <- c()
    cps <- c(cps, runif(K-k+1, min=phi+delta, max=1))
    cps <- c(cps, runif(k-1, min=0, max=phi+delta))
    sort(cps)
}

gen.u.rand.under <- function(k, K=5, phi=0.3, delta=0.1){
    cps <- c()
    cps <- c(cps, runif(k, min=0, max=phi-delta))
    cps <- c(cps, runif(K-k, min=phi-delta, max=1))
    sort(cps)
}

# generate the scenarios for Pr(Yn|M_n, Ak)
gen.mu.rand.mul <- function(k, J, K=5, phi=0.3, delta=0.1, typ=2){
    if (typ == 1){
        pss <- lapply(1:J, function(i)gen.u.rand.under(k, K, phi, delta))
    }else if(typ==2){
        pss <- lapply(1:J, function(i)gen.u.rand(k, K, phi, delta))
    }else{
        pss <- lapply(1:J, function(i)gen.u.rand.over(k, K, phi, delta))
    }
    
    pssMat <- do.call(rbind, pss)
    pssMat
}

p.true <- c(0.2, 0.45, 0.5, 0.6, 0.8)
tns <- c(9, 9, 9, 9, 9)
phi <- 0.3
delta <- 0.05
J <- 1000
ndose <- length(tns)
pss <- lapply(1:ndose, function(k)gen.mu.rand.mul(k, J=J, K=ndose, phi=phi, delta=delta, typ=2))
pss.under <- lapply(1:ndose, function(k)gen.mu.rand.mul(k, J=J, K=ndose, phi=phi, delta=delta, typ=1))
pss.over <- lapply(1:ndose, function(k)gen.mu.rand.mul(k, J=J, K=ndose, phi=phi, delta=delta, typ=3))

res <- list()
tys <- rbinom(5, 9, p.true)
pks <- sapply(1:ndose, function(k)prob.k(pss[[k]], tys, tns))
pks.over <- sapply(1:ndose, function(k)prob.k(pss.over[[k]], tys, tns))
pks.under <- sapply(1:ndose, function(k)prob.k(pss.under[[k]], tys, tns))
pkss <- cbind(pks.under, pks, pks.over)
tb <- apply(pkss, 1,FUN=function(i)which.max(i)) 
tb






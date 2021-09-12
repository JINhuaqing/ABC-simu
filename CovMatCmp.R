# num.sps <- 10000
# # Method 1
# m1.x2 <- rbeta(num.sps, 2, 2)
# m1.x1 <- runif(num.sps, 0, m1.x2)
# m1.x3 <- runif(num.sps, m1.x2, 1)
# sps1 <- t(rbind(m1.x1, m1.x2, m1.x3))
# 
# 
# # Method 2
# m2.x1 <- runif(num.sps)
# m2.x2 <- runif(num.sps)
# m2.x3 <- runif(num.sps)
# m2.xs <- rbind(m2.x1, m2.x2, m2.x3)
# m2.xs.sorted <- apply(m2.xs, 2, sort)
# m2.x1 <- m2.xs.sorted[1, ]
# m2.x2 <- m2.xs.sorted[2, ]
# m2.x3 <- m2.xs.sorted[3, ]
# sps2 <- t(rbind(m2.x1, m2.x2, m2.x3))
# 
# # Method 4
# m4.x2 <- rbeta(num.sps, 1, 1)
# m4.x1 <- runif(num.sps, 0, m4.x2)
# m4.x3 <- runif(num.sps, m4.x2, 1)
# sps4 <- t(rbind(m4.x1, m4.x2, m4.x3))
# 
# 
# det(cov(sps1))
# det(cov(sps2))
# det(cov(sps4))

# geneerate random scenarios
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

# generate the scenarios for Pr(Yn|M_n, Ak)
gen.mu.rand <- function(k, J, K=5, phi=0.3, delta=0.1){
    pss <- lapply(1:J, function(i)gen.u.rand(k, K, phi, delta))
    pssMat <- do.call(rbind, pss)
    pssMat
}


gen.u.rand2 <- function(k, K=5, phi=0.3, delta=0.1){
    cps <- rep(0, K)
    cps[k] <- runif(1, phi-delta, phi+delta)
    if (k > 1){
        for (j in (k-1):1){
            cps[j] <- runif(1, min=0, max=min(cps[j+1], phi-delta))
        }
    }
    if (k < K){
        for (j in (k+1):K){
            cps[j] <- runif(1, min=max(cps[j-1], phi+delta), max=0.8)
        }
    }
    cps
}

gen.mu.rand2 <- function(k, J, K=5, phi=0.3, delta=0.1){
    pss <- lapply(1:J, function(i)gen.u.rand2(k, K, phi, delta))
    pssMat <- do.call(rbind, pss)
    pssMat
}

deltas <- seq(0.02, 0.2, 0.02)

k <- 1
K <- 8
phi <- 0.3
J <- 10000
trss <- list()
detss  <- list()

flag <- 1
for (delta in deltas){
    sps1 <- gen.mu.rand(k, J, K, phi, delta)
    sps2 <- gen.mu.rand2(k, J, K, phi, delta)
    trs <- c(sum(diag(cov(sps1))), sum(diag(cov(sps2))))
    dets <- c(det(cov(sps1)),  det(cov(sps2)))
    trss[[flag]] <- trs
    detss[[flag]] <- dets
    flag <- 1 + flag
}


par(mfrow=c(1, 2))
trssMat <- do.call(rbind, trss)
logTrs <- log(trssMat)
plot(deltas, logTrs[, 1], type="l", ylim=c(min(logTrs), max(logTrs)), ylab="Log of Trace for CovMat")
lines(deltas, logTrs[, 2], lty=2, col=2)
legend("topright", c("MCA", "NOC"), col=1:2, lty=1:2)

detssMat <- do.call(rbind, detss)
logDets  <- log(detssMat)
plot(deltas, logDets[, 1], type="l", ylab="Log of Det for CovMat", ylim=c(min(logDets), max(logDets)))
lines(deltas, logDets[, 2], lty=2, col=2)
legend("topright", c("MCA", "NOC"), col=1:2, lty=1:2)

# MCA2.simu.fn is the method we use in the paper

#setwd("C:/Users/Dell/Documents/ProjectCode/phaseI/Rcode")
#setwd("/Users/jinhuaqing/Documents/Projects_Code/phaseI/Rcode")
#setwd("C:/Users/JINHU/Documents/ProjectCode/MCA")
source("utilities.R")
library(magrittr)
library(BOIN)
library(arrApply)


gen.u.rand <- function(k, K=5, phi=0.3, delta=0.1){
    #cps <- c(phi)
    if (k > 0) {
        cps <- c(runif(1, phi-delta, phi+delta))
        if (k > 1){
            cps <- c(cps, runif(k-1, min=0, max=phi-delta))
        }
        if (k < K){
            cps <- c(cps, runif(K-k, min=phi+delta, 1))
        }
    }else{
            cps <- runif(K, min=phi+delta, 1)
    }
    sort(cps)
}

# generate the scenarios for Pr(Yn|M_n, Ak)
gen.mu.rand <- function(k, J, K=5, phi=0.3, delta=0.1){
    pss <- lapply(1:J, function(i)gen.u.rand(k, K, phi, delta))
    pssMat <- do.call(rbind, pss)
    pssMat
}

gen.prior <- function(K, phi, J=1e3, delta=0.05){
    pss <- lapply(0:K, function(k)gen.mu.rand(k, J=J, K=K, phi=phi, delta=delta))
    #pss.prior <- t(apply(matrix(runif(K*J), ncol=K), 1, sort))
    pss.prior <- do.call(rbind, pss)
    pss.prior
}


# kpidx.fn <- function(pss.prior, tys, tns){
#     K <- length(tys)
#     Num <- dim(pss.prior)[1]
#     dist.fn <- function(i){
#          ps.gen <- pss.prior[i, ]
#          tys.gen <- rbinom(K, tns, ps.gen)
#          data.dist <- sum((tys - tys.gen)**2)
#          data.dist <= log(sum(tns))
#     }
#     kp.idx <- sapply(1:Num, dist.fn)
#     kp.idx
# }

# better way to generate psudo-data
kpidx.fn <- function(pss.prior, tys, tns){
    K <- length(tys)
    Num <- dim(pss.prior)[1]
    pss.prior.vec <- as.vector(t(pss.prior))
    tns.vec <- rep(tns, Num)
    tys.gen <- rbinom(length(tns.vec), tns.vec, pss.prior.vec)
    tys.vec <- rep(tys, Num)
    tns.mat <- matrix(tns.vec, ncol=K, byrow=T)
    diff.mat <- matrix(tys.vec - tys.gen, ncol=K, byrow=T)
    # avoid dividing-by-zero problem
    tns.mat[tns.mat==0] <- 0.1
    rate.diff.mat <- diff.mat / tns.mat
    kp.idx <- rowSums(rate.diff.mat**2) <= 0.05
    kp.idx
}



# Simulation function for MCA
MCAABC.simu.fn <- function(phi, p.true, ncohort=12, init.level=1, 
                              cohortsize=1, add.args=list()){
    # phi: Target DIL rate
    # p.true: True DIL rates under the different dose levels
    # ncohort: The number of cohorts
    # cohortsize: The sample size in each cohort
    # alp.prior, bet.prior: prior parameters
    #set.seed(2)
    earlystop <- 0
    ndose <- length(p.true)
    cidx <- init.level
    
    tys <- rep(0, ndose) # number of responses for different doses.
    tns <- rep(0, ndose) # number of subject for different doses.
    tover.doses <- rep(0, ndose) # Whether each dose is overdosed or not, 1 yes
    pss.prior <- gen.prior(ndose, phi=phi, J=add.args$J, delta=add.args$delta)

    
    
    
    
    for (i in 1:ncohort){
        pc <- p.true[cidx] 
        
        # sample from current dose
        cres <- rbinom(cohortsize, 1, pc)
        
        # update results
        tys[cidx] <- tys[cidx] + sum(cres)
        tns[cidx] <- tns[cidx] + cohortsize
        
        
        
        cy <- tys[cidx]
        cn <- tns[cidx]
        
        add.args <- c(list(y=cy, n=cn, tys=tys, tns=tns, cidx=cidx), add.args)
        
        # if (overdose.fn(phi, add.args)){
        #     tover.doses[cidx:ndose] <- 1
        # }
        # 
        # if (tover.doses[1] == 1){
        #     earlystop <- 1
        #     break()
        # }
        
        

        kp.idx <- kpidx.fn(pss.prior, tys, tns)
        post.sps <- pss.prior[kp.idx, ]
        over.prob <- mean(post.sps[, 1]>phi)
        if (over.prob >= 0.90){
            earlystop <- 1
            break()
        }

        cor.ps <- c()
        for (j in 1:ndose){
            c.sps <- post.sps[, j]
            sps.idx <- (c.sps >= phi - 1.5*add.args$delta) & (c.sps <= phi+add.args$delta)
            cor.ps[j] <- mean(sps.idx)
        
        }

        #print(tns)
        #print(post.ms)
        cMTD <- which.max(cor.ps)

        #post.ms <- colMeans(post.sps)
        #cMTD <- which.min(abs(post.ms-phi))
        if (cidx > cMTD){
            cidx <- cidx - 1
        }else if (cidx == cMTD){
            cidx <- cidx
        }else {
            cidx <- cidx + 1
        }
        
        
    }
    
    
    if (earlystop==0){
        #MTD <- select.mtd(phi, tns, tys, cutoff.eli=2)$MTD
        MTD <- cMTD
    }else{
        MTD <- 99
    }
    list(MTD=MTD, dose.ns=tns, DLT.ns=tys, p.true=p.true, target=phi)
}

# target <- 0.30
# ncohort <- 10
# cohortsize <- 3
# init.level <- 1
# p.true <- c(0.1, 0.2, 0.30, 0.4, 0.5)
# print(p.true)
# add.args2 <- list(alp.prior=0.5, bet.prior=0.5, J=1e4, delta=0.05, cutoff.eli=0.95, cutoff.num=3)
# 
# MCAnew.res <- MCAABC.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, init.level=init.level,  add.args=add.args2)

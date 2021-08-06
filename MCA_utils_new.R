source("utilities.R")
library(magrittr)
library(BOIN)
library(arrApply)

sel.mtd.fn <- function(cidx, pkss, tns){ 
    ndose <- dim(pkss)[1]
    pks <- pkss[, 2]
    seq <- apply(pkss, 1,FUN=function(i)which.max(i)) 
    #print(seq)
    seq <- ascending.fn(seq)
#    print(seq)
    mask.idx <- seq == 3 & tns >= 9
    low.idx <- which(mask.idx)[1]
    if (!is.na(low.idx)){
        mask.idx[low.idx:ndose] <- TRUE
    }

    if (sum(mask.idx) == ndose){
        return(99)
    }

    # if (sum(seq==2) !=0) {
    #     mask.idx <- which(seq!=2)
    # }else if (sum(seq==1)!=0){
    #     mask.idx <- which(seq!=1)
    # }else{
    #     if (nsub < 30*ndose){
    #         mask.idx <- which(seq==3)
    #     }else{
    #         return(99)
    #     }
    # }
    pks[mask.idx] <- -1
    cMTD <- which.max(pks)
    return(cMTD)

}

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



# calculate Pr(Yn|Ak, M_n) given generated scenarios
prob.k <- function(pss, dlts, nums){
    J <- dim(pss)[1]
    numsMat <- matrix(rep(nums, J), nrow=J, byrow=T)
    dltsMat <- matrix(rep(dlts, J), nrow=J, byrow=T)
    pps <- pss ** dltsMat * (1-pss)**(numsMat-dltsMat)
    ws <- arrApply(pps, 2, "prod")
    w <- mean(ws)
    w
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


# Simulation function for MCA
MCAnew.simu.fn <- function(phi, p.true, ncohort=12, init.level=1, 
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
    pss <- lapply(1:ndose, function(k)gen.mu.rand.mul(k, J=add.args$J, K=ndose, phi=phi, delta=add.args$delta, typ=2))
    pss.under <- lapply(1:ndose, function(k)gen.mu.rand.mul(k, J=add.args$J, K=ndose, phi=phi, delta=add.args$delta, typ=1))
    pss.over <- lapply(1:ndose, function(k)gen.mu.rand.mul(k, J=add.args$J, K=ndose, phi=phi, delta=add.args$delta, typ=3))
    
    
    
    
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
        
        
        # calculate the Pr(Y_n|A_k, M_n), unnormalized
        pks <- sapply(1:ndose, function(k)prob.k(pss[[k]], tys, tns))
        pks.over <- sapply(1:ndose, function(k)prob.k(pss.over[[k]], tys, tns))
        pks.under <- sapply(1:ndose, function(k)prob.k(pss.under[[k]], tys, tns))
        #print("---------")
        #print(rbind(tys, tns))
        pkss <- cbind(pks.under, pks, pks.over)
        #print(log(pkss))
        cMTD <- sel.mtd.fn(cidx, pkss, tns)
        #print(c(which.max(pks), cMTD))
        #cMTD <- which.max(pks)
        if (cMTD  == 99){
            earlystop <- 1
            break()
        }
        if (cidx > cMTD){
            cidx <- cidx - 1
        }else if (cidx == cMTD){
            cidx <- cidx
        }else {
            cidx <- cidx + 1
        }
        
        
    }
    
    
    if (earlystop==0){
        MTD <- cMTD
        #MTD <- select.mtd(phi, tns, tys, cutoff.eli=2)$MTD
    }else{
        MTD <- 99
    }
    list(MTD=MTD, dose.ns=tns, DLT.ns=tys, p.true=p.true, target=phi)
}



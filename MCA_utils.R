# MCA2.simu.fn is the method we use in the paper

#setwd("C:/Users/Dell/Documents/ProjectCode/phaseI/Rcode")
#setwd("/Users/jinhuaqing/Documents/Projects_Code/phaseI/Rcode")
#setwd("C:/Users/JINHU/Documents/ProjectCode/MCA")
source("utilities.R")
library(magrittr)
library(BOIN)
library(arrApply)

gen.u.rand <- function(k, K=5, phi=0.3, delta=0.1){
    cps <- c(phi)
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




# Simulation function for MCA
MCA.simu.fn <- function(phi, p.true, ncohort=12, init.level=1, 
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
        
        if (overdose.fn(phi, add.args)){
            tover.doses[cidx:ndose] <- 1
        }
        
        if (tover.doses[1] == 1){
            earlystop <- 1
            break()
        }
        
        
        # calculate the Pr(Y_n|A_k, M_n), unnormalized
        pss <- lapply(1:ndose, function(k)gen.mu.rand(k, J=add.args$J, K=ndose, phi=phi, delta=add.args$delta))
        pks <- sapply(1:ndose, function(k)prob.k(pss[[k]], tys, tns))
        cMTD <- which.max(pks)
        if (cidx > cMTD){
            cidx <- cidx - 1
        }else if (cidx == cMTD){
            cidx <- min(cidx, sum(1-tover.doses))
        }else {
            cidx <- min(cidx+1, sum(1-tover.doses))
        }
        
        
    }
    
    
    if (earlystop==0){
        MTD <- select.mtd(phi, tns, tys, cutoff.eli=add.args$cutoff.eli)$MTD
    }else{
        MTD <- 99
    }
    list(MTD=MTD, dose.ns=tns, DLT.ns=tys, p.true=p.true, target=phi, over.doses=tover.doses)
}


MCA2.simu.fn <- function(phi, p.true, ncohort=12, init.level=1, 
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
        
        if (overdose.fn(phi, add.args)){
            tover.doses[cidx:ndose] <- 1
        }
        
        if (tover.doses[1] == 1){
            earlystop <- 1
            break()
        }
        
        
        # calculate the Pr(Y_n|A_k, M_n), unnormalized
        pss <- lapply(1:ndose, function(k)gen.mu.rand(k, J=add.args$J, K=ndose, phi=phi, delta=add.args$delta))
        pks <- sapply(1:ndose, function(k)prob.k(pss[[k]], tys, tns))
        cMTD <- which.max(pks)
        if (cidx > cMTD){
            cidx <- cidx - 1
        }else if (cidx == cMTD){
            cidx <- min(cidx, sum(1-tover.doses))
        }else {
            cidx <- min(cidx+1, sum(1-tover.doses))
        }
        
        
    }
    
    
    if (earlystop==0){
        # maxD <- sum(1-tover.doses)
        # pksp <- sapply(1:maxD, function(k)prob.k(pss[[k]], tys, tns))
        # MTD <- which.max(pksp)
        MTD <- min(cMTD, sum(1-tover.doses))
    }else{
        MTD <- 99
    }
    list(MTD=MTD, dose.ns=tns, DLT.ns=tys, p.true=p.true, target=phi, over.doses=tover.doses)
}

#phi <- 0.3
#p.true <- c(0.2, 0.3, 0.4, 0.5, 0.6)
#add.args <- list(alp.prior=1, bet.prior=1, J=1000, delta=0.1)
#MCA.simu.fn(phi, p.true, ncohort=10, cohortsize=3, add.args=add.args)
#MCA2.simu.fn(phi, p.true, ncohort=10, cohortsize=3, add.args=add.args)
#

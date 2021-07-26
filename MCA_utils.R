# MCA2.simu.fn is the method we use in the paper

#setwd("C:/Users/Dell/Documents/ProjectCode/phaseI/Rcode")
#setwd("/Users/jinhuaqing/Documents/Projects_Code/phaseI/Rcode")
#setwd("C:/Users/JINHU/Documents/ProjectCode/MCA")
source("utilities.R")
library(magrittr)
library(BOIN)
library(arrApply)

iso.reg <- function(y, n)
{
    ## isotonic transformation using the pool adjacent violator algorithm (PAVA)
    pava <- function (x, wt = rep(1, length(x))) 
    {
        n <- length(x)
        if (n <= 1) 
            return(x)
        if (any(is.na(x)) || any(is.na(wt))) {
            stop("Missing values in 'x' or 'wt' not allowed")
        }
        lvlsets <- (1:n)
        repeat {
            viol <- (as.vector(diff(x)) < 0)
            if (!(any(viol))) 
                break
            i <- min((1:(n - 1))[viol])
            lvl1 <- lvlsets[i]
            lvl2 <- lvlsets[i + 1]
            ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
            x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
            lvlsets[ilvl] <- lvl1
        }
        x
    }
    
  ## poster mean and variance of toxicity probabilities using beta(0.005, 0.005) as the prior 
  phat = (y+0.005)/(n+0.01); 
  phat.var = (y+0.005)*(n-y+0.005)/((n+0.01)^2*(n+0.01+1))
   
  ## perform the isotonic transformation using PAVA
  phat = pava(phat, wt=1/phat.var) 
  phat
}

mtd.sel.post.fn <- function(tys, tns, phi, delta, alp.prior, bet.prior){
   # tys: num of dlts for each dose level, vector of K  
   # tns: num of patients for each dose level, vector of K  
   # phi: the target prob
   # delta: win size
   # alp.prior, bet.prior: parameters of the prior
    
    post.alps <- tys + alp.prior
    post.bets <- tns - tys + bet.prior
    vs <- post.alps*post.bets/((post.alps+post.bets)**2)/(post.alps+post.bets+1)
    ms <- iso.reg(tys, tns)
    
    alps <- (ms**2*(1-ms)-ms*vs)/vs
    bets <- alps/ms - alps
    
    K <- length(alps)
    low <- phi - delta
    up <- phi + delta
    probs <- c()
    for (k in 1:K){
        alp <- alps[k]
        bet <- bets[k]
        cur.prob <- pbeta(up, alp, bet) - pbeta(low, alp, bet)
        probs[k] <- cur.prob
    }
    mtd <- which.max(probs)
    res <- list(mtd=mtd, probs=probs)
    res
}



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
    pss <- lapply(1:ndose, function(k)gen.mu.rand(k, J=add.args$J, K=ndose, phi=phi, delta=add.args$delta))

    
    
    
    
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
        pks <- sapply(1:ndose, function(k)prob.k(pss[[k]], tys, tns))
        cMTD <- which.max(pks)
        if (cidx > cMTD){
            cidx <- cidx - 1
        }else if (cidx == cMTD){
            cidx <- cidx
        }else {
            cidx <- cidx + 1
        }
        
        
    }
    
    
    if (earlystop==0){
        MTD <- select.mtd(phi, tns, tys, cutoff.eli=2)$MTD
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
    pss <- lapply(1:ndose, function(k)gen.mu.rand(k, J=add.args$J, K=ndose, phi=phi, delta=add.args$delta))
    
    
    
    
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
        pks <- sapply(1:ndose, function(k)prob.k(pss[[k]], tys, tns))
        cMTD <- which.max(pks)
        if (cidx > cMTD){
            cidx <- cidx - 1
        }else if (cidx == cMTD){
            cidx <- cidx
        }else {
            cidx <- cidx + 1
        }
        
        
    }
    
    
    if (earlystop==0){
        MTD <- mtd.sel.post.fn(tys, tns, phi, add.args$delta, add.args$alp.prior, add.args$bet.prior)$mtd
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

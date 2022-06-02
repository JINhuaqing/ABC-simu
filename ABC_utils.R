#setwd("C:/Users/Dell/Documents/ProjectCode/phaseI/Rcode")
#setwd("/Users/jinhuaqing/Documents/Projects_Code/phaseI/Rcode")
#setwd("C:/Users/JINHU/Documents/ProjectCode/MCA")
source("utilities.R")
library(magrittr)
library(BOIN)
library(arrApply)
library(spatstat)


gen.u.rand <- function(k, K=5, phi=0.3, delta=0.1){
    #cps <- c(phi)
    if (k==(K+1)){
            cps <- runif(K, 0, max=phi-1*delta)
    }else if (k > 0) {
        cps <- c(runif(1, phi-1*delta, phi+1*delta))
        #cps <- c(phi)
        if (k > 1){
            cps <- c(cps, runif(k-1, min=0, max=phi-1*delta))
        }
        if (k < K){
            cps <- c(cps, runif(K-k, min=phi+1*delta, max=2*phi))
        }
    }else if (k==0){
            cps <- runif(K, min=phi+1*delta, max=2*phi)
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
    pss <- lapply(0:K, function(k)gen.mu.rand(k, J=J*(1+as.numeric(k==-1)), K=K, phi=phi, delta=delta))
    pss.prior <- do.call(rbind, pss)
    #pss.prior <- t(apply(matrix(runif(K*J, 0, 2*phi), ncol=K), 1, sort))
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
    #print(mean(kp.idx))
    kp.idx
}

kpws.fn <- function(pss.prior, tys, tns, h=0.01){
    K <- length(tys)
    Num <- dim(pss.prior)[1]
    pss.prior.vec <- as.vector(t(pss.prior))
    tns.vec <- rep(tns, Num)
    tys.gen <- rbinom(length(tns.vec), tns.vec, pss.prior.vec)
    tys.vec <- rep(tys, Num)
    tns.mat <- matrix(tns.vec, ncol=K, byrow=T)
    diff.mat <- matrix(tys.vec - tys.gen, ncol=K, byrow=T)
    tns.mat[tns.mat==0] <- 0.1
    rate.diff.mat <- diff.mat / tns.mat
    if (is.null(h))
        h <- 0.01
    ws <- exp(-rowSums(rate.diff.mat**2)/h)# kernel fn exp(-x^2/h), rigirously, it should be exp(-x^2/2/h^2)
    ws 
}

# take weighted mode
weighted.mode <- function(vs, ws){
    #vs: the values
    #ws: the weights
    bins <- seq(min(vs), max(vs), length.out=100)
    binWs <- c()
    for (ix in 1:(length(bins)-1)){
        binWs <- c(binWs, sum(ws[vs>=bins[ix] & vs<bins[ix+1]]))
    }
    mIdx <- which.max(binWs)
    w.mode <- (bins[mIdx]+bins[mIdx+1])/2
    return(w.mode)
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
    ps.name <- paste0("./pssprior-ndose-", ndose, "-phi-", 100*phi, "-J-", add.args$J, "-delta-", 100*add.args$delta, ".RData")
    if (file.exists(ps.name)){
        load(ps.name)
    }else{
        pss.prior <- gen.prior(ndose, phi=phi, J=add.args$J, delta=add.args$delta)
        save(pss.prior, file=ps.name)
    }

    
    
    
    
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
        
        sel.cMTD.fn <- function(phi, pss.prior, kp.ws){

           ndose <- dim(pss.prior)[2]
           post.ms <- sapply(1:ndose, function(i)weighted.median(pss.prior[, i], w=kp.ws))
           cMTD <- which.min(abs(post.ms-phi))
           cMTD

        }
        

        kp.ws <- kpws.fn(pss.prior, tys, tns, h=add.args$h)
             
        cMTD <- sel.cMTD.fn(phi, pss.prior, kp.ws)

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
    }else{
        MTD <- 99
    }
    list(MTD=MTD, dose.ns=tns, DLT.ns=tys, p.true=p.true, target=phi)
}


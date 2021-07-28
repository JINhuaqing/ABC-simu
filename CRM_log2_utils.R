library(dfcrm)
library(statmod)
source("utilities.R")

logit.fn <- function(x){
    rv <- log(x/(1-x))
    rv
}

logistic.fn <- function(x){
    exp(x)/(1+exp(x))
}


ske2xs.fn <- function(p.prior, alp.m, bet.m){
    xs <- (logit.fn(p.prior) - alp.m)/exp(bet.m)
    xs
}

posterior <- function(alpha, beta, xs, y, d) {
    alp.v <- 2
    bet.v <- 2
    lik=1;
    for(i in 1:length(y))
    {
        logit.pi <- alpha + exp(beta)*xs[d[i]]
        pi = logistic.fn(logit.pi)
        lik = lik*pi^y[i]*(1-pi)^(1-y[i]);
    }
    return(lik*exp(-0.5*alpha*alpha/alp.v)*exp(-0.5*beta*beta/bet.v));
}

posterior.p1 <- function(alpha, beta, xs, y, d, target){
    v <- logit.fn(target) - xs[1]* exp(beta)
    fcts <- as.numeric(alpha > v)
    
    posterior(alpha, beta, xs, y, d)*fcts
}

# used to calculate the posterior mean of pi
posttoxf <- function(alpha, beta, xs, y, d, j) { 
    logit.pi <- alpha + exp(beta)*xs[j]
    pi = logistic.fn(logit.pi)
    rv <- pi*posterior(alpha, beta, xs, y, d)
    return(rv)
}

inte.fn <- function(fun, numds=50, ...){
    herm <- gauss.quad(numds, kind="hermite")
    nds <- herm$nodes
    ws <- herm$weights
    
    W.inv <- exp(nds**2)
    W.inv.mat1 <- matrix(rep(W.inv, numds), nrow=numds)
    W.inv.mat2 <- matrix(rep(W.inv, numds), nrow=numds, byrow=T)
    
    ws.mat1 <- matrix(rep(ws, numds), nrow=numds)
    ws.mat2 <- matrix(rep(ws, numds), nrow=numds, byrow=T)
    
    fv.list <- lapply(1:numds, function(i)fun(nds, nds[i], ...))
    fv.mat <- do.call(rbind, fv.list)
    
    vs <- W.inv.mat1 * W.inv.mat2 *ws.mat1 * ws.mat2*fv.mat
    sum(vs)
}

CRM.log2.simu.fn <-function(target = 0.30, ## Target toxicity pr
              p.true, p.prior, init.level=1, cohortsize=1, ncohort=12, add.args=list()){
    
    
    
    
    ndose <- length(p.true)
    if (missing(p.prior)){
        p.prior <- getprior(0.05, target, ceiling(ndose/2), ndose)
    }
    xs <- ske2xs.fn(p.prior, alp.m=0, bet.m=0)
    pts=rep(0,ndose);
    dlt=rep(0,ndose);
    tover.doses <- rep(0, ndose); #overdose vec
    pi.hat = numeric(ndose); # estimate of toxicity prob
    nstop = 0; # number of trial stopped due to high toxicity 
    t.start=Sys.time();
        
    y=NULL;  #binary outcome
    d=NULL;  #dose level
    dose.curr = init.level
    stop=0; #indicate if trial stops early
    for(i in 1:ncohort)
    {
        if (i == ncohort){
            cohortsize <- 1
        }
        
        # generate data for the new patient
        y = c(y, rbinom(cohortsize, 1, p.true[dose.curr]));
        d = c(d, rep(dose.curr, cohortsize));
        
        # num of dlt and subjs in current dose
        cdlt <- sum(y[d==dose.curr])
        cpts <- sum(d==dose.curr)
        
        
        
        # calculate posterior mean of toxicity probability at each dose leavel
        marginal=inte.fn(posterior,numds=50,xs=xs,y=y,d=d)
        for(j in 1:ndose) { pi.hat[j] = inte.fn(posttoxf, xs=xs,y=y,d=d,j=j)/marginal;}
        
        
        p.overtox <- inte.fn(posterior.p1,numds=50,xs=xs, y=y, d=d, target=target)/marginal;p.overtox	
        if (p.overtox > add.args$cutoff.eli){
            stop <- 1
            break()
        }

        diff = abs(pi.hat-target);
        cMTD = which.min(diff);
        if (dose.curr> cMTD){
            dose.curr <- dose.curr - 1
        }else if (dose.curr== cMTD){
            dose.curr <- dose.curr
        }else {
            dose.curr <- dose.curr + 1
        }
    }
    if(stop==1) 
    {
        MTD <- 99 
    }
    else 
    { 
        MTD <- cMTD
    }

    for (j in 1:ndose){
        pts[j] <- sum(d==j)
        dlt[j] <- sum(y[d==j])
    }	
        
        
    dose.ns <- pts
    dlt.ns <- dlt
    list(MTD=MTD, dose.ns=dose.ns, DLT.ns=dlt.ns, p.true=p.true, target=target, over.doses=tover.doses)
}




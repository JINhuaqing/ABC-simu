library(dfcrm)
source("utilities.R")



posterior <- function(alpha, p, y, d) {
    sigma2 = 2;
    lik=1;
    for(i in 1:length(y))
    {
        pi = p[d[i]]^(exp(alpha));
        lik = lik*pi^y[i]*(1-pi)^(1-y[i]);
    }
    return(lik*exp(-0.5*alpha*alpha/sigma2));
}

# used to calculate the posterior mean of pi
posttoxf <- function(alpha, p, y, d, j) { 
    p[j]^(exp(alpha))*posterior(alpha, p, y, d); 
}

CRM.simu.fn <-function(target = 0.30, ## Target toxicity pr
              p.true, p.prior, init.level=1, cohortsize=1, ncohort=12, add.args=list()){
    
    
    
    ndose <- length(p.true)
    if (missing(p.prior)){
        p.prior <- getprior(0.05, target, ceiling(ndose/2), ndose)
    }
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
        
        # generate data for the new patient
        y = c(y, rbinom(cohortsize, 1, p.true[dose.curr]));
        d = c(d, rep(dose.curr, cohortsize));
        
        # num of dlt and subjs in current dose
        cdlt <- sum(y[d==dose.curr])
        cpts <- sum(d==dose.curr)
        
        
        
        # calculate posterior mean of toxicity probability at each dose leavel
        marginal=integrate(posterior,lower=-Inf,upper=Inf,p.prior,y,d)$value
        for(j in 1:ndose) { pi.hat[j] = integrate(posttoxf,lower=-Inf,upper=Inf,p.prior,y,d,j)$value/marginal;}
        
        
        p.overtox <- integrate(posterior,lower=-Inf,upper=log(log(target)/log(p.prior[1])),p.prior,y,d)$value/marginal;	
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




#setwd("C:/Users/JINHU/Documents/ProjectCode/MCA")
setwd("/home/r6user2/Documents/TQ/MCA")
library(magrittr)
library(parallel)

source("utilities.R")
source("CRM_utils.R")
source("MCA_utils_real.R")

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

target <- 0.25
ncohort <- 13
cohortsize <- 3
init.level <- 1
nsimu <- 5000
seeds <- 1:nsimu

add.args <- list(alp.prior=1, bet.prior=1, J=10000, delta=0.10, cutoff.eli=0.95, cutoff.num=3)
p.trues <- list()
yns <- c(24, 10, 3)
tns <- c(3, 4, 2)
alp <- 1
bet <- 1
p.trues[[1]] <- (tns+alp)/(yns+bet+alp)


idx <- 1
p.true <- p.trues[[idx]]
tmtd <- MTD.level(target, p.true)


run.fn <- function(i){
    print(i)
    set.seed(seeds[i])
    MCA.res <- MCA.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, init.level=init.level,  add.args=add.args)
    CRM.res <- CRM.simu.fn(target=target, p.true=p.true, init.level=init.level, cohortsize=cohortsize, ncohort=ncohort, add.args=add.args)
    ress <- list(
                 MCA = MCA.res,
                 #MCA2 = MCA2.res,
                 CRM = CRM.res, 
                 paras=list(p.true=p.true, 
                             mtd=tmtd, 
                             add.args=add.args,
                             target=target,
                             ncohort=ncohort,
                             cohortsize=cohortsize)
        )
    ress
    
}

ncores <- 40
m.names <- c("MCA", "CRM")
results <- mclapply(1:nsimu, run.fn, mc.cores=ncores)
file.name <- paste0("./results/", "RealDataMCA_NoELiLJ", 100*add.args$cutoff.eli, "_", nsimu, "_ncohort_", ncohort, ".RData")
save(results, file=file.name)


sum.all <- list()
for (m.name in m.names){
   sum.all[[m.name]] <-  phase1.post.fn(lapply(1:nsimu, function(i)results[[i]][[m.name]]))
}
print(tmtd)
phase.I.pretty.tb(sum.all)


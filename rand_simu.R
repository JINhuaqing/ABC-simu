#setwd("C:/Users/JINHU/Documents/ProjectCode/MCA")
setwd("/home/r6user2/Documents/TQ/MCA")
library(magrittr)
library(parallel)

source("utilities.R")
source("CRM_utils.R")
source("intv_utils.R")
source("MCA_utils_ABC.R")


target <- 0.3
ncohort <- 10
cohortsize <- 3
init.level <- 1
ndose <- 5

#add.args <- list(alp.prior=0.5, bet.prior=0.5, J=1000, delta=0.05, cutoff.eli=0.95, cutoff.num=3)
add.args <- list(alp.prior=0.5, bet.prior=0.5, J=2e4, delta=0.10, cutoff.eli=0.95, cutoff.num=3, h=0.01)
nsimu <- 5000
seeds <- 1:nsimu

## Target = 0.3
# dose 3, mu1=0.55, mu2=0.40, 0.1
# dose 3, mu1=mu2=0.30, 0.07
# dose 3, mu1=mu2=0.46, 0.1
# dose 3, mu1=mu2=0.64, 0.15

# dose 5, mu1=0.60, mu2=0.50, 0.1
# dose 5, mu1=mu2=0.23, 0.05
# dose 5, mu1=mu2=0.38, 0.07
# dose 5, mu1=mu2=0.53, 0.1
# dose 5, mu1=mu2=0.71, 0.15

# dose 7, mu1=0.70, mu2=0.50, 0.1
# dose 7, mu1=mu2=0.42, 0.07
# dose 7, mu1=mu2=0.56, 0.1
# dose 7, mu1=mu2=0.74, 0.15

# Target = 0.25
# dose 8, mu1=mu2=0.23, 0.05
# dose 8, mu1=mu2=0.41, 0.07
# dose 8, mu1=mu2=0.58, 0.1
# dose 8, mu1=mu2=0.77, 0.15

set.seed(1)
ps.name <- paste0("./pssprior-ndose-", ndose, "-phi-", 100*target, "-J-", add.args$J, "-delta-", 100*add.args$delta, ".RData")
if (F){
#if (file.exists(ps.name)){
        load(ps.name)
}else{
        pss.prior <- gen.prior(ndose, phi=target, J=add.args$J, delta=add.args$delta)
        save(pss.prior, file=ps.name)
}

mu <- 0.23
Delta <- 0.05
mus <- c(0.23, 0.38, 0.53, 0.71)
Deltas <- c(0.05, 0.07, 0.10, 0.15)
for (jj in 1:4){
    mu <- mus[jj]
    Delta <- Deltas[jj]
    ndose <- ndose
    run.fn <- function(k){
        print(k)
        set.seed(seeds[k])
        p.true.all <- gen.rand.doses(ndose, target, mu1=mu, mu2=mu)
        p.true <- p.true.all$p.true
        tmtd <- p.true.all$mtd.level
        #print(p.true)
    
        CRM.res <- CRM.simu.fn(target=target, p.true=p.true, init.level=init.level, cohortsize=cohortsize, ncohort=ncohort, add.args=add.args)
       #(1--CCD, 2--mTPI, 3--BOIN, 4--Keyboard, 5--UMPBI) \n")
        MCAnew.res <- MCAABC.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, init.level=init.level,  add.args=add.args)
        
        CCD.res   <- intv.simu.fn(target=target, p.true=p.true, ncohort=ncohort,  cutoff.eli=add.args$cutoff.eli, init.level=init.level, cohortsize=cohortsize, design=1)
        mTPI.res  <- intv.simu.fn(target=target, p.true=p.true, ncohort=ncohort,  cutoff.eli=add.args$cutoff.eli, init.level=init.level, cohortsize=cohortsize, design=2)
        BOIN.res  <- intv.simu.fn(target=target, p.true=p.true, ncohort=ncohort,  cutoff.eli=add.args$cutoff.eli, init.level=init.level, cohortsize=cohortsize, design=3)
        keyB.res  <- intv.simu.fn(target=target, p.true=p.true, ncohort=ncohort,  cutoff.eli=add.args$cutoff.eli, init.level=init.level, cohortsize=cohortsize, design=4)
        UMPBI.res <- intv.simu.fn(target=target, p.true=p.true, ncohort=ncohort,  cutoff.eli=add.args$cutoff.eli, init.level=init.level, cohortsize=cohortsize, design=5)
        ress <- list(
                     MCAnew = MCAnew.res,
                     BOIN = BOIN.res, 
                     CCD = CCD.res, 
                     keyB = keyB.res, 
                     mTPI= mTPI.res, 
                     UMPBI= UMPBI.res, 
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
    
    
    file.name <- paste0("./results/", "SimuMCA_ABC_NoEliLJ", 100*add.args$cutoff.eli, "_", nsimu, "_ncohort_", ncohort, "_random_", Delta,  "_priorDelta_", 100*add.args$delta, ".RData")
    results <- mclapply(1:nsimu, run.fn, mc.cores=50)
    save(results, file=file.name)
    print(post.process.random(results))
}



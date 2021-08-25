#setwd("C:/Users/JINHU/Documents/ProjectCode/MCA")
setwd("/home/r6user2/Documents/TQ/MCA")
library(magrittr)
library(parallel)

source("utilities.R")
source("CRM_utils.R")
#source("MCA_utils.R")
source("MCA_utils_ABC.R")
source("intv_utils.R")


target <- 0.25
ncohort <- 16
cohortsize <- 3
init.level <- 1
nsimu <- 5000
seeds <- 1:nsimu

add.args <- list(alp.prior=0.5, bet.prior=0.5, J=2e4, delta=0.10, cutoff.eli=0.95, cutoff.num=3)
p.trues <- list()
p.trues[[1]] <- c(0.05, 0.25, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)
p.trues[[2]] <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.25, 0.50, 0.60)
p.trues[[3]] <- c(0.01, 0.05, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)
p.trues[[4]] <- c(0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99)
p.trues[[5]] <- c(0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85)
p.trues[[6]] <- c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75)


idx <- 4
for (idx in 1:6){

p.true <- p.trues[[idx]]
tmtd <- MTD.level(target, p.true)
print(p.true)

ndose <- length(p.true)
set.seed(1)
ps.name <- paste0("./pssprior-ndose-", ndose, "-phi-", 100*target, "-J-", add.args$J, "-delta-", 100*add.args$delta, ".RData")
if (F){
#if (file.exists(ps.name)){
        load(ps.name)
}else{
        pss.prior <- gen.prior(ndose, phi=target, J=add.args$J, delta=add.args$delta)
        save(pss.prior, file=ps.name)
}


run.fn <- function(i){
    print(c(idx, i))
    #set.seed(seeds[i])
    MCAnew.res <- MCAABC.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, init.level=init.level,  add.args=add.args)
    CRM.res <- CRM.simu.fn(target=target, p.true=p.true, init.level=init.level, cohortsize=cohortsize, ncohort=ncohort, add.args=add.args)
        CCD.res   <- intv.simu.fn(target=target, p.true=p.true, ncohort=ncohort,  cutoff.eli=add.args$cutoff.eli, init.level=init.level, cohortsize=cohortsize, design=1)
        mTPI.res  <- intv.simu.fn(target=target, p.true=p.true, ncohort=ncohort,  cutoff.eli=add.args$cutoff.eli, init.level=init.level, cohortsize=cohortsize, design=2)
        BOIN.res  <- intv.simu.fn(target=target, p.true=p.true, ncohort=ncohort,  cutoff.eli=add.args$cutoff.eli, init.level=init.level, cohortsize=cohortsize, design=3)
        keyB.res  <- intv.simu.fn(target=target, p.true=p.true, ncohort=ncohort,  cutoff.eli=add.args$cutoff.eli, init.level=init.level, cohortsize=cohortsize, design=4)
        UMPBI.res <- intv.simu.fn(target=target, p.true=p.true, ncohort=ncohort,  cutoff.eli=add.args$cutoff.eli, init.level=init.level, cohortsize=cohortsize, design=5)
    
    ress <- list(
                 MCAnew = MCAnew.res,
                 CRM = CRM.res, 
                 CCD = CCD.res, 
                 mTPI= mTPI.res, 
                 BOIN = BOIN.res, 
                 keyB = keyB.res, 
                 UMPBI= UMPBI.res, 
                 paras=list(p.true=p.true, 
                             mtd=tmtd, 
                             add.args=add.args,
                             target=target,
                             ncohort=ncohort,
                             cohortsize=cohortsize)
        )
    ress
    
}

ncores <- 50
m.names <- c("MCAnew", "BOIN", "CCD", "CRM", "keyB", "mTPI","UMPBI")
results <- mclapply(1:nsimu, run.fn, mc.cores=ncores)
file.name <- paste0("./results/", "SimuMCA_NoEliLJ", 100*add.args$cutoff.eli, "_", nsimu, "_ncohort_", ncohort, "_fix2_",  idx, ".RData")
save(results, file=file.name)


sum.all <- list()
for (m.name in m.names){
   sum.all[[m.name]] <-  phase1.post.fn(lapply(1:nsimu, function(i)results[[i]][[m.name]]))
}
print(tmtd)
print(phase.I.pretty.tb(sum.all))
}

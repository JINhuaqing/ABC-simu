#setwd("C:/Users/JINHU/Documents/ProjectCode/MCA")
setwd("/home/r6user2/Documents/TQ/MCA")
library(magrittr)
library(parallel)

source("utilities.R")
source("CRM_utils.R")
source("MCA_utils.R")
source("intv_utils.R")


#set.seed(2021, kind = "L'Ecuyer-CMRG") 
target <- 0.20
ncohort <- 12
cohortsize <- 3
init.level <- 1
nsimu <- 5000
seeds <- 1:nsimu

add.args <- list(alp.prior=1, bet.prior=1, J=1000, delta=0.15, cutoff.eli=0.85, cutoff.num=3)
p.trues <- list()
p.trues[[1]] <- c(0.05, 0.10, 0.20, 0.30, 0.50, 0.70)
p.trues[[2]] <- c(0.30, 0.40, 0.52, 0.61, 0.76, 0.87)
p.trues[[3]] <- c(0.05, 0.06, 0.08, 0.11, 0.19, 0.34)
p.trues[[4]] <- c(0.06, 0.08, 0.12, 0.18, 0.40, 0.71)
p.trues[[5]] <- c(0.00, 0.00, 0.03, 0.05, 0.11, 0.22)

idx <- 5
p.true <- p.trues[[idx]]
tmtd <- MTD.level(target, p.true)


run.fn <- function(i){
    print(i)
    set.seed(seeds[i])
    MCA2.res <- MCA2.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, init.level=init.level,  add.args=add.args)
    CRM.res <- CRM.simu.fn(target=target, p.true=p.true, init.level=init.level, cohortsize=cohortsize, ncohort=ncohort, add.args=add.args)
    #MCA.res <- MCA.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, init.level=init.level,  add.args=add.args)
   #(1--CCD, 2--mTPI, 3--BOIN, 4--Keyboard, 5--UMPBI) \n")
    CCD.res   <- intv.simu.fn(target=target, p.true=p.true, ncohort=ncohort,  init.level=init.level, cohortsize=cohortsize, design=1, cutoff.eli=add.args$cutoff.eli)
    mTPI.res  <- intv.simu.fn(target=target, p.true=p.true, ncohort=ncohort,  init.level=init.level, cohortsize=cohortsize, design=2, cutoff.eli=add.args$cutoff.eli)
    BOIN.res  <- intv.simu.fn(target=target, p.true=p.true, ncohort=ncohort,  init.level=init.level, cohortsize=cohortsize, design=3, cutoff.eli=add.args$cutoff.eli)
    keyB.res  <- intv.simu.fn(target=target, p.true=p.true, ncohort=ncohort,  init.level=init.level, cohortsize=cohortsize, design=4, cutoff.eli=add.args$cutoff.eli)
    UMPBI.res <- intv.simu.fn(target=target, p.true=p.true, ncohort=ncohort,  init.level=init.level, cohortsize=cohortsize, design=5, cutoff.eli=add.args$cutoff.eli)
    
    ress <- list(
    #             MCA = MCA.res,
                 MCA2 = MCA2.res,
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

ncores <- 30
m.names <- c("MCA2", "BOIN", "CCD", "CRM", "keyB", "mTPI","UMPBI")
#m.names <- c("MCA", "MCA2", "BOIN", "CCD", "CRM", "keyB", "mTPI","UMPBI")
results <- mclapply(1:nsimu, run.fn, mc.cores=ncores)
file.name <- paste0("./results/", "Simu", 100*add.args$cutoff.eli, "_", nsimu, "_ncohort_", ncohort, "_fix1_",  idx, ".RData")
save(results, file=file.name)


sum.all <- list()
for (m.name in m.names){
   sum.all[[m.name]] <-  phase1.post.fn(lapply(1:nsimu, function(i)results[[i]][[m.name]]))
}
print(tmtd)
phase.I.pretty.tb(sum.all)


#setwd("C:/Users/JINHU/Documents/ProjectCode/MCA")
setwd("/home/r6user2/Documents/TQ/MCA")
library(magrittr)
library(parallel)

source("utilities.R")
source("CRM_utils.R")
source("MCA_utils.R")
source("MCA_utils_new.R")


target <- 0.25
ncohort <- 16
cohortsize <- 3
init.level <- 1
nsimu <- 1000
seeds <- 1:nsimu

add.args <- list(alp.prior=1, bet.prior=1, J=1000, delta=0.05, cutoff.eli=0.95, cutoff.num=3)
p.trues <- list()
p.trues[[1]] <- c(0.05, 0.25, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)
p.trues[[2]] <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.25, 0.50, 0.60)
p.trues[[3]] <- c(0.01, 0.05, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95)
p.trues[[4]] <- c(0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 0.99)
p.trues[[5]] <- c(0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85)
p.trues[[6]] <- c(0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75)


idx <- 6
p.true <- p.trues[[idx]]
tmtd <- MTD.level(target, p.true)
print(p.true)


run.fn <- function(i){
    print(i)
    #set.seed(seeds[i])
    MCA.res <- MCA.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, init.level=init.level,  add.args=add.args)
        MCAnew.res <- MCAnew.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, init.level=init.level,  add.args=add.args)
    CRM.res <- CRM.simu.fn(target=target, p.true=p.true, init.level=init.level, cohortsize=cohortsize, ncohort=ncohort, add.args=add.args)
    
    ress <- list(
                 MCA = MCA.res,
                 MCAnew = MCAnew.res,
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

ncores <- 70
m.names <- c("MCA", "MCAnew", "CRM")
results <- mclapply(1:nsimu, run.fn, mc.cores=ncores)
#file.name <- paste0("./results/", "SimuMCA_NoEliLJ", 100*add.args$cutoff.eli, "_", nsimu, "_ncohort_", ncohort, "_fix2_",  idx, ".RData")
#save(results, file=file.name)


sum.all <- list()
for (m.name in m.names){
   sum.all[[m.name]] <-  phase1.post.fn(lapply(1:nsimu, function(i)results[[i]][[m.name]]))
}
print(tmtd)
print(phase.I.pretty.tb(sum.all))

#setwd("C:/Users/JINHU/Documents/ProjectCode/MCA")
setwd("/home/r6user2/Documents/TQ/MCA")
library(magrittr)
library(parallel)

source("utilities.R")
source("CRM_log2_utils.R")
source("MCA_utils_real.R")


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
alp <- 0
bet <- 0
p.trues[[1]] <- (tns+alp)/(yns+bet+alp)


idx <- 1
p.true <- p.trues[[idx]]
tmtd <- MTD.level(target, p.true)


run.fn <- function(i){
    print(i)
    set.seed(seeds[i])
    MCA.res <- MCA.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, init.level=init.level,  add.args=add.args)
    CRM.res <- CRM.log2.simu.fn(target=target, p.true=p.true, init.level=init.level, cohortsize=cohortsize, ncohort=ncohort, add.args=add.args)
    ress <- list(
                 MCA = MCA.res,
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

ncores <- 50
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


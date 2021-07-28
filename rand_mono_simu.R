#setwd("C:/Users/JINHU/Documents/ProjectCode/MCA")
setwd("/home/r6user2/Documents/TQ/MCA")
library(magrittr)
library(parallel)

source("utilities.R")
source("MCA_utils_no_mono.R")
source("MCA_utils.R")


target <- 0.25
ncohort <- 16
cohortsize <- 3
init.level <- 1

add.args <- list(alp.prior=1, bet.prior=1, J=1000, delta=0.1, cutoff.eli=0.95, cutoff.num=3)
nsimu <- 5000
seeds <- 1:nsimu


mu <- 0.23
Delta <- 0.05
mus <- c(0.23, 0.41, 0.58, 0.77)
Deltas <- c(0.05, 0.07, 0.10, 0.15)
for (jj in 4){
    mu <- mus[jj]
    Delta <- Deltas[jj]
    ndose <- 8
    run.fn <- function(k){
        print(k)
        set.seed(seeds[k])
        p.true.all <- gen.rand.doses(ndose, target, mu1=mu, mu2=mu)
        p.true <- p.true.all$p.true
        tmtd <- p.true.all$mtd.level
    
        MCA.res <- MCA.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, init.level=init.level,  add.args=add.args)
        MCAno.res <- MCAno.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, init.level=init.level,  add.args=add.args)
        ress <- list(
                     MCA = MCA.res,
                     MCAno = MCAno.res, 
                     paras=list(p.true=p.true, 
                                 mtd=tmtd, 
                                 add.args=add.args,
                                 target=target,
                                 ncohort=ncohort,
                                 cohortsize=cohortsize)
            )
        ress
    }
    
    
    file.name <- paste0("./results/", "SimuNomono", 100*add.args$cutoff.eli, "_", nsimu, "_ncohort_", ncohort, "_random_", Delta,  ".RData")
    results <- mclapply(1:nsimu, run.fn, mc.cores=40)
    print(post.process.random(results))
    save(results, file=file.name)
}



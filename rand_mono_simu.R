#setwd("C:/Users/JINHU/Documents/ProjectCode/MCA")
setwd("/home/r6user2/Documents/TQ/MCA")
library(magrittr)
library(parallel)

source("utilities.R")
source("MCA_utils_ABC_no_mono.R")
source("MCA_utils_ABC.R")


target <- 0.3
ncohort <- 10
cohortsize <- 3
init.level <- 1

add.args <- list(alp.prior=0.5, bet.prior=0.5, J=20000, delta=0.1, cutoff.eli=0.95, cutoff.num=3, h=0.01)
nsimu <- 5000
seeds <- 1:nsimu


mu <- 0.23
Delta <- 0.05
mus <- c(0.23, 0.38, 0.53, 0.71)
Deltas <- c(0.05, 0.07, 0.10, 0.15)
for (jj in 4){
    mu <- mus[jj]
    Delta <- Deltas[jj]
    ndose <- 5
    run.fn <- function(k){
        print(k)
        set.seed(seeds[k])
        p.true.all <- gen.rand.doses(ndose, target, mu1=mu, mu2=mu)
        p.true <- p.true.all$p.true
        tmtd <- p.true.all$mtd.level
    
        MCA.res <- MCAABC.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, init.level=init.level,  add.args=add.args)
        MCAno.res <- MCAABCno.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, init.level=init.level,  add.args=add.args)
        ress <- list(
                     MCAnew = MCA.res,
                     MCAnewno = MCAno.res, 
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
    results <- mclapply(1:nsimu, run.fn, mc.cores=80)
    print(post.process.random(results))
    save(results, file=file.name)
}



#setwd("C:/Users/JINHU/Documents/ProjectCode/MCA")
setwd("/home/huaqingj/MyResearch/ABC-simu")
library(magrittr)
library(parallel)

source("utilities.R")
source("ABC_utils.R")


target <- 0.2
mus <- c(0.33, 0.46, 0.63, 0.83) #0.2
#mus <- c(0.23, 0.38, 0.53, 0.71) #0.3
ncohort <- 10
cohortsize <- 3
init.level <- 1
ndose <- 5

#add.args <- list(alp.prior=0.5, bet.prior=0.5, J=1000, delta=0.05, cutoff.eli=0.95, cutoff.num=3)
add.args <- list(alp.prior=0.5, bet.prior=0.5, J=2e4, delta=0.10, cutoff.eli=0.95, cutoff.num=3, h=0.01, type="median")
add.args.mean <- list(alp.prior=0.5, bet.prior=0.5, J=2e4, delta=0.10, cutoff.eli=0.95, cutoff.num=3, h=0.01, type="mean")
add.args.mode <- list(alp.prior=0.5, bet.prior=0.5, J=2e4, delta=0.10, cutoff.eli=0.95, cutoff.num=3, h=0.01, type="mode")
nsimu <- 5000
seeds <- 1:nsimu

## Target = 0.3
# dose 5, mu1=mu2=0.23, 0.05
# dose 5, mu1=mu2=0.38, 0.07
# dose 5, mu1=mu2=0.53, 0.1
# dose 5, mu1=mu2=0.71, 0.15


# Target = 0.20
# dose 5, mu1=mu2=0.33, 0.05
# dose 5, mu1=mu2=0.46, 0.07
# dose 5, mu1=mu2=0.63, 0.1
# dose 5, mu1=mu2=0.83, 0.15

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
Deltas <- c(0.05, 0.07, 0.10, 0.15)
for (jj in 3){
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
    
       #(1--CCD, 2--mTPI, 3--BOIN, 4--Keyboard, 5--UMPBI) \n")
        Median.res <- MCAABC.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, init.level=init.level,  add.args=add.args)
        Mean.res <- MCAABC.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, init.level=init.level,  add.args=add.args.mean)
        Mode.res <- MCAABC.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, init.level=init.level,  add.args=add.args.mode)
        
        ress <- list(
                     median = Median.res,
                     mean = Mean.res,
                     mode = Mode.res,
                     paras=list(p.true=p.true, 
                                 mtd=tmtd, 
                                 add.args=add.args,
                                 target=target,
                                 ncohort=ncohort,
                                 cohortsize=cohortsize)
            )
        ress
    }
    
    
    file.name <- paste0("./results/", "SimuMCA_mmmABC_NoEliLJ", 100*add.args$cutoff.eli, "_", nsimu, "_ncohort_", ncohort, "_random_", Delta,  "_priorDelta_", 100*add.args$delta, "_target_", 100*target, ".RData")
    print(file.name)
    results <- mclapply(1:nsimu, run.fn, mc.cores=10)
    save(results, file=file.name)
    print(post.process.random(results))
}



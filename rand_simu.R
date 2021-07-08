setwd("../")
library(magrittr)
library(parallel)
source("utilities.R")
source("./phaseI/crm_utils.R")
source("./phaseI/boin_utils.R")
source("ORM_utils.R")


#set.seed(10)
target <- 0.3
ncohort <- 10
cohortsize <- 3

add.args <- list(alp.prior=target, bet.prior=1-target)

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

mu <- 0.23
run.fn <- function(k){
    print(k)
    p.true.all <- gen.rand.doses(5, target, mu1=mu, mu2=mu)
    p.true <- p.true.all$p.true
    tmtd <- p.true.all$mtd.level

    orm.res <- ORM.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, add.args=add.args)
    crm.res <- crm.simu.fn(target=target, p.true=p.true, cohortsize=cohortsize, ncohort=ncohort)
    boin.res <- boin.simu.fn(target=target, p.true=p.true, ncohort=ncohort, cohortsize)

    ress <- list(
                 orm=orm.res,
                 crm = crm.res, 
                 boin = boin.res, 
                 paras=list(p.true=p.true, 
                         mtd=tmtd, 
                         add.args=add.args,
                         target=target,
                         ncohort=ncohort,
                         cohortsize=cohortsize)
                 )
    ress
}


nsimu <- 5000
file.name <- paste0("../results/", "MTDSimu_", nsimu, "random_0.05",  ".RData")
results <- mclapply(1:nsimu, run.fn, mc.cores=20)
post.process.random(results)
save(results, file=file.name)


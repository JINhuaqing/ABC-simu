setwd("../")
library(magrittr)
library(parallel)
source("utilities.R")
source("./phaseI/crm_utils.R")
source("./phaseI/boin_utils.R")
source("./phaseI/anova_settings.R")
source("ORM_utils.R")


nsimu <- 1000

# R code to find mu for random scenario
#target <- 0.25
#can.mus <- seq(0.1, 0.85, 0.01)
#means <- c()
#for (mu in can.mus){
#    ress <- lapply(1:10000, function(i)gen.rand.doses(9, target, mu1=mu, mu2=mu))
#    mean.v <- sapply(ress, prob.diff.fn, target=target) %>% mean
#    means <- c(means, mean.v)
#}
#print(paste(can.mus, means))
#fasdf

mus.list <- list()
# target 0.25, 0.30, 0.33
mus.list[[1]] <- list()
mus.list[[2]] <- list()
mus.list[[3]] <- list()

# target 0.25
# dose 3, 5, 7, 9
mus.list[[1]][[1]] <- c(0.18, 0.34, 0.49, 0.68)
mus.list[[1]][[2]] <- c(0.27, 0.42, 0.57, 0.76)
mus.list[[1]][[3]] <- c(0.31, 0.45, 0.60, 0.79)
mus.list[[1]][[4]] <- c(0.33, 0.47, 0.62, 0.81)

# target 0.30
# dose 3, 5, 7, 9
mus.list[[2]][[1]] <- c(0.13, 0.31, 0.46, 0.64)
mus.list[[2]][[2]] <- c(0.23, 0.38, 0.53, 0.71)
mus.list[[2]][[3]] <- c(0.28, 0.42, 0.56, 0.74)
mus.list[[2]][[4]] <- c(0.29, 0.43, 0.57, 0.76)

# target 0.33
# dose 3, 5, 7, 9
mus.list[[3]][[1]] <- c(0.09, 0.29, 0.44, 0.62)
mus.list[[3]][[2]] <- c(0.21, 0.37, 0.51, 0.69)
mus.list[[3]][[3]] <- c(0.25, 0.40, 0.54, 0.72)
mus.list[[3]][[4]] <- c(0.27, 0.42, 0.56, 0.74)

flag <- 0
for (i1 in 4:length(nlevels)){
    for (i2 in 1:length(csizes)){
        for (i3 in c(2, 6, 8, 9)){ # sample size
            for (i4 in 1:length(diff.probs)){
                for (i5 in 1:length(targets)){
                    flag <- flag + 1
                    nlevel <- nlevels[i1]
                    cohortsize <- csizes[i2]
                    cur.mu <- mus.list[[i5]][[i1]][i4]
                    sample.size <- sample.sizes[i3]
                    target <- targets[i5]
                    ncohort <- sample.size/cohortsize
                    paras <- c(nlevel, cohortsize, cur.mu, sample.size, target, ncohort)
                    names(paras) <- c("nlevel", "chortesize", "mu", "sample size", "target", "num of cohort")
                    print(paras)
                    add.args <- list(alp.prior=target, bet.prior=1-target)

                    run.fn <- function(k){
                        print(c(flag, k))
                        p.true.all <- gen.rand.doses(nlevel, target, mu1=cur.mu, mu2=cur.mu)
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
                    
                    
                    file.name <- paste0("../results/", "Anova_MTDSimu_", i1, i2, i3, i4, i5, "_", nsimu, ".RData")
                    results <- mclapply(1:nsimu, run.fn, mc.cores=20)
                    #post.process.random(results)
                    save(results, file=file.name)
                }
            }
        }
    }
}


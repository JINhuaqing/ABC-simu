setwd("/home/r6user2/Documents/TQ/MCA")
library(magrittr)
library(parallel)

source("utilities.R")
source("MCA_utils_ABC.R")
#source("MCA_utils.R")
source("anova_settings.R")


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

cohortsize <- 3
flag <- 0
for (i1 in 3:length(nlevels)){
    for (i2 in 5:length(deltas)){
        for (i3 in 9:length(sample.sizes)){ # sample size
            for (i4 in 1:length(diff.probs)){
                for (i5 in 1:length(targets)){
                    for (i6 in 1:length(hs)){
                    flag <- flag + 1
                    nlevel <- nlevels[i1]
                    delta <- deltas[i2]
                    cur.mu <- mus.list[[i5]][[i1]][i4]
                    sample.size <- sample.sizes[i3]
                    target <- targets[i5]
                    h <- hs[i6]
                    ncohort <- sample.size/cohortsize
                    paras <- c(nlevel, cohortsize, cur.mu, sample.size, target, ncohort, delta, h)
                    names(paras) <- c("nlevel", "chortesize", "mu", "sample size", "target", "num of cohort", "delta", "h")
                    print(paras)

                    run.fn <- function(k){
                        print(c(flag, k))
                        p.true.all <- gen.rand.doses(nlevel, target, mu1=cur.mu, mu2=cur.mu)
                        p.true <- p.true.all$p.true
                        tmtd <- p.true.all$mtd.level
                        if (delta == 1){
                            cdelta <- round(runif(1, 0, 0.2), 2)
                        }else{
                            cdelta <- delta
                        }
                        add.args <- list(alp.prior=0.5, bet.prior=0.5, J=2e4, delta=cdelta, cutoff.eli=0.95, cutoff.num=3, h=h)

          ps.name <- paste0("./pssprior-ndose-", nlevel, "-phi-", 100*target, "-J-", add.args$J, "-delta-", 100*add.args$delta, ".RData")
          if (file.exists(ps.name)){
                  load(ps.name)
          }else{
                  pss.prior <- gen.prior(nlevel, phi=target, J=add.args$J, delta=add.args$delta)
                  save(pss.prior, file=ps.name)
          }
                        MCA.res <- MCAABC.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, init.level=1,  add.args=add.args)
                    
                    
                        ress <- list(
                                     MCA = MCA.res, 
                                     paras=list(p.true=p.true, 
                                             mtd=tmtd, 
                                             add.args=add.args,
                                             target=target,
                                             ncohort=ncohort,
                                             cohortsize=cohortsize)
                                     )
                        ress
                    }
                    
                    
                    file.name <- paste0("./results/", "AnovaNoEli_Simu_", i1, i2, i3, i4, i5, i6, "_", nsimu, ".RData")
                    results <- mclapply(1:nsimu, run.fn, mc.cores=60)
                    print(post.process.random(results))
                    save(results, file=file.name)
                  }
                }
            }
        }
    }
}


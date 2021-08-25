setwd("/home/r6user2/Documents/TQ/MCA")
library(magrittr)
library(parallel)

source("utilities.R")
source("MCA_utils_ABC.R")
#source("MCA_utils.R")
source("anova_settings.R")


nsimu <- 1000
deltas <- seq(0, 0.2, 0.01)


flag <- 0
for (i1 in 2:length(nlevels)){
    for (i2 in 1:length(deltas)){
       for (i3 in 1:length(targets)){
        flag <- flag + 1
        print(c(flag, i1, i2, i3))
        nlevel <- nlevels[i1]
        delta <- deltas[i2]
        cdelta <- delta
        target <- targets[i3]
        add.args <- list(alp.prior=0.5, bet.prior=0.5, J=2e4, delta=cdelta, cutoff.eli=0.95, cutoff.num=3)

        ps.name <- paste0("./pssprior-ndose-", nlevel, "-phi-", 100*target, "-J-", add.args$J, "-delta-", 100*add.args$delta, ".RData")
        if (F){
                  load(ps.name)
         }else{
                  pss.prior <- gen.prior(nlevel, phi=target, J=add.args$J, delta=add.args$delta)
                  save(pss.prior, file=ps.name)
         }
                    
                    
        }
    }
}

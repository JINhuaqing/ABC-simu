setwd("C:/Users/JINHU/Documents/ProjectCode/MCA")
source("utilities.R")
source("MCA_utils_ABC_real.R")


# the results for 5000 repititions
load("./results/RealDataMCA_NoELiLJ95_5000_ncohort_13.RData")
m.names <- c("MCAnew", "CRM")
nsimu <- 5000
sum.all <- list()
for (m.name in m.names){
   sum.all[[m.name]] <-  phase1.post.fn(lapply(1:nsimu, function(i)results[[i]][[m.name]]))
}
phase.I.pretty.tb(sum.all)


# select examples
target <- 0.25
ncohort <- 13
cohortsize <- 3
init.level <- 1
nsimu <- 5000
add.args <- list(alp.prior=0.5, bet.prior=0.5, J=2e4, delta=0.10, cutoff.eli=0.95, cutoff.num=3, h=0.01)
p.true <- c(0.125, 0.40, 0.667)
tmtd <- MTD.level(target, p.true)
ndose <- length(p.true)

set.seed(1) 
ps.name <- paste0("./pssprior-ndose-", ndose, "-phi-", 100*target, "-J-", add.args$J, "-delta-", 100*add.args$delta, ".RData")
if (F){
    load(ps.name)
}else{
    pss.prior <- gen.prior(ndose, phi=target, J=add.args$J, delta=add.args$delta)
    save(pss.prior, file=ps.name)
}



set.seed(10)
ex.res <- MCAABCreal.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, init.level=init.level,  add.args=add.args)
ex.res$MTD
pksMat <- do.call(rbind, ex.res$tk.pkss)
resMat <- cbind(pksMat, ex.res$tk.dlts, ex.res$tk.idxs, c(rep(3, 12), 1), ex.res$tk.over)
colnames(resMat) <- c(paste("Dose Level", 1:3), "num DLT", "CurDose idx", "num Patients", "OverDose")
round(resMat, 2)

# DLTs:
ex.res$tk.dlts
# dose level 
ex.res$tk.idxs
# patient number:
c(rep(3, 12), 1)

# DLTs:        0 2 0 1 2 1 1 0 0 1 0 0 0
# dose level : 1 2 1 2 2 1 1 1 1 1 1 1 1
# Subj number: 3 3 3 3 3 3 3 3 3 3 3 3 1

setwd("C:/Users/JINHU/Documents/ProjectCode/MCA")
source("utilities.R")
source("MCA_utils_real.R")


# the results for 5000 repititions
load("./results/finalRes/RealData95_5000_ncohort_13_fix3_1.RData")
m.names <- c("MCA", "MCA2", "CRM")
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
add.args <- list(alp.prior=1, bet.prior=1, J=1000, delta=0.10, cutoff.eli=0.95, cutoff.num=3)
p.true <- c(0.12, 0.40, 0.67)
tmtd <- MTD.level(target, p.true)

set.seed(10) 
ex.res <- MCA2.simu.fn(target, p.true, ncohort=ncohort, cohortsize=cohortsize, init.level=init.level,  add.args=add.args)
ex.res$MTD
ex.res$tk.idxs
ex.res$tk.dlts
ex.res$tk.pkss
ex.res$over.doses
ex.res$tk.over
pksMat <- do.call(rbind, ex.res$tk.pkss)
resMat <- cbind(pksMat, ex.res$tk.dlts, ex.res$tk.idxs, ex.res$tk.over)
colnames(resMat) <- c(paste("Dose Level", 1:3), "DLT", "CurDose", "OverDose")
resMat


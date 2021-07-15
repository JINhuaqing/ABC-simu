setwd("/home/r6user2/Documents/TQ/MCA")
source("utilities.R")


cutoff.eli <- 0.95
nsimu <- 5000
ncohort <- 16
idx <- 1

file.name <- paste0("./results/", "Simu", 100*cutoff.eli, "_", nsimu, "_ncohort_", ncohort, "_fix2_",  idx, ".RData")
load(file.name)
m.names <- c("MCA", "MCA2", "BOIN", "CCD", "CRM", "keyB", "mTPI","UMPBI")
sum.all <- list()
for (m.name in m.names){
   sum.all[[m.name]] <-  phase1.post.fn(lapply(1:nsimu, function(i)results[[i]][[m.name]]))
}
phase.I.pretty.tb(sum.all)

file.name <- paste0("./results/", "SimuLJ", 100*cutoff.eli, "_", nsimu, "_ncohort_", ncohort, "_fix2_",  idx, ".RData")
load(file.name)
m.names <- c("MCA2", "BOIN", "CCD", "CRM", "keyB", "mTPI","UMPBI")
sum.all <- list()
for (m.name in m.names){
   sum.all[[m.name]] <-  phase1.post.fn(lapply(1:nsimu, function(i)results[[i]][[m.name]]))
}
phase.I.pretty.tb(sum.all)

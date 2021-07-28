setwd("/home/r6user2/Documents/TQ/MCA")
#setwd("C:/Users/JINHU/Documents/ProjectCode/MCA")
source("utilities.R")


cutoff.eli <- 0.95
nsimu <- 5000
ncohort <- 10
idx <- 1
#
file.name <- paste0("./results/", "SimuMCA_NoEliLJ", 100*cutoff.eli, "_", nsimu, "_ncohort_", ncohort, "_fix3_",  idx, ".RData")
load(file.name)
m.names <- c("MCA", "BOIN", "CCD", "CRM", "keyB", "mTPI","UMPBI")
sum.all <- list()
for (m.name in m.names){
   sum.all[[m.name]] <-  phase1.post.fn(lapply(1:nsimu, function(i)results[[i]][[m.name]]))
}
phase.I.pretty.tb(sum.all)

file.name <- "./results/RealDataMCA_NoELiLJ95_5000_ncohort_13.RData"
load(file.name)
m.names <- c("MCA", "CRM")
sum.all <- list()
for (m.name in m.names){
   sum.all[[m.name]] <-  phase1.post.fn(lapply(1:nsimu, function(i)results[[i]][[m.name]]))
}
phase.I.pretty.tb(sum.all)

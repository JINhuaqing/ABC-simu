setwd("/home/r6user2/Documents/TQ/MCA")
#setwd("C:/Users/JINHU/Documents/ProjectCode/MCA")
source("utilities.R")


nsimu <- 5000
file.name1 <- "./results/SimuMCA_NoEliLJ95_5000_ncohort_12_fix1_5.RData"
file.name2 <- "./results/SimuMCA_NoEliLJ95_5000_ncohort_12_fix1_5.RData"
load(file.name1)
m.names <- c("MCAnew", "BOIN", "CCD", "CRM", "keyB", "mTPI","UMPBI")
sum.all <- list()
for (m.name in m.names){
   sum.all[[m.name]] <-  phase1.post.fn(lapply(1:nsimu, function(i)results[[i]][[m.name]]))
}
phase.I.pretty.tb(sum.all)

load(file.name2)
m.names <- c("MCAnew", "BOIN", "CCD", "CRM", "keyB", "mTPI","UMPBI")
sum.all <- list()
for (m.name in m.names){
   sum.all[[m.name]] <-  phase1.post.fn(lapply(1:nsimu, function(i)results[[i]][[m.name]]))
}
phase.I.pretty.tb(sum.all)


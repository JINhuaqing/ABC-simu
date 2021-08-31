setwd("/home/r6user2/Documents/TQ/MCA")
#setwd("C:/Users/JINHU/Documents/ProjectCode/MCA")
source("utilities.R")


nsimu <- 5000
file.name1 <- "./results/SimuMCA_NoEliLJ95_5000_ncohort_10_fix3_1.RData"
file.name2 <- "./results/RealDataMCA_NoELiLJ95_5000_ncohort_13.RData"
load(file.name1)
m.names <- c("MCA", "BOIN", "CCD", "CRM", "keyB", "mTPI","UMPBI")
sum.all <- list()
for (m.name in m.names){
   sum.all[[m.name]] <-  phase1.post.fn(lapply(1:nsimu, function(i)results[[i]][[m.name]]))
}
phase.I.pretty.tb(sum.all)

load(file.name2)
m.names <- c("MCAnew", "CRM")
sum.all <- list()
for (m.name in m.names){
   sum.all[[m.name]] <-  phase1.post.fn(lapply(1:nsimu, function(i)results[[i]][[m.name]]))
}
phase.I.pretty.tb(sum.all)


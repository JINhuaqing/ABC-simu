#setwd("C:/Users/JINHU/Documents/ProjectCode/MCA")
setwd("/home/huaqingj/MyResearch/ABC-simu")
source("utilities.R")


#cutoff.eli <- 0.95
#nsimu <- 5000
#ncohort <- 16
#Delta <- 0.05
#
#file.name <- paste0("./results/", "Simu25LJ", 100*cutoff.eli, "_", nsimu, "_ncohort_", ncohort, "_random_", Delta,  ".RData")
#load(file.name)
#print(post.process.random(results))
#
#file.name <- paste0("./results/", "Simu25", 100*cutoff.eli, "_", nsimu, "_ncohort_", ncohort, "_random_", Delta,  ".RData")
#load(file.name)
#print(post.process.random(results))
#file.name <- "./results/SimuC195_5000_ncohort_30_random_0.15.RData"
#load(file.name)
#print(post.process.random(results))
#fads


file.name <- "./results/SimuMCA_ABC_NoEliLJ95_5000_ncohort_10_random_0.05_priorDelta_10_target_20.RData"
load(file.name)
print(post.process.random(results))
file.name <- "./results/SimuMCA_ABC_NoEliLJ95_5000_ncohort_10_random_0.07_priorDelta_10_target_20.RData"
load(file.name)
print(post.process.random(results))
file.name <- "./results/SimuMCA_ABC_NoEliLJ95_5000_ncohort_10_random_0.1_priorDelta_10_target_20.RData"
load(file.name)
print(post.process.random(results))
file.name <- "./results/SimuMCA_ABC_NoEliLJ95_5000_ncohort_10_random_0.15_priorDelta_10_target_20.RData"
load(file.name)
print(post.process.random(results))
#file.name <- "./results/SimuMCA_NoEliLJ95_5000_ncohort_16_random_0.15.RData"
#load(file.name)
#print(post.process.random(results))

#setwd("C:/Users/JINHU/Documents/ProjectCode/MCA")
setwd("/home/r6user2/Documents/TQ/MCA")
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

file.name <- "./results/SimuNomono95_5000_ncohort_10_random_0.15.RData"
load(file.name)
print(post.process.random(results))

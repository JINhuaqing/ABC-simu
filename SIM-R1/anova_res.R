#library(bestNormalize)
library(dplyr)

setwd("C:/Users/Dell/Documents/ProjectCode/MCA/SIM-R1")
setwd("C:/Users/JINHU/Documents/ProjectCode/MCA/SIM-R1")
res <- read.csv("./result4ANOVA_SIM_R1B.csv")
res <- res[, -5]
names(res) <- c("Targets", "nLevels", "delta", "DiffProbs", "SS", "Y")

res$delta <- as.character(res$delta)
res <- filter(res, !delta=="CRM", SS>=12)

res$nLevels <- as.factor(res$nLevels)
res$SS <- as.factor(res$SS)
res$DiffProbs <- as.factor(res$DiffProbs)
res$Targets <- as.factor(res$Targets)
res$delta <- as.factor(res$delta)
res$Y <- res$Y/100
summary(res$delta)


fit <- aov(Y~(nLevels+SS+DiffProbs+Targets+delta)*(nLevels+SS+DiffProbs+Targets+delta), 
           data=res)
summary(fit)



# fit diagnostic
errs <- fit$residuals
norm.errs <- errs/sd(errs)
ks.test(norm.errs, pnorm)
plot(norm.errs)
norm.errs %>% density %>% plot
rnorm(10000) %>% density %>% lines(col="red")


dat <- group_by(res, DiffProbs, delta)
tb <- summarise(dat, Y=mean(Y))
tb
x.vs <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25)

plot(x.vs, tb$Y[tb$DiffProbs==0.05], ylim=c(0.2, 0.8), type="l", ylab = "Prob Diff", xlab=expression(delta))
for (dp in c(0.07, 0.1, 0.15)){
    lines(x.vs, tb$Y[tb$DiffProbs==dp])
}


SST <- sum((res$Y - mean(res$Y))**2)
sum.fit <- summary(fit)
SSE.delta <- sum(sum.fit[[1]]$`Sum Sq`[c(5)])
SSE.delta/SST * 100

4*4*3*17

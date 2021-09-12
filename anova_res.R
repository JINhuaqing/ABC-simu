library(dplyr)
setwd("C:/Users/JINHU/Documents/ProjectCode/MCA")
source("./utilities.R")
source("./anova_settings.R")
anova.ress <- dir("./results/anovaNoEliABC", pattern="*.RData", full.names = TRUE)

sep.paths <- strsplit(anova.ress, "_")
raw.labs <- sapply(sep.paths, function(i)i[3])
tmp.fn <- function(i){
    if (length(i)==5){
        return(as.numeric(i))
    }else{
       num.i <- as.numeric(i) 
       rv <- c(num.i[1:6])
       #rv <- c(num.i[1:2], 10*num.i[3]+num.i[4], num.i[5:6])
       return(rv)
    }
}
labs <- lapply(strsplit(raw.labs, ""), tmp.fn)
mat.labs <- do.call(rbind, labs)

anova.data.list <- list()
for (i in 1:length(anova.ress)){
    print(i)
    anova.fil <- anova.ress[i]
    s.labs <- mat.labs[i, ]
    s.true.labs <- rep(0, length(s.labs))
    s.true.labs[1] <- nlevels[s.labs[1]]
    s.true.labs[2] <- deltas[s.labs[2]]
    s.true.labs[3] <- sample.sizes[s.labs[3]]
    s.true.labs[4] <- diff.probs[s.labs[4]]
    s.true.labs[5] <- targets[s.labs[5]]
    s.true.labs[6] <- hs[s.labs[6]]
    
    load(anova.fil)
    res <- post.process.random(results)
    com.res <- c(s.true.labs,  unlist(res))
    
    anova.data.list[[i]] <- com.res
}


anova.data <- do.call(rbind, anova.data.list)
colnames(anova.data) <- c(c("nLevels", "delta", "SampleSizes", 
                         "DiffProbs", "Targets", "hs"), 
                       colnames(anova.data)[7:13])
anova.data.df <- as.data.frame(anova.data)
rownames(anova.data.df) <- NULL
head(anova.data.df)

anova.data.df$nLevels <- as.factor(anova.data.df$nLevels)
anova.data.df$DiffProbs <- as.factor(anova.data.df$DiffProbs)
anova.data.df$Targets <- as.factor(anova.data.df$Targets)
anova.data.df$SampleSizes<- as.factor(anova.data.df$SampleSizes)
anova.data.df$delta<- as.factor(anova.data.df$delta)
anova.data.df$hs<- as.factor(anova.data.df$hs)

# remove the cases when sample size is 12, too many outliers
anova.data.df <- anova.data.df[anova.data.df$SampleSizes!=12, ]


fit <- aov(MTD.Sel~(nLevels+delta+SampleSizes+DiffProbs+Targets+hs)*(nLevels+delta+SampleSizes+DiffProbs+Targets+hs),
           data=anova.data.df)
summary(fit)


## SST 
SST <- sum((anova.data.df$MTD.Sel - mean(anova.data.df$MTD.Sel))**2);SST
sum.fit <- summary(fit)
SSE.delta <- sum(sum.fit[[1]]$`Sum Sq`[c(2)])
SSE.delta/SST * 100
# main factors, DiffProbs, nLevels, SampleSize

# fit diagnostic
errs <- fit$residuals
norm.errs <- errs/sd(errs)
ks.test(norm.errs, pnorm)
plot(norm.errs)
norm.errs %>% density %>% plot
rnorm(10000) %>% density %>% lines(col="red")


png(filename="plots/anova_noEli.png", unit="in", height=6, width=12, res=300)
par(mfrow=c(1, 2))

# 
dat <- group_by(anova.data.df, DiffProbs, delta)
tb <- summarise(dat, Y=mean(MTD.Sel))
x.vs <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25)

plot(x.vs, tb$Y[tb$DiffProbs==0.05], ylim=c(0.2, 0.8), xaxt="n",
     type="b", ylab = "MTD selection (%)", xlab=expression(delta), col=1, lty=1, pch=1)
flag <- 2
for (dp in c(0.07, 0.1, 0.15)){
    lines(x.vs, tb$Y[tb$DiffProbs==dp], col=flag, lty=flag, pch=flag, type="b")
    flag <- flag + 1
}
legend("topright", c(expression(Delta ~'= 0.05'~"  "), 
                     expression(Delta ~'= 0.07'~"  "),
                     expression(Delta ~'= 0.10'~"  "),
                     expression(Delta ~'= 0.15'~"  ")
                     ), col=1:4, lty=1:4,, pch=1:4, ncol=2)
axis(side=1, at=x.vs, labels=c(0, 0.05, 0.10, 0.15, 0.20, "random"))

dat <- group_by(anova.data.df, DiffProbs, hs)
tb <- summarise(dat, Y=mean(MTD.Sel))
x.vs <- rev(hs)
plot(x.vs, tb$Y[tb$DiffProbs==0.05], ylim=c(0.2, 0.8), #xaxt="n",
     type="b", ylab = "MTD selection (%)", xlab="h", col=1, lty=1, pch=1, log="x")
flag <- 2
for (dp in c(0.07, 0.1, 0.15)){
    lines(x.vs, tb$Y[tb$DiffProbs==dp], col=flag, lty=flag, pch=flag, type="b")
    flag <- flag + 1
}
legend("topright", c(expression(Delta ~'= 0.05'~"  "), 
                     expression(Delta ~'= 0.07'~"  "),
                     expression(Delta ~'= 0.10'~"  "),
                     expression(Delta ~'= 0.15'~"  ")
                     ), col=1:4, lty=1:4,, pch=1:4, ncol=2)

dev.off()

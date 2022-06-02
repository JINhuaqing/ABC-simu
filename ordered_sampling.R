setwd("C:/Users/JINHU/OneDrive - connect.hku.hk/文档/ProjectCode/ABC-simu")
library(magrittr)

# Method 1
m1.x2 <- rbeta(100000, 2, 2)
m1.x1 <- runif(100000, 0, m1.x2)
m1.x3 <- runif(100000, m1.x2, 1)

# Method 2
m2.x1 <- runif(100000)
m2.x2 <- runif(100000)
m2.x3 <- runif(100000)
m2.xs <- rbind(m2.x1, m2.x2, m2.x3)
m2.xs.sorted <- apply(m2.xs, 2, sort)
m2.x1 <- m2.xs.sorted[1, ]
m2.x2 <- m2.xs.sorted[2, ]
m2.x3 <- m2.xs.sorted[3, ]


# Method 3
m3.x1 <- c()
m3.x2 <- c()
m3.x3 <- c()
while (length(m3.x1)<=100000){
    xs <- runif(3)
    if (xs[1] < xs[2] & xs[2] < xs[3]){
        m3.x1 <- c(m3.x1, xs[1])
        m3.x2 <- c(m3.x2, xs[2])
        m3.x3 <- c(m3.x3, xs[3])
    }
}

# Method 4
m4.x2 <- rbeta(100000, 1, 1)
m4.x1 <- runif(100000, 0, m4.x2)
m4.x3 <- runif(100000, m4.x2, 1)


dats1 <- list(m1.x1, m2.x1, m3.x1, m4.x1)
dats2 <- list(m1.x2, m2.x2, m3.x2, m4.x2)
dats3 <- list(m1.x3, m2.x3, m3.x3, m4.x3)
legend.fn <- function(dats){
    x1stats <- round(do.call(rbind, lapply(dats, quantile, probs=c(0.25, 0.5, 0.75))), 3);x1stats
    labs <- apply(x1stats, 1, function(curStats)sprintf("(Median:%.3f, IQR:[%.3f, %.3f])", curStats[2], curStats[1], curStats[3]))
    paste0("Method ", 1:4, " ", labs)
}
legend.fn(dats2)



png(filename="plots/toy.png", unit="in", height=6, width=12, res=300)
par(mfrow=c(1, 3))
plot(density(m1.x1), main=expression('X'[1]), col=2, lty=2, lwd=2, ylim=c(0, 4), xlab="")
lines(density(m2.x1), col=3, lty=3, lwd=2)
lines(density(m3.x1),  col=4, lty=4, lwd=2)
lines(density(m4.x1),  col=5, lty=5, lwd=2)
legend("topright", legend=legend.fn(dats1), col=2:5, lty=2:5, lwd=rep(2, 4))

plot(density(m1.x2), main=expression('X'[2]), col=2, lty=2, lwd=2, ylim=c(0, 2), xlab="")
lines(density(m2.x2), col=3, lty=3, lwd=2)
lines(density(m3.x2),  col=4, lty=4, lwd=2)
lines(density(m4.x2),  col=5, lty=5, lwd=2)
legend("topright", legend=legend.fn(dats2), col=2:5, lty=2:5, lwd=rep(2, 4))

plot(density(m1.x3), main=expression('X'[3]), col=2, lty=2, lwd=2, ylim=c(0, 4), xlab="")
lines(density(m2.x3), col=3, lty=3, lwd=2)
lines(density(m3.x3),  col=4, lty=4, lwd=2)
lines(density(m4.x3),  col=5, lty=5, lwd=2)
legend("topleft", legend=legend.fn(dats3), col=2:5, lty=2:5, lwd=rep(2, 4))
par(mfrow=c(1, 1))
dev.off()
rm(list=ls())
source("ABC_utils.R")


load("./pssprior-ndose-5-phi-20-J-20000-delta-10.RData")


tys <- c(0, 1, 3, 0, 0)
tns <- c(3, 6, 12, 0, 0)
ws <- kpws.fn(pss.prior, tys, tns, h=0.1)
nws <- ws/sum(ws)

idx <- 3
vs <- pss.prior[, idx]
kde <- density(vs, weights=nws)
jpeg("./plots/wSpSKDE.jpg", width=6, height=6, unit="in", res=500)
plot(kde, main="Density of weighted samples at dose level 2", xlab = "")
dev.off()

weighted.mean(vs, nws)
weighted.median(vs, nws)
weighted.mode(vs, nws)



weighted.mode <- function(vs, ws){
    bins <- seq(min(vs), max(vs), length.out=100)
    binWs <- c()
    for (ix in 1:(length(bins)-1)){
        binWs <- c(binWs, sum(ws[vs>=bins[ix] & vs<bins[ix+1]]))
    }
    mIdx <- which.max(binWs)
    w.mode <- (bins[mIdx]+bins[mIdx+1])/2
    return(w.mode)
}



tL <- list(x=1)

if (is.null(tL$y) | (tL$y=="T")){
    print(1)
}else if (tL$y =="F"){
    print(0)
}

T | fdsa

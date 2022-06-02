rm(list=ls())
library(ggplot2)
setwd("C:/Users/JINHU/OneDrive - connect.hku.hk/文档/ProjectCode/ABC-simu")
#setwd("/Users/jinhuaqing/Documents/Projects_Code/phaseI/")
source("utilities.R")

fils <- dir("results/", pattern="SimuMCA.+10_random.*20.RData", full.names=T);fils
for (fil in fils){
load(fil)

grp.names <- c("MTD selection", "MTD allocation", "Overdose selection", "Overdose allocation")

#grp.names <- c("MTD selection", "MTD allocation", "Overdose selection", "Overdose allocation" 
#               , "Risk of high toxicity",  "Average DLT rate")
m.names <- c("ABC", "BOIN", "CCD", "CRM", "Keyboard", "mTPI", "UMPBI")


# Bootstrap CI
numRep <- 1000
bootFil <-  paste0(strsplit(fil, ".R", T)[[1]][1], "_Boot", numRep, ".RData");bootFil
if (!file.exists(bootFil)){
    nSimu <- length(results)
    resList <- list()
    for (ix in 1:numRep){
        print(ix)
        idxs <- sample.int(nSimu, nSimu, T)
        resultsNew <- results[idxs]
        tbNew <- post.process.random(resultsNew)
        resList[[ix]] <- tbNew
    }
    nams <- names(tbNew)[1:4]
    CIsList <- list()
    CIsDiffList <- list()
    for (nam in nams){
        bt.sps <- do.call(cbind, lapply(resList, function(x)x[nam]))
        CIs <- apply(bt.sps, 1, function(x)quantile(x, c(0.025, 0.975)))
        CIs.diff <- apply(bt.sps, 1, function(x)c(-1.96*sd(x), 1.96*sd(x)))
        CIsList[[nam]] <- CIs
        CIsDiffList[[nam]] <- CIs.diff
    }

    save(CIsList, CIsDiffList, file=bootFil)
    
}
load(bootFil)

tb.names <- c("MCAnew", "BOIN", "CCD", "CRM", "keyB", "mTPI", "UMPBI")

tb <- post.process.random(results);tb
tb <- tb[, 1:4]
v.var <- c(as.vector(unlist(tb["MCAnew", ])),
           as.vector(unlist(tb["BOIN", ])), 
           as.vector(unlist(tb["CCD", ])), 
           as.vector(unlist(tb["CRM", ])), 
           as.vector(unlist(tb["keyB", ])), 
           as.vector(unlist(tb["mTPI", ])), 
           as.vector(unlist(tb["UMPBI", ]))
           ) * 100
g.var <- rep(grp.names, times=length(m.names))
m.var <- rep(m.names, each=4)

#lws <- as.vector(do.call(rbind, lapply(CIsList, function(x)x[1, ]))) * 100
#ups <- as.vector(do.call(rbind, lapply(CIsList, function(x)x[2, ]))) * 100

lws <- v.var + as.vector(do.call(rbind, lapply(CIsDiffList, function(x)x[1, ]))) * 100
ups <- v.var + as.vector(do.call(rbind, lapply(CIsDiffList, function(x)x[2, ]))) * 100


data <- data.frame(g=factor(g.var, levels=grp.names), m=factor(m.var, levels=m.names), v=v.var, vl=lws, vu=ups)
ggplot(data = data, mapping = aes(x = g, y = v, fill = m)) + 
    geom_bar(stat = 'identity', position = position_dodge(0.9)) + ylim(c(0, 75)) + 
    geom_errorbar(aes(ymax=vu, ymin=vl), position = position_dodge(0.9), width = 0.5) +
    theme(legend.position = c(0.90, 0.7) , plot.title = element_text(hjust=0.5)
          #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)
          ) +  xlab("") + ylab("Percentage (%)") + 
    guides(fill=guide_legend(title='Methods')) 
    #ggtitle("Average probability difference around the target = 0.05")
    
figPath <- paste0("plots/", strsplit(strsplit(fil, "/", T)[[1]][2], ".R", T)[[1]][1], ".jpg");figPath
ggsave(figPath, width=5, height = 4.5, units="in")

    
}

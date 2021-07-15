rm(list=ls())
library(ggplot2)
setwd("C:/Users/JINHU/Documents/ProjectCode/MCA")
#setwd("/Users/jinhuaqing/Documents/Projects_Code/phaseI/")
source("utilities.R")

fils <- dir("results/finalRes", pattern="Simu95.+random.+15");fils

fil <- paste0("results/finalRes/", fils[1])
load(fil)

grp.names <- c("MTD selection", "MTD allocation", "Overdose selection", "Overdose allocation" 
               , "Risk of high toxicity",  "Average DLT rate")
m.names <- c("MCA", "BOIN", "CCD", "CRM", "KeyBoard", "mTPI", "UMBPI")

g.var <- rep(grp.names, times=length(m.names))
m.var <- rep(m.names, each=6)

tb <- post.process.random(results);tb
tb <- tb[-1, ]
v.var <- c(as.vector(unlist(tb["MCA2", ])),
           as.vector(unlist(tb["BOIN", ])), 
           as.vector(unlist(tb["CCD", ])), 
           as.vector(unlist(tb["CRM", ])), 
           as.vector(unlist(tb["keyB", ])), 
           as.vector(unlist(tb["mTPI", ])), 
           as.vector(unlist(tb["UMPBI", ]))
           ) * 100

data <- data.frame(g=factor(g.var, levels=grp.names), m=factor(m.var, levels=m.names), v=v.var)
ggplot(data = data, mapping = aes(x = g, y = v, fill = m)) + geom_bar(stat = 'identity', position = 'dodge') +
    theme(legend.position = "bottom", plot.title = element_text(hjust=0.5)) +  xlab("") + ylab("Percentage (%)") + 
    guides(fill=guide_legend(title='Methods')) 
    #ggtitle("Average probability difference around the target = 0.05")
    
ggsave("plots/random_15.jpg", width=10, height = 4.5, units="in")

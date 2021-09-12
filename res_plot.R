rm(list=ls())
library(ggplot2)
setwd("C:/Users/JINHU/Documents/ProjectCode/MCA")
#setwd("/Users/jinhuaqing/Documents/Projects_Code/phaseI/")
source("utilities.R")

fils <- dir("results/", pattern="SimuMCA.+10_random.+07_prior");fils

fil <- paste0("results/", fils[1]);fil
load(fil)

grp.names <- c("MTD selection", "MTD allocation", "Overdose selection", "Overdose allocation")

#grp.names <- c("MTD selection", "MTD allocation", "Overdose selection", "Overdose allocation" 
#               , "Risk of high toxicity",  "Average DLT rate")
m.names <- c("ABC", "BOIN", "CCD", "CRM", "Keyboard", "mTPI", "UMPBI")

g.var <- rep(grp.names, times=length(m.names))
m.var <- rep(m.names, each=4)

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

data <- data.frame(g=factor(g.var, levels=grp.names), m=factor(m.var, levels=m.names), v=v.var)
ggplot(data = data, mapping = aes(x = g, y = v, fill = m)) + geom_bar(stat = 'identity', position = 'dodge') + ylim(c(0, 75)) + 
    theme(legend.position = c(0.90, 0.7) , plot.title = element_text(hjust=0.5)
          #axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)
          ) +  xlab("") + ylab("Percentage (%)") + 
    guides(fill=guide_legend(title='Methods')) 
    #ggtitle("Average probability difference around the target = 0.05")
    
ggsave("plots/ABC_random_c16_07.jpg", width=5, height = 4.5, units="in")

library(reshape)
library(ggplot2)

rm(list = ls())
setwd("~/Dropbox/GitHub/wgs2/")

log <- read.table("./Manu/R/log.txt", header = T)
colSums(log)

log <- melt(data = log, id.vars = "CHR")
log$variable <- factor(log$variable, levels = c("ALL", "PASS", "MAF"))

pdf("./Manu/pdfs/variants_log.pdf", width = 9, height = 5)

ggplot(log, aes(factor(CHR), value, fill = variable)) + geom_bar(stat="identity", position = "dodge") + 
  scale_fill_brewer(palette = "Set1") + theme_bw() + xlab("") + ylab("Count")

dev.off()

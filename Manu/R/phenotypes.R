library(scales)
library(ggplot2)
library(dplyr)

rm(list = ls())
setwd("~/Dropbox/GitHub/wgs2")
load("./Manu/R/mdata.rdt")

mdata$APOE = factor(mdata$APOE, levels = c("33", "22", "23", "24", "34"))

col4 <- c("#1f78b4", "#b2df8a", "#33a02c", "#fb9a99")
col6 <- c("grey70", "firebrick1", "dodgerblue3", "gold1", "chartreuse3", "darkorchid2")

blank_theme <- theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    legend.title = element_blank(),
    plot.title=element_text(size=14, face="bold")
)

tbl <- table(mdata$AD1) # AD
pct <- round(tbl/sum(tbl)*100)
lbs <- paste(names(tbl), ":", sep = "")
lbs <- paste(paste(lbs, pct), "%", sep = "")

pdf("./Manu/pdfs/phenotypes_AD.pdf")
par(mar = c(5, 4, 4, 4), cex = 1.8, col = "grey30")
pie(tbl, labels = lbs, col = col4, border = F, init.angle = 45)
dev.off()

tbl <- table(mdata$APOE) # APOE
names(tbl) = c("e3/e3", "e2/e2", "e2/e3", "e2/e4", "e3/e4")
pct <- round(tbl/sum(tbl)*100)
lbs <- paste(names(tbl), ":", sep = "")
lbs <- paste(paste(lbs, pct), "%", sep = "")

pdf("./Manu/pdfs/phenotypes_APOE.pdf")
par(cex = 1.8, col = "grey30")
pie(tbl, label = lbs, col = col6, border = F)
dev.off()

tbl <- table(mdata$Sex)  # Sex
pct <- round(tbl/sum(tbl)*100)
lbs <- c("Male", "Female")
lbs <- paste(paste(lbs, pct), "%", sep = "")

pdf("./Manu/pdfs/phenotypes_Sex.pdf")
par(cex = 1.8, col = "grey30")
pie(tbl, label = lbs, col = c("#1f78b4", "#fb9a99"), border = F)
dev.off()

pdf("./Manu/pdfs/phenotypes_APOE_bar.pdf", width = 4, height = 2)
ggplot(mdata, aes(x = AD1, fill = APOE)) + geom_bar(position = "fill", width = 0.65) +
  theme_bw() + xlab("") + ylab("Genotype Frequency") + coord_flip() +
  scale_fill_manual(values = col6, name="APOE", breaks=c("33", "22", "23", "24", "34"), 
                    labels=c("e3/e3", "e2/e2", "e2/e3", "e2/e4", "e3/e4")) +
  theme(legend.title = element_blank(), legend.key = element_blank()) 
dev.off()

pdf("./Manu/pdfs/phenotypes_Sex_bar.pdf", width = 4, height = 2)
ggplot(mdata, aes(x = AD1, fill = as.factor(Sex))) + geom_bar(position = "fill", width = 0.65) +
  theme_bw() + xlab("") + ylab("Sex Frequency") + coord_flip() +
  scale_fill_manual(values = c("#1f78b4", "#fb9a99"), name="", breaks=c("0", "1"), labels=c("M", "F")) +
  theme(legend.title = element_blank(), legend.key = element_blank()) 
dev.off()

pdf("./Manu/pdfs/phenotypes_Age.pdf", width = 3, height = 2)
ggplot(mdata, aes(x = AD1, y = Age, fill = AD1)) + geom_boxplot(width = 0.65) +
  scale_fill_manual(values = col4) + theme_bw() + xlab("") + ylab("Age") + coord_flip() +
  theme(legend.position = "none")
dev.off()

library(rstan)
library(dplyr)

rm(list = ls())
setwd("~/gitHub/Genetics/priors")

load("./mdata.rdt") 
load("./igap_t1k.rdt") # ./igap/igap.Rmd for more details
igap = igap_t1k$igap
geno = igap_t1k$geno

geno = geno[, mdata$SRR]
data <- list(N = 570, K = 4, D = 2, x = mdata[c("Age", "Sex")], y = as.numeric(mdata$AD1))

glm <- prior0 <- stan_model("./prior0.stan")
glm <- prior1 <- stan_model("./prior1.stan")
glm <- prior2 <- stan_model("./prior2.stan")
glm <- prior3 <- stan_model("./prior3.stan")
glm <- prior4 <- stan_model("./prior4.stan")

mymle <- lapply(1:nrow(geno), function(g1) { cat(g1, "\n")
  data$g = geno[g1, ]
  igap1 = igap[rownames(geno)[g1], ]
# data$p1 = igap1$Beta; data$p2 = igap1$SE # p1
# data$p1 = igap1$Beta; data$p2 = 1/igap1$SE -1 # p2
  data$prior = igap1$Beta/igap1$SE # p3, p4
  mle1 = optimizing(glm, data = data, hessian = TRUE, algorithm = "LBFGS")
  mle1.se = tryCatch(sqrt(diag(solve(-mle1$hessian)))["p"], error=function(e) NULL)
  data.frame(effect = mle1$par["p"], effect.se = mle1.se)
}); names(mymle) = rownames(geno)

mymle = as.data.frame(do.call(rbind, mymle))
mymle$Pval <- pnorm(abs(mymle$effect), sd = mymle$effect.se, lower.tail = F) * 2

mymle0 = mymle
mymle1 = mymle
mymle2 = mymle
mymle3 = mymle
mymle4 = mymle

save(mymle0, mymle1, mymle2, mymle3, mymle4, file = "./mymle.rdt")

# igap$Pval <- pnorm(abs(igap$Beta), sd = igap$SE, lower.tail = F) * 2
# igap$Pval[igap$Pval < 1e-100] = 1e-100

plot(igap$Beta, mymle0$effect, xlim = c(-1, 2), ylim = c(-1, 2), 
     xlab = "Effect in IGAP", ylab = "Effect in Bayes-GLMM", main = "No prior")
abline(a = 0, b = 1, col = "red")
plot(igap$SE, mymle0$effect.se, xlim = c(0, 0.2), ylim = c(0, 0.6), 
     xlab = "Standard Error in IGAP", ylab = "Standard Error in Bayes-GLMM", main = "No prior")
abline(a = 0, b = 1, col = "red")
plot(-log10(igap$Pvalue), -log10(mymle0$Pval), xlim = c(0, 150), ylim = c(0, 6), 
     xlab = "-log10(Pval) in IGAP", ylab = "-log10(Pval) in Bayes-GLMM", main = "No prior")
abline(a = 0, b = 1, col = "red")

par(mfrow = c(1, 2))

plot(-log10(mymle1$Pval), -log10(igap$Pvalue),
     xlab = "ADSP - Prior by method 1", ylab = "IGAP")
abline(a = 0, b = 1, col = "red")
plot(-log10(mymle1$Pval), -log10(mymle0$Pval),
     xlab = "ADSP - Prior by method 1", ylab = "ADSP - No prior")
abline(a = 0, b = 1, col = "red")

plot(-log10(mymle2$Pval), -log10(igap$Pvalue),
     xlab = "ADSP - Prior by method 2", ylab = "IGAP")
abline(a = 0, b = 1, col = "red")
plot(-log10(mymle2$Pval), -log10(mymle0$Pval),
     xlab = "ADSP - Prior by method 2", ylab = "ADSP - No prior")
abline(a = 0, b = 1, col = "red")

plot(-log10(mymle3$Pval), -log10(igap$Pvalue),
     xlab = "ADSP - Prior by method 3", ylab = "IGAP")
abline(a = 0, b = 1, col = "red")
plot(-log10(mymle3$Pval), -log10(mymle0$Pval),
     xlab = "ADSP - Prior by method 3", ylab = "ADSP - No prior")
abline(a = 0, b = 1, col = "red")

plot(-log10(mymle4$Pval), -log10(igap$Pvalue),
     xlab = "ADSP - Prior by method 4", ylab = "IGAP")
abline(a = 0, b = 1, col = "red")
plot(-log10(mymle4$Pval), -log10(mymle0$Pval),
     xlab = "ADSP - Prior by method 4", ylab = "ADSP - No prior")
abline(a = 0, b = 1, col = "red")

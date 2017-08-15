####
# Describe: use affy packages to preprocess .CEL
# Environment: R 3.3.2
# Author: Yu Chen (SKLMG)
# Input:.CEL file
# Output:Expression matrix
# R Package requirement: affy, tcltk, arulesViz, simpleaffy
# Data: 2017-08-15
####

library(affy)# or library(oligo)
library(tcltk)
library(arulesViz)
library(simpleaffy)
dir <- tk_choose.dir(caption = "Select folder")
cel.files <- list.files(path = dir, pattern = ".+\\.cel$", ignore.case = TRUE,
    full.names = TRUE, recursive = TRUE)
basename(cel.files)
data.raw <- ReadAffy(filenames = cel.files)
sampleNames(data.raw)
sampleNames(data.raw) <- paste("CHIP", 1:length(cel.files), sep = "-")
sampleNames(data.raw)
pData(data.raw)

################查看perfect match probe 以及 mismatch probe的情况
pm.data <- pm(data.raw)
head(pm.data)
mm.data <- mm(data.raw)
head(mm.data)

##################芯片灰度图
n.cel <- length(cel.files)
par(mfrow = c(ceiling(n.cel/2), 2))
par(mar = c(0.5, 0.5, 2, 0.5))
pallette.gray <- c(rep(gray(0:10/10), times = seq(1, 41, by = 4)))
for (i in 1:n.cel) image(data.raw[, i], col = pallette.gray)

#################boxplot
par(mfrow = c(1, 1))
par(mar = c(4, 4, 3, 0.5))
par(cex = 0.7)
if (n.cel > 40) par(cex = 0.5)
cols <- rainbow(n.cel * 1.2)
boxplot(data.raw, col = cols, xlab = "Sample", ylab = "Log intensity")
#histogram
par(mar = c(4, 4, 0.5, 0.5))
hist(data.raw, lty = 1:3, col = cols)
legend("topright", legend = sampleNames(data.raw), lty = 1:3, col = cols, box.col = "transparent",
    xpd = TRUE)
#################MA-plot

###################RIN
par(mfrow = c(1, 1))
par(mar = c(4, 4, 3, 0.5))
RNAdeg <- AffyRNAdeg(data.raw)
summaryAffyRNAdeg(RNAdeg)
plotAffyRNAdeg(RNAdeg, cols = cols)
legend("topleft", legend = sampleNames(data.raw), lty = 1, col = cols, box.col = "transparent",
    xpd = TRUE)
box()#理想状况下各样品的线（分段）是平行的。



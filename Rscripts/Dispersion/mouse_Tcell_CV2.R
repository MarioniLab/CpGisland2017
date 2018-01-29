# script for calculating alpha-overdispersions for supplementary figures
library(e1071)
library(limSolve)
library(statmod)

# NB: change these to the appropriate paths
tcell.cells <- read.table("~/Dropbox/Tcell/Tcell_SFnorm.tsv",
                          sep="\t", h=T, stringsAsFactors=F)
rownames(tcell.cells) <- tcell.cells$gene_id
tcell.cells <- tcell.cells[grepl(rownames(tcell.cells), pattern="ENS"), ]

tcell.means <- rowMeans(tcell.cells[, 1:(dim(tcell.cells)[2]-1)])

# variance should be calculated on the linear normalized counts, not the log2 of the normalized counts
tcell.vars <- apply(tcell.cells[, 1:(dim(tcell.cells)[2]-1)],
                    1, FUN=function(Q) var(Q))

# get the variance mean of the counts on the linear scale
tcell.count_var <- apply(tcell.cells[, 1:(dim(tcell.cells)[2]-1)],
                         1, FUN=function(Q) var(2**Q))
tcell.count_mean <- apply(tcell.cells[, 1:(dim(tcell.cells)[2]-1)],
                          1, FUN=function(Q) mean(2**Q))

# create gene expression groups based on average expression over cells
tcell.gene.summary <- as.data.frame(cbind(tcell.means, tcell.vars, log2(tcell.count_mean), log2(tcell.count_var)))
colnames(tcell.gene.summary) <- c("Mean", "Var", "CountVar", "CountMean")
tcell.gene.summary$CV2 <- tcell.gene.summary$Var/(tcell.gene.summary$Mean ** 2)
tcell.gene.summary$GENE <- rownames(tcell.gene.summary)
tcell.gene.summary <- tcell.gene.summary[(!tcell.gene.summary$CountMean < 0), ]


# find the minimum mean prior to fitting
minMeanForFit <- unname(quantile(tcell.gene.summary$Mean[which(tcell.gene.summary$CV2 > 0.2)], 0.8))

# select genes with mean value greater than min value for fitting
useForFit <- tcell.gene.summary$Mean <= 0.1

# fit with a gamma-distributed GLM
tcell.fit <- glmgam.fit(cbind(a0 = 1, a1tilde=1/tcell.gene.summary$Mean[!useForFit]), 
                        tcell.gene.summary$CV2[!useForFit])

tcell.gene.summary$Residual.CV2[!useForFit] <- abs(tcell.gene.summary$CV2[!useForFit] - fitted.values(tcell.fit))
tcell.gene.summary$Recip.means[!useForFit] <- 1/tcell.gene.summary$Mean[!useForFit]

################################
## plot mean-CV2 relationship ##
################################
# vector of points that follow mean points in order
xg <- seq(0, max(tcell.gene.summary$Mean[tcell.gene.summary$Mean != Inf]),
          length.out=100000)

a0 <- unname(tcell.fit$coefficients["a0"])
a1 <- unname(tcell.fit$coefficients["a1tilde"])
vfit <- (a1/xg) + a0

png("~/Dropbox/Noise_genomics/Figures/ms_figures/Supplementary_mouse-Tcell_CV2-fit.png",
    height=3.75, width=8.75, res=300, units="in")
par(mfrow=c(1, 2), mar=c(4.4, 4.4, 1.1, 2.1))
smoothScatter(x=tcell.gene.summary$Mean,
              y=tcell.gene.summary$CV2,
              xlab=expression(paste("log"[2], " Mean Expression")),
              ylab=expression(paste("Coefficient of Variation"^2)))
lines(xg, vfit, lwd=3, col='orange')

plot(x=tcell.gene.summary$Mean,
     y=tcell.gene.summary$Residual.CV2,
     pch='.',
     xlab=expression(paste("log"[2], " Mean Expression")),
     ylab=expression(paste("Residual Coefficient of Variation"^2)))
dev.off()

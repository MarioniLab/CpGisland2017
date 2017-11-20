# script for calculating alpha-overdispersions for supplementary figures
library(e1071)
library(limSolve)
library(statmod)

# NB: change these to the appropriate paths
mesc.cells <- read.table("~/Dropbox/mESC/mESC_SFnorm.tsv",
                         sep="\t", h=T, stringsAsFactors=F)
rownames(mesc.cells) <- mesc.cells$gene_id
mesc.cells <- mesc.cells[grepl(rownames(mesc.cells), pattern="ENS"), ]

# select just the G1 cells
mesc.cells <- mesc.cells[, c(colnames(mesc.cells)[grepl(colnames(mesc.cells), pattern="G1")], "gene_id")]

mesc.means <- rowMeans(mesc.cells[, 1:(dim(mesc.cells)[2]-1)])

# variance should be calculated on the linear normalized counts, not the log2 of the normalized counts
mesc.vars <- apply(mesc.cells[, 1:(dim(mesc.cells)[2]-1)],
                    1, FUN=function(Q) var(Q))

# get the variance mean of the counts on the linear scale
mesc.count_var <- apply(mesc.cells[, 1:(dim(mesc.cells)[2]-1)],
                         1, FUN=function(Q) var(2**Q))
mesc.count_mean <- apply(mesc.cells[, 1:(dim(mesc.cells)[2]-1)],
                          1, FUN=function(Q) mean(2**Q))

# create gene expression groups based on average expression over cells
mesc.gene.summary <- as.data.frame(cbind(mesc.means, mesc.vars, log2(mesc.count_mean), log2(mesc.count_var)))
colnames(mesc.gene.summary) <- c("Mean", "Var", "CountVar", "CountMean")
mesc.gene.summary$CV2 <- mesc.gene.summary$Var/(mesc.gene.summary$Mean ** 2)
mesc.gene.summary$GENE <- rownames(mesc.gene.summary)
mesc.gene.summary <- mesc.gene.summary[(!mesc.gene.summary$CountMean < 0), ]


# find the minimum mean prior to fitting
minMeanForFit <- unname(quantile(mesc.gene.summary$Mean[which(mesc.gene.summary$CV2 > 0.2)], 0.8))

# select genes with mean value greater than min value for fitting
useForFit <- 1/mesc.gene.summary$Mean <= 0.05

# fit with a gamma-distributed GLM
mesc.fit <- glmgam.fit(cbind(a0 = 1, a1tilde=1/mesc.gene.summary$Mean[!useForFit]), 
                        mesc.gene.summary$CV2[!useForFit])

mesc.gene.summary$Residual.CV2[!useForFit] <- fitted.values(mesc.fit) - mesc.gene.summary$CV2[!useForFit]
mesc.gene.summary$Recip.means[!useForFit] <- 1/mesc.gene.summary$Mean[!useForFit]

################################
## plot mean-CV2 relationship ##
################################
# vector of points that follow mean points in order
xg <- seq(0, max(mesc.gene.summary$Mean[mesc.gene.summary$Mean != Inf]),
          length.out=100000)

a0 <- unname(mesc.fit$coefficients["a0"])
a1 <- unname(mesc.fit$coefficients["a1tilde"])
vfit <- (a1/xg) + a0

png("~/Dropbox/Noise_genomics/Figures/ms_figures/Supplementary_mESC_CV2-fit.png",
    height=3.75, width=8.75, res=300, units="in")
par(mfrow=c(1, 2), mar=c(4.4, 4.4, 1.1, 2.1))
smoothScatter(x=mesc.gene.summary$Mean,
              y=mesc.gene.summary$CV2,
              xlab=expression(paste("log"[2], " Mean Expression")),
              ylab=expression(paste("Coefficient of Variation"^2)))
lines(xg, vfit, lwd=3, col='orange')

plot(x=mesc.gene.summary$Mean,
       y=mesc.gene.summary$Residual.CV2,
       pch='.',
       xlab=expression(paste("log"[2], " Mean Expression")),
       ylab=expression(paste("Residual Coefficient of Variation"^2)))
dev.off()


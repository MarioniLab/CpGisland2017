# script for calculating alpha-overdispersions for supplementary figures
library(e1071)
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

##########################
# plot the residual CV^2 #
##########################
# calculate the residual CV^2
# select genes with mean value greater than min value for fitting
useForFit <- mesc.gene.summary$Mean <= 0.1

# fit with a gamma-distributed GLM
fit <- glmgam.fit(cbind(a0 = 1, a1tilde=1/mesc.gene.summary$Mean[!useForFit]), 
                  mesc.gene.summary$CV2[!useForFit])

linear.fit <- glm.fit(cbind(a0 = 1, a1tilde=1/mesc.gene.summary$Mean[!useForFit]), 
                      mesc.gene.summary$CV2[!useForFit], family=gaussian)
# assess the quality of the fit by checking how much of the variance of the log CV^2
# is explained by the fit, and how much by the sampling variance of the estimator
# as recommended in Brennecke et al.
residualFit <- var(log(fitted.values(fit)) - log(mesc.gene.summary$CV2[!useForFit]))
totalVar <- var(log(mesc.gene.summary$CV2[!useForFit]))

# explained variance
1 - (residualFit/totalVar)

mesc.gene.summary$Residual.CV2[!useForFit] <- abs(mesc.gene.summary$CV2[!useForFit] - fitted.values(fit))

png("~/Dropbox/Noise_genomics/Figures/ms_figures/Supplementary_mESC_rCV2_density.png",
    height=5.75, width=5.75, res=300, units="in")
par(mar=c(4.6, 4.6, 2.1, 1.1))
plot(density(na.omit(mesc.gene.summary$Residual.CV2)), xlab=expression(paste("Residual CV"^2)), main="", cex.axis=1.7, lwd=2,
     cex.lab=2)
dev.off()

png("~/Dropbox/Noise_genomics/Figures/ms_figures/Supplementary_mESC_dispersion-CV2Vsalpha.png",
     height=7.25, width=5.75, res=300, units="in")
par(mfrow=c(2, 1), mar=c(4.6, 4.6, 2.1, 1.1))
plot(x=mesc.gene.summary$CV2,
     y=mesc.gene.summary$Alpha.loess,
     ylab=expression(paste(alpha, " Overdispersion")),
     xlab=expression(paste("CV"^2)),
     cex=1, pch='.')

plot(x=mesc.gene.summary$Residual.CV2,
     y=mesc.gene.summary$Alpha.loess,
     ylab=expression(paste("Loess ", alpha, " Overdispersion")),
     xlab=expression(paste("Residual CV"^2)),
     cex=1, pch='.')
dev.off()


# script for calculating alpha-overdispersions for supplementary figures
library(e1071)

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

# estimate the over dispersion paramer, alpha, using support vector regression and local regression
set.seed(42)
mesc.svm <- svm(CountMean ~ CountVar, mesc.gene.summary)
mesc.gene.summary$Alpha.svr <- residuals(mesc.svm)

mesc.loess <- loess(CountMean ~ CountVar, data=mesc.gene.summary, span=0.25)
mesc.gene.summary$Alpha.loess <- residuals(mesc.loess)

#########################################################
# plot mean-variance dependence for loess and SVR fits ##
#########################################################
png("~/Dropbox/Noise_genomics/Figures/ms_figures/Supplementary_mESC_dispersion-loess.png",
    width=3.75, height=3.75, res=300, units="in")
par(mar=c(4.6, 4.6, 2.1, 1.1))
smoothScatter(y=mesc.gene.summary$CountVar,
              x=mesc.gene.summary$CountMean,
              ylab=expression(paste("mESC Expression ", log[2], " Variance")),
              xlab=expression(paste("mESC Expression ", log[2], " Mean")))
lines(loess.smooth(y=mesc.gene.summary$CountVar, x=mesc.loess$fitted), col='red', lwd=2)
dev.off()


png("~/Dropbox/Noise_genomics/Figures/ms_figures/Supplementary_mESC_dispersion-svr.png",
    width=3.75, height=3.75, res=300, units="in")
par(mar=c(4.6, 4.6, 2.1, 1.1))
smoothScatter(y=mesc.gene.summary$CountVar,
              x=mesc.gene.summary$CountMean,
              xlab=expression(paste("mESC Expression ", log[2], " Variance")),
              ylab=expression(paste("mESC Expression ", log[2], " Mean")))
lines(loess.smooth(y=mesc.gene.summary$CountVar, x=fitted(mesc.svm)), col='purple', lwd=3)
dev.off()

###################################################
# plot the loess and SVR residuals for comparison #
###################################################

# also plot the Pearson correlation between them
alpha.cor <- cor(mesc.gene.summary$Alpha.loess,
                 mesc.gene.summary$Alpha.svr,
                 method="pearson")

png("~/Dropbox/Noise_genomics/Figures/ms_figures/Supplementary_mESC_dispersion_loessVsSVR.png",
    height=3.75, width=3.75, res=300, units="in")
par(mar=c(4.6, 4.6, 2.1, 1.1))
smoothScatter(x=mesc.gene.summary$Alpha.loess,
              y=mesc.gene.summary$Alpha.svr,
              xlab=expression(paste("Loess ", alpha, " Overdispersion")),
              ylab=expression(paste("SVR ", alpha, " Overdispersion")),
              cex=1.4)
text(x=1, y=3, labels=paste("Pearson\nr =", round(alpha.cor, 2)),
    font=3)
dev.off()

################################################################
# plot the loess and SVR residuals against the mean expression #
################################################################
res.loess <- loess(Alpha.loess ~ Mean, data=mesc.gene.summary, span=0.25)
mesc.gene.summary$Alpha_r.loess <- residuals(res.loess)

res.loess <- loess(Alpha.svr ~ Mean, data=mesc.gene.summary, span=0.25)
mesc.gene.summary$Alpha_r.svr <- residuals(res.loess)

# also add a loess fit to show the lack of relationship
png("~/Dropbox/Noise_genomics/Figures/ms_figures/Supplementary_mESC_dispersion-meanVsalpha.png",
    height=2.5, width=7.75, res=300, units="in")
par(mfrow=c(1, 2), mar=c(4.6, 4.6, 2.1, 1.1))

smoothScatter(x=mesc.gene.summary$Mean,
       y=mesc.gene.summary$Alpha.loess,
       ylab=expression(paste("Loess ", alpha, " Overdispersion")),
       xlab=expression(paste("log"[2], " Mean Expression")),
       cex=1, pch='.')
lines(loess.smooth(x=mesc.gene.summary$Mean,
                   y=mesc.gene.summary$Alpha_r.loess), col='red', lwd=2)

smoothScatter(x=mesc.gene.summary$Mean,
     y=mesc.gene.summary$Alpha.svr,
     ylab=expression(paste("SVR ", alpha, " Overdispersion")),
     xlab=expression(paste("log"[2], " Mean Expression")),
     cex=1, pch='.')
lines(loess.smooth(x=mesc.gene.summary$Mean,
                   y=mesc.gene.summary$Alpha_r.svr), col='purple', lwd=2)

dev.off()
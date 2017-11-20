# script for calculating alpha-overdispersions for supplementary figures
library(e1071)

hesc.cells <- read.table("~/Dropbox/hESC/hESC_norm.tsv",
                         sep="\t", h=T, stringsAsFactors=F)

rownames(hesc.cells) <- hesc.cells$ensembl_gene_id
hesc.cells <- hesc.cells[grepl(rownames(hesc.cells), pattern="ENS"), ]

hesc.means <- rowMeans(hesc.cells[, 1:(dim(hesc.cells)[2]-1)])

# variance should be calculated on the linear normalized counts, not the log2 of the normalized counts
hesc.vars <- apply(hesc.cells[, 1:(dim(hesc.cells)[2]-1)],
                    1, FUN=function(Q) var(Q))

# get the variance mean of the counts on the linear scale
hesc.count_var <- apply(hesc.cells[, 1:(dim(hesc.cells)[2]-1)],
                         1, FUN=function(Q) var(2**Q))
hesc.count_mean <- apply(hesc.cells[, 1:(dim(hesc.cells)[2]-1)],
                          1, FUN=function(Q) mean(2**Q))

# create gene expression groups based on average expression over cells
hesc.gene.summary <- as.data.frame(cbind(hesc.means, hesc.vars, log2(hesc.count_mean), log2(hesc.count_var)))
colnames(hesc.gene.summary) <- c("Mean", "Var", "CountVar", "CountMean")
hesc.gene.summary$CV2 <- hesc.gene.summary$Var/(hesc.gene.summary$Mean ** 2)
hesc.gene.summary$GENE <- rownames(hesc.gene.summary)
hesc.gene.summary <- hesc.gene.summary[(!hesc.gene.summary$CountMean < 0), ]

# estimate the over dispersion paramer, alpha, using support vector regression and local regression
set.seed(42)
hesc.svm <- svm(CountMean ~ CountVar, hesc.gene.summary)
hesc.gene.summary$Alpha.svr <- residuals(hesc.svm)

hesc.loess <- loess(CountMean ~ CountVar, data=hesc.gene.summary, span=0.25)
hesc.gene.summary$Alpha.loess <- residuals(hesc.loess)

#########################################################
# plot mean-variance dependence for loess and SVR fits ##
#########################################################
png("~/Dropbox/Noise_genomics/Figures/ms_figures/Supplementary-hESC_dispersion-loess.png",
    width=3.75, height=3.75, res=300, units="in")
par(mar=c(4.6, 4.6, 2.1, 1.1))
smoothScatter(y=hesc.gene.summary$CountVar,
              x=hesc.gene.summary$CountMean,
              ylab=expression(paste("Cd4+ T cell Expression ", log[2], " Variance")),
              xlab=expression(paste("Cd4+ T cell Expression ", log[2], " Mean")))
lines(loess.smooth(y=hesc.gene.summary$CountVar, x=hesc.loess$fitted), col='red', lwd=2)
dev.off()


png("~/Dropbox/Noise_genomics/Figures/ms_figures/Supplementary-hESC_dispersion-svr.png",
    width=3.75, height=3.75, res=300, units="in")
par(mar=c(4.6, 4.6, 2.1, 1.1))
smoothScatter(y=hesc.gene.summary$CountVar,
              x=hesc.gene.summary$CountMean,
              xlab=expression(paste("Cd4+ T cell Expression ", log[2], " Variance")),
              ylab=expression(paste("Cd4+ T cell Expression ", log[2], " Mean")))
lines(loess.smooth(y=hesc.gene.summary$CountVar, x=fitted(hesc.svm)), col='purple', lwd=3)
dev.off()

###################################################
# plot the loess and SVR residuals for comparison #
###################################################

# also plot the Pearson correlation between them
alpha.cor <- cor(hesc.gene.summary$Alpha.loess,
                 hesc.gene.summary$Alpha.svr,
                 method="pearson")

png("~/Dropbox/Noise_genomics/Figures/ms_figures/Supplementary-hESC_dispersion_loessVsSVR.png",
    height=3.75, width=3.75, res=300, units="in")
par(mar=c(4.6, 4.6, 2.1, 1.1))
smoothScatter(x=hesc.gene.summary$Alpha.loess,
              y=hesc.gene.summary$Alpha.svr,
              xlab=expression(paste("Loess ", alpha, " Overdispersion")),
              ylab=expression(paste("SVR ", alpha, " Overdispersion")),
              cex=1.4)
text(x=1, y=3, labels=paste("Pearson\nr =", round(alpha.cor, 2)),
     font=3)
dev.off()

################################################################
# plot the loess and SVR residuals against the mean expression #
################################################################
res.loess <- loess(Alpha.loess ~ Mean, data=hesc.gene.summary, span=0.25)
hesc.gene.summary$Alpha_r.loess <- residuals(res.loess)

res.loess <- loess(Alpha.svr ~ Mean, data=hesc.gene.summary, span=0.25)
hesc.gene.summary$Alpha_r.svr <- residuals(res.loess)

# also add a loess fit to show the lack of relationship
png("~/Dropbox/Noise_genomics/Figures/ms_figures/Supplementary-hESC_dispersion-meanVsalpha.png",
    height=2.5, width=7.75, res=300, units="in")
par(mfrow=c(1, 2), mar=c(4.6, 4.6, 2.1, 1.1))

smoothScatter(x=hesc.gene.summary$Mean,
              y=hesc.gene.summary$Alpha.loess,
              ylab=expression(paste("Loess ", alpha, " Overdispersion")),
              xlab=expression(paste("log"[2], " Mean Expression")),
              cex=1, pch='.')
lines(loess.smooth(x=hesc.gene.summary$Mean,
                   y=hesc.gene.summary$Alpha_r.loess), col='red', lwd=2)

smoothScatter(x=hesc.gene.summary$Mean,
              y=hesc.gene.summary$Alpha.svr,
              ylab=expression(paste("SVR ", alpha, " Overdispersion")),
              xlab=expression(paste("log"[2], " Mean Expression")),
              cex=1, pch='.')
lines(loess.smooth(x=hesc.gene.summary$Mean,
                   y=hesc.gene.summary$Alpha_r.svr), col='purple', lwd=2)

dev.off()
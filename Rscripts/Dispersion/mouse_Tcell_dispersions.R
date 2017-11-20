# script for calculating alpha-overdispersions for supplementary figures
library(e1071)

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

# estimate the over dispersion paramer, alpha, using support vector regression and local regression
set.seed(42)
tcell.svm <- svm(CountMean ~ CountVar, tcell.gene.summary)
tcell.gene.summary$Alpha.svr <- residuals(tcell.svm)

tcell.loess <- loess(CountMean ~ CountVar, data=tcell.gene.summary, span=0.25)
tcell.gene.summary$Alpha.loess <- residuals(tcell.loess)

#########################################################
# plot mean-variance dependence for loess and SVR fits ##
#########################################################
png("~/Dropbox/Noise_genomics/Figures/ms_figures/Supplementary_mouse-Tcell_dispersion-loess.png",
    width=3.75, height=3.75, res=300, units="in")
par(mar=c(4.6, 4.6, 2.1, 1.1))
smoothScatter(y=tcell.gene.summary$CountVar,
              x=tcell.gene.summary$CountMean,
              ylab=expression(paste("Cd4+ T cell Expression ", log[2], " Variance")),
              xlab=expression(paste("Cd4+ T cell Expression ", log[2], " Mean")))
lines(loess.smooth(y=tcell.gene.summary$CountVar, x=tcell.loess$fitted), col='red', lwd=2)
dev.off()


png("~/Dropbox/Noise_genomics/Figures/ms_figures/Supplementary_mouse-Tcell_dispersion-svr.png",
    width=3.75, height=3.75, res=300, units="in")
par(mar=c(4.6, 4.6, 2.1, 1.1))
smoothScatter(y=tcell.gene.summary$CountVar,
              x=tcell.gene.summary$CountMean,
              xlab=expression(paste("Cd4+ T cell Expression ", log[2], " Variance")),
              ylab=expression(paste("Cd4+ T cell Expression ", log[2], " Mean")))
lines(loess.smooth(y=tcell.gene.summary$CountVar, x=fitted(tcell.svm)), col='purple', lwd=3)
dev.off()

###################################################
# plot the loess and SVR residuals for comparison #
###################################################

# also plot the Pearson correlation between them
alpha.cor <- cor(tcell.gene.summary$Alpha.loess,
                 tcell.gene.summary$Alpha.svr,
                 method="pearson")

png("~/Dropbox/Noise_genomics/Figures/ms_figures/Supplementary_mouse-Tcell_dispersion_loessVsSVR.png",
    height=3.75, width=3.75, res=300, units="in")
par(mar=c(4.6, 4.6, 2.1, 1.1))
smoothScatter(x=tcell.gene.summary$Alpha.loess,
              y=tcell.gene.summary$Alpha.svr,
              xlab=expression(paste("Loess ", alpha, " Overdispersion")),
              ylab=expression(paste("SVR ", alpha, " Overdispersion")),
              cex=1.4)
text(x=1, y=3, labels=paste("Pearson\nr =", round(alpha.cor, 2)),
     font=3)
dev.off()

################################################################
# plot the loess and SVR residuals against the mean expression #
################################################################
res.loess <- loess(Alpha.loess ~ Mean, data=tcell.gene.summary, span=0.25)
tcell.gene.summary$Alpha_r.loess <- residuals(res.loess)

res.loess <- loess(Alpha.svr ~ Mean, data=tcell.gene.summary, span=0.25)
tcell.gene.summary$Alpha_r.svr <- residuals(res.loess)

# also add a loess fit to show the lack of relationship
png("~/Dropbox/Noise_genomics/Figures/ms_figures/Supplementary_mouse-Tcell_dispersion-meanVsalpha.png",
    height=2.5, width=7.75, res=300, units="in")
par(mfrow=c(1, 2), mar=c(4.6, 4.6, 2.1, 1.1))

smoothScatter(x=tcell.gene.summary$Mean,
              y=tcell.gene.summary$Alpha.loess,
              ylab=expression(paste("Loess ", alpha, " Overdispersion")),
              xlab=expression(paste("log"[2], " Mean Expression")),
              cex=1, pch='.')
lines(loess.smooth(x=tcell.gene.summary$Mean,
                   y=tcell.gene.summary$Alpha_r.loess), col='red', lwd=2)

smoothScatter(x=tcell.gene.summary$Mean,
              y=tcell.gene.summary$Alpha.svr,
              ylab=expression(paste("SVR ", alpha, " Overdispersion")),
              xlab=expression(paste("log"[2], " Mean Expression")),
              cex=1, pch='.')
lines(loess.smooth(x=tcell.gene.summary$Mean,
                   y=tcell.gene.summary$Alpha_r.svr), col='purple', lwd=2)

dev.off()
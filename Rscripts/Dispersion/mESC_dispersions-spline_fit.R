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

##############################################
# smooth spline fitting between mean and CV2 #
##############################################
# calculate the residual CV^2
spline.fit <- smooth.spline(x=mesc.gene.summary$Mean,
                            y=mesc.gene.summary$CV2)

# ~ 60 values are dropped across the whole expression spectrum
keep.values <- mesc.gene.summary$Mean %in% spline.fit$x
mesc.gene.summary$Spline.Fit <- NA
mesc.gene.summary[keep.values, ]$Spline.Fit <- spline.fit$y

# fit a predicted spline at the missing points
mesc.gene.summary[!keep.values, ]$Spline.Fit <- predict(spline.fit, x=mesc.gene.summary$Mean[!keep.values])$y

smoothScatter(mesc.gene.summary$Mean,
              mesc.gene.summary$CV2,
              xlab=expression(paste("Mean Log"[2], "  Normalized Expression")),
              ylab=expression("CV"^2))
lines(spline.fit)


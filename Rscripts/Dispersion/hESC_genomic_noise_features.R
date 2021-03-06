## Investigating genomic factors that influence gene expression noise
library(ggplot2)
library(reshape2)
library(Rtsne)
library(limSolve)
library(statmod)
library(biomaRt)
library(e1071)
source("~/Dropbox/R_sessions/GGMike/palette_256.R")
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

hesc.cells <- read.table("~/Dropbox/hESC/hESC_norm.tsv",
                         sep="\t", h=T, stringsAsFactors=F)
rownames(hesc.cells) <- hesc.cells$ensembl_gene_id
hesc.cells <- hesc.cells[grepl(rownames(hesc.cells), pattern="ENS"), ]

hesc.meta <- read.table("~/Dropbox/hESC/hESC-meta.tsv",
                        h=T, stringsAsFactors=F, sep="\t")
hesc.meta$Sample <- paste0("X", hesc.meta$Sample)

hesc.means <- rowMeans(hesc.cells[, 1:(dim(hesc.cells)[2]-1)])
hesc.vars <- apply(hesc.cells[, 1:(dim(hesc.cells)[2]-1)],
                   1, var)
hesc.median <- apply(hesc.cells[, 1:(dim(hesc.cells)[2]-1)],
                     1, median)
hesc.mad <- apply(hesc.cells[, 1:(dim(hesc.cells)[2]-1)],
                     1, FUN=function(M) mad(M, constant=1))

# get the variance mean of the counts on the linear scale
hesc.count_var <- apply(hesc.cells[, 1:(dim(hesc.cells)[2]-1)],
                        1, FUN=function(Q) var(2**Q))
hesc.count_mean <- apply(hesc.cells[, 1:(dim(hesc.cells)[2]-1)],
                         1, FUN=function(Q) mean(2**Q))

# create gene expression groups based on average expression over cells
hesc.exprs.groups <- as.factor(cut_number(hesc.means, n=10))
hesc.gene.summary <- as.data.frame(cbind(hesc.means, hesc.vars, hesc.median, hesc.mad,
                                         hesc.exprs.groups, log2(hesc.count_mean), log2(hesc.count_var)))
colnames(hesc.gene.summary) <- c("Mean", "Var", "Median", "MAD", "Group", "CountMean", "CountVar")
hesc.gene.summary$CV2 <- hesc.gene.summary$Var/(hesc.gene.summary$Mean** 2)
hesc.gene.summary$CV2[is.na(hesc.gene.summary$CV2)] <- 0
hesc.gene.summary$GENE <- rownames(hesc.gene.summary)
hesc.gene.summary <- hesc.gene.summary[(!hesc.gene.summary$CountMean < 0), ]

# estimate the over dispersion paramer, alpha, using support vector regression
set.seed(42)
hesc.svm <- svm(CountVar ~ CountMean, hesc.gene.summary)
hesc.gene.summary$Alpha <- residuals(hesc.svm)

## calculate the residual overdispersion, removing any remaining mean effect
hesc.loess <- loess(Alpha ~ Mean, data=hesc.gene.summary, span=0.2)
hesc.gene.summary$Alpha_r <- residuals(hesc.loess)

# calculate the residual CV^2
# find the minimum mean prior to fitting
minMeanForFit <- unname(quantile(hesc.gene.summary$Mean[which(hesc.gene.summary$CV2 > 0.2)], 0.8))

# select genes with mean value greater than min value for fitting
useForFit <- 1/hesc.gene.summary$Mean <= 0.05

# fit with a gamma-distributed GLM
fit <- glmgam.fit(cbind(a0 = 1, a1tilde=1/hesc.gene.summary$Mean[!useForFit]), 
                  hesc.gene.summary$CV2[!useForFit])

hesc.gene.summary$Residual.CV2[!useForFit] <- abs(hesc.gene.summary$CV2[!useForFit] - fitted.values(fit))

## Investigating genomic factors that influence gene expression noise
library(ggplot2)
library(reshape2)
library(Rtsne)
library(biomaRt)
library(e1071)
library(limSolve)
library(statmod)
library(scran)
source("~/Dropbox/R_sessions/GGMike/palette_256.R")
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

# # load the black list of pseudogenes to remove
# blacklist <- read.table("~/Dropbox/ENSEMBL/mm10/pseudogene_blacklist.tsv",
#                         h=FALSE, sep="\t", stringsAsFactors=FALSE)

# mesc.cells <- read.table("~/Dropbox/mESC/mESC-Ziegenhain_SFnorm.tsv",
#                            sep="\t", h=T, stringsAsFactors=F)

mesc.cells <- read.table("~/Dropbox/mESC/mESC_SFnorm.tsv",
                         sep="\t", h=T, stringsAsFactors=F)

rownames(mesc.cells) <- mesc.cells$gene_id
mesc.cells <- mesc.cells[grepl(rownames(mesc.cells), pattern="ENS"), ]
# mesc.cells <- mesc.cells[!rownames(mesc.cells) %in% blacklist$V1, ]

# select just the G1 cells
mesc.cells <- mesc.cells[, c(colnames(mesc.cells)[grepl(colnames(mesc.cells), pattern="G1")], "gene_id")]

mesc.means <- rowMeans(mesc.cells[, 1:(dim(mesc.cells)[2]-1)])

# variance should be calculated on the linear normalized counts, not the log2 of the normalized counts
mesc.vars <- apply(mesc.cells[, 1:(dim(mesc.cells)[2]-1)],
                   1, FUN=function(Q) var(Q))
mesc.median <- apply(mesc.cells[, 1:(dim(mesc.cells)[2]-1)],
                   1, median)
mesc.mad <- apply(mesc.cells[, 1:(dim(mesc.cells)[2]-1)],
                   1, FUN=function(M) mad(M, constant=1))

linear.counts <- (2**mesc.cells[, 1:(dim(mesc.cells)[2]-1)]) - 1
linear.counts[linear.counts < 0] <- 0

# get the variance mean of the counts on the linear scale
mesc.count_var <- apply(linear.counts,
                        1, FUN=var)
mesc.count_mean <- apply(linear.counts,
                         1, FUN=mean)
mesc.count_cv2 <- mesc.count_var/(mesc.count_mean ** 2)
useGenes <- log2(mesc.count_mean) > 0.0

# what is the proportion of cells that express any given gene?  This should be linked to the inherent
# variability, i.e. low average expression in 100 cells vs moderate average expression in 10 cells
prop.cells <- apply(mesc.cells[, 1:(dim(mesc.cells)[2]-1)] >= 1,
                    1, sum)/(ncol(mesc.cells)-1)

# does it make more sense to think about what the expression of each gene is where there is non-zero expression?
nz.mean <- mesc.median <- apply(mesc.cells[, 1:(dim(mesc.cells)[2]-1)],
                                1, function(Q) mean(Q[Q > 1]))

# use the distance to median as an alternative mean-independent measure of variability
# DM also needs a length normalization, maybe they all do?
mesc.dm <- DM(mesc.count_mean[useGenes], mesc.count_cv2[useGenes], win.size=100)

# calculate correlation between SP1 expression and all other genes.
sp1.exprs <- mesc.cells["ENSMUSG00000001280", 1:(dim(mesc.cells)[2]-1)]
sp3.exprs <- mesc.cells["ENSMUSG00000027109", 1:(dim(mesc.cells)[2]-1)]
sp1.cor <- apply(mesc.cells[, 1:(dim(mesc.cells)[2]-1)],
                 1, FUN=function(Q) cor(Q, t(sp1.exprs)))

sp1.sp3.ratio <- sp1.exprs/sp3.exprs
sp.ratio.cor <- apply(mesc.cells[, 1:(dim(mesc.cells)[2]-1)],
                      1, FUN=function(Q) cor(Q, t(sp1.sp3.ratio)))

# create gene expression groups based on average expression over cells
mesc.exprs.groups <- as.factor(cut_number(mesc.means, n=5))
mesc.gene.summary <- as.data.frame(cbind(mesc.means[useGenes], mesc.vars[useGenes], mesc.median[useGenes], mesc.mad[useGenes],
                                         mesc.exprs.groups[useGenes], sp1.cor[useGenes], sp.ratio.cor[useGenes],
                                         log2(mesc.count_mean)[useGenes], log2(mesc.count_var)[useGenes], mesc.dm, prop.cells[useGenes]))
colnames(mesc.gene.summary) <- c("Mean", "Var", "Median", "MAD", "Group", "SP1.COR", "SP.RATIO", "CountMean", "CountVar", 
                                 "DM", "PropExprs")
mesc.gene.summary$CV2 <- mesc.gene.summary$Var/(mesc.gene.summary$Mean ** 2)
mesc.gene.summary$CV <- mesc.gene.summary$Var/mesc.gene.summary$Mean
mesc.gene.summary$GENE <- rownames(mesc.gene.summary)
mesc.gene.summary <- mesc.gene.summary[(!mesc.gene.summary$CountMean < 0), ]

# calculate the residual CV^2
# select genes with mean value greater than min value and CV lower than max value for fitting
useForFit <- mesc.gene.summary$Mean <= 0.1

# fit with a gamma-distributed GLM
fit <- glmgam.fit(cbind(a0 = 1, a1tilde=1/mesc.gene.summary$Mean[!useForFit]), 
                  mesc.gene.summary$CV2[!useForFit], maxit = 1000)

mesc.gene.summary$Residual.CV2[!useForFit] <- abs(mesc.gene.summary$CV2[!useForFit] - fitted.values(fit))
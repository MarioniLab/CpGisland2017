## Investigating genomic factors that influence gene expression noise
library(ggplot2)
library(reshape2)
library(Rtsne)
library(glmnet)
library(biomaRt)
library(e1071)
source("~/Dropbox/R_sessions/GGMike/palette_256.R")
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

panc.cell <- read.table("~/Dropbox/pancreas/GSE86473-norm.tsv.gz",
                        sep="\t", h=T, stringsAsFactors=F)
ensembl <- useEnsembl(biomart='ensembl', dataset='hsapiens_gene_ensembl')
gene_symbol <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                     filters='external_gene_name', mart=ensembl,
                     values=panc.cell$hgnc_symbol)

panc.merge <- merge(panc.cell, gene_symbol, by.x='hgnc_symbol',
                    by.y='external_gene_name')
panc.cells <- panc.merge[, -1]
rownames(panc.cells) <- panc.cells$ensembl_gene_id
panc.cells <- panc.cells[grepl(rownames(panc.cells), pattern="ENS"), ]

panc.meta <- read.table("~/Dropbox/pancreas/GSE86473_marker_metadata.tsv",
                        h=T, stringsAsFactors=F, sep="\t")

panc.alpha <- panc.cells[, colnames(panc.cells) %in% panc.meta$Sample[panc.meta$CellType == "alpha"]]


panc.means <- rowMeans(panc.alpha[, 1:(dim(panc.alpha)[2]-1)])
panc.vars <- apply(panc.alpha[, 1:(dim(panc.alpha)[2]-1)],
                   1, var)
panc.median <- apply(panc.alpha[, 1:(dim(panc.alpha)[2]-1)],
                     1, median)
panc.mad <- apply(panc.alpha[, 1:(dim(panc.alpha)[2]-1)],
                  1, FUN=function(M) mad(M, constant=1))

# get the variance mean of the counts on the linear scale
panc.count_var <- apply(panc.alpha[, 1:(dim(panc.alpha)[2]-1)],
                        1, FUN=function(Q) var(2**Q))
panc.count_mean <- apply(panc.alpha[, 1:(dim(panc.alpha)[2]-1)],
                         1, FUN=function(Q) mean(2**Q))

# create gene expression groups based on average expression over cells
panc.exprs.groups <- as.factor(cut_number(panc.means, n=10))
panc.gene.summary <- as.data.frame(cbind(panc.means, panc.vars, panc.median, panc.mad,
                                         panc.exprs.groups, log2(panc.count_mean), log2(panc.count_var)))
colnames(panc.gene.summary) <- c("Mean", "Var", "Median", "MAD", "Group", "CountMean", "CountVar")
panc.gene.summary$CV2 <- panc.gene.summary$Var/(panc.gene.summary$Mean** 2)
panc.gene.summary$CV2[is.na(panc.gene.summary$CV2)] <- 0
panc.gene.summary$GENE <- rownames(panc.gene.summary)
panc.gene.summary <- panc.gene.summary[(!panc.gene.summary$CountMean < 0), ]

# estimate the over dispersion paramer, alpha, using support vector regression
set.seed(42)
panc.svm <- svm(CountMean ~ CountVar, panc.gene.summary)
panc.gene.summary$Alpha <- residuals(panc.svm)

## calculate the residual overdispersion, removing any remaining mean effect
panc.loess <- loess(Alpha ~ Mean, data=panc.gene.summary, span=0.2)
panc.gene.summary$Alpha_r <- residuals(panc.loess)

## Investigating genomic factors that influence gene expression noise
library(ggplot2)
library(reshape2)
library(Rtsne)
library(glmnet)
library(biomaRt)
source("~/Dropbox/R_sessions/GGMike/palette_256.R")
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

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
mesc.median <- apply(mesc.cells[, 1:(dim(mesc.cells)[2]-1)],
                   1, median)
mesc.mad <- apply(mesc.cells[, 1:(dim(mesc.cells)[2]-1)],
                   1, FUN=function(M) mad(M, constant=1))

# calculate correlation between SP1 expression and all other genes.
sp1.exprs <- mesc.cells["ENSMUSG00000001280", 1:(dim(mesc.cells)[2]-1)]
sp3.exprs <- mesc.cells["ENSMUSG00000027109", 1:(dim(mesc.cells)[2]-1)]
sp1.cor <- apply(mesc.cells[, 1:(dim(mesc.cells)[2]-1)],
                 1, FUN=function(Q) cor(Q, t(sp1.exprs)))

sp1.sp3.ratio <- sp1.exprs/sp3.exprs
sp.ratio.cor <- apply(mesc.cells[, 1:(dim(mesc.cells)[2]-1)],
                      1, FUN=function(Q) cor(Q, t(sp1.sp3.ratio)))

# create gene expression groups based on average expression over cells
mesc.exprs.groups <- as.factor(cut_number(mesc.means, n=10))
mesc.gene.summary <- as.data.frame(cbind(mesc.means, mesc.vars, mesc.median, mesc.mad, mesc.exprs.groups, sp1.cor, sp.ratio.cor))
colnames(mesc.gene.summary) <- c("Mean", "Var", "Median", "MAD", "Group", "SP1.COR", "SP.RATIO")
mesc.gene.summary$CV2 <- mesc.gene.summary$Var/(mesc.gene.summary$Mean ** 2)
mesc.gene.summary$GENE <- rownames(mesc.gene.summary)
mesc.gene.summary <- mesc.gene.summary[(!mesc.gene.summary$Mean <= 0.05), ]

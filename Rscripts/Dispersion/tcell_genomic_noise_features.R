## Investigating genomic factors that influence gene expression noise
library(ggplot2)
library(reshape2)
library(Rtsne)
library(e1071)
library(limSolve)
library(statmod)

tcell.cells <- read.table("~/Dropbox/Tcell/Tcell_SFnorm.tsv",
                          sep="\t", h=T, stringsAsFactors=F)
rownames(tcell.cells) <- tcell.cells$gene_id
tcell.cells <- tcell.cells[grepl(rownames(tcell.cells), pattern="ENS"), ]
# tcell.cells <- tcell.cells[!rownames(tcell.cells) %in% blacklist$V1, ]

tcell.means <- rowMeans(tcell.cells[, 1:(dim(tcell.cells)[2]-1)])
tcell.vars <- apply(tcell.cells[, 1:(dim(tcell.cells)[2]-1)],
                    1, var)
tcell.median <- apply(tcell.cells[, 1:(dim(tcell.cells)[2]-1)],
                      1, median)
tcell.mad <- apply(tcell.cells[, 1:(dim(tcell.cells)[2]-1)],
                      1, mad)

# calculate correlation between SP1 expression and all other genes.
sp1.exprs <- tcell.cells["ENSMUSG00000001280", 1:(dim(tcell.cells)[2]-1)]
sp3.exprs <- tcell.cells["ENSMUSG00000027109", 1:(dim(tcell.cells)[2]-1)]
sp1.cor <- apply(tcell.cells[, 1:(dim(tcell.cells)[2]-1)],
                 1, FUN=function(Q) cor(Q, t(sp1.exprs)))

# get the variance mean of the counts on the linear scale
tcell.count_var <- apply(tcell.cells[, 1:(dim(tcell.cells)[2]-1)],
                         1, FUN=function(Q) var(2**Q))
tcell.count_mean <- apply(tcell.cells[, 1:(dim(tcell.cells)[2]-1)],
                          1, FUN=function(Q) mean(2**Q))

sp1.sp3.ratio <- sp1.exprs/sp3.exprs
sp.ratio.cor <- apply(tcell.cells[, 1:(dim(tcell.cells)[2]-1)],
                      1, FUN=function(Q) cor(Q, t(sp1.sp3.ratio)))

# create gene expression groups based on average expression over cells
tcell.exprs.groups <- as.factor(cut_number(tcell.means, n=10))
tcell.gene.summary <- as.data.frame(cbind(tcell.means, tcell.vars, tcell.median, tcell.mad, tcell.exprs.groups, sp1.cor, sp.ratio.cor,
                                          log2(tcell.count_mean), log2(tcell.count_var)))
colnames(tcell.gene.summary) <- c("Mean", "Var", "Median", "MAD", "Group", "SP1.COR", "SP.RATIO", "CountMean", "CountVar")
tcell.gene.summary$CV2 <- tcell.gene.summary$Var/(tcell.gene.summary$Mean ** 2)
tcell.gene.summary$CV <- tcell.gene.summary$Var/tcell.gene.summary$Mean
tcell.gene.summary$GENE <- rownames(tcell.gene.summary)
tcell.gene.summary <- tcell.gene.summary[(!tcell.gene.summary$CountMean < 0), ]

# calculate the residual CV^2
# find the minimum mean prior to fitting
minMeanForFit <- unname(quantile(tcell.gene.summary$Mean[which(tcell.gene.summary$CV > 2)], 0.2))

# select genes with mean value greater than min value for fitting
useForFit <- tcell.gene.summary$Mean > minMeanForFit

# fit with a gamma-distributed GLM
tcell.fit <- glmgam.fit(cbind(a0 = 1, a1tilde=1/tcell.gene.summary$Mean[useForFit]), 
                        tcell.gene.summary$CV2[useForFit])
tcell.gene.summary$Residual.CV2[useForFit] <- abs(tcell.gene.summary$CV2[useForFit] - fitted.values(tcell.fit))
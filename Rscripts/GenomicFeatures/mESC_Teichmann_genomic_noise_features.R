## Investigating genomic factors that influence gene expression noise on mESCs in serum or 2i media + LIF
library(ggplot2)
library(reshape2)
library(Rtsne)
library(glmnet)
library(biomaRt)
source("~/Dropbox/R_sessions/GGMike/palette_256.R")
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

mesc.cells <- read.table("~/Dropbox/mESC/Teichman-Serum_2i-mESC_SFnorm.tsv",
                         sep="\t", h=T, stringsAsFactors=F)
rownames(mesc.cells) <- mesc.cells$gene_id
mesc.cells <- mesc.cells[grepl(rownames(mesc.cells), pattern="ENS"), ]

# split into separate 2i, 2ia and serum1 (there's a batch effect within the serum data)
serum.cells <- colnames(mesc.cells)[grepl(colnames(mesc.cells), pattern="serum1")]
a2i.cells <- colnames(mesc.cells)[grepl(colnames(mesc.cells), pattern="a2i")]
cells.2i <- colnames(mesc.cells)[grepl(colnames(mesc.cells), pattern="_2i")]

#########################
## serum cultured cells #
#########################

serum.mesc.means <- rowMeans(mesc.cells[, serum.cells])
serum.mesc.vars <- apply(mesc.cells[, serum.cells],
                   1, var)
serum.mesc.median <- apply(mesc.cells[, serum.cells],
                     1, median)
serum.mesc.mad <- apply(mesc.cells[, serum.cells],
                  1, FUN=function(M) mad(M, constant=1))

# create gene expression groups based on average expression over cells
serum.mesc.exprs.groups <- as.factor(cut_number(serum.mesc.means, n=10))
serum.mesc.gene.summary <- as.data.frame(cbind(serum.mesc.means, serum.mesc.vars, 
                                               serum.mesc.median, serum.mesc.mad, 
                                               serum.mesc.exprs.groups))
colnames(serum.mesc.gene.summary) <- c("serum.Mean", "serum.Var", "serum.Median", "serum.MAD", "serum.Group")
serum.mesc.gene.summary$serum.CV2 <- serum.mesc.gene.summary$serum.Var/(serum.mesc.gene.summary$serum.Mean ** 2)
serum.mesc.gene.summary$GENE <- rownames(serum.mesc.gene.summary)
serum.mesc.gene.summary <- serum.mesc.gene.summary[(!serum.mesc.gene.summary$serum.Mean <= 0.05), ]

#######################
## a2i cultured cells #
#######################

a2i.mesc.means <- rowMeans(mesc.cells[, a2i.cells])
a2i.mesc.vars <- apply(mesc.cells[, a2i.cells],
                         1, var)
a2i.mesc.median <- apply(mesc.cells[, a2i.cells],
                           1, median)
a2i.mesc.mad <- apply(mesc.cells[, a2i.cells],
                        1, FUN=function(M) mad(M, constant=1))

# create gene expression groups based on average expression over cells
a2i.mesc.exprs.groups <- as.factor(cut_number(a2i.mesc.means, n=10))
a2i.mesc.gene.summary <- as.data.frame(cbind(a2i.mesc.means, a2i.mesc.vars, 
                                             a2i.mesc.median, a2i.mesc.mad, 
                                             a2i.mesc.exprs.groups))
colnames(a2i.mesc.gene.summary) <- c("a2i.Mean", "a2i.Var", "a2i.Median", "a2i.MAD", "a2i.Group")
a2i.mesc.gene.summary$a2i.CV2 <- a2i.mesc.gene.summary$a2i.Var/(a2i.mesc.gene.summary$a2i.Mean ** 2)
a2i.mesc.gene.summary$GENE <- rownames(a2i.mesc.gene.summary)
a2i.mesc.gene.summary <- a2i.mesc.gene.summary[(!a2i.mesc.gene.summary$a2i.Mean <= 0.05), ]

#######################
## 2i cultured cells #
#######################

n2i.mesc.means <- rowMeans(mesc.cells[, cells.2i])
n2i.mesc.vars <- apply(mesc.cells[, cells.2i],
                       1, var)
n2i.mesc.median <- apply(mesc.cells[, cells.2i],
                         1, median)
n2i.mesc.mad <- apply(mesc.cells[, cells.2i],
                      1, FUN=function(M) mad(M, constant=1))

# create gene expression groups based on average expression over cells
n2i.mesc.exprs.groups <- as.factor(cut_number(n2i.mesc.means, n=10))
n2i.mesc.gene.summary <- as.data.frame(cbind(n2i.mesc.means, n2i.mesc.vars, 
                                             n2i.mesc.median, n2i.mesc.mad, 
                                             n2i.mesc.exprs.groups))
colnames(n2i.mesc.gene.summary) <- c("n2i.Mean", "n2i.Var", "n2i.Median", "n2i.MAD", "n2i.Group")
n2i.mesc.gene.summary$n2i.CV2 <- n2i.mesc.gene.summary$n2i.Var/(n2i.mesc.gene.summary$n2i.Mean ** 2)
n2i.mesc.gene.summary$GENE <- rownames(n2i.mesc.gene.summary)
n2i.mesc.gene.summary <- n2i.mesc.gene.summary[(!n2i.mesc.gene.summary$n2i.Mean <= 0.05), ]

mesc.gene.summary <- Reduce(x=list("serum"=serum.mesc.gene.summary,
                                   "a2i"=a2i.mesc.gene.summary,
                                   "n2i"=n2i.mesc.gene.summary),
                            f=function(x, y) merge(x, y, by='GENE'))


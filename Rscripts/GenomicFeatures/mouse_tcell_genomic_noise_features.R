## Investigating genomic factors that influence gene expression noise over age, and between mouse strains
library(ggplot2)
library(reshape2)
library(Rtsne)
library(glmnet)
library(cowplot)
library(GGally)
source("~/Dropbox/R_sessions/GGMike/palette_256.R")
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

tcell.cells <- read.table("~/Dropbox/Tcell/Nils/tcell_SFnorm.tsv",
                         sep="\t", h=T, stringsAsFactors=F)

rownames(tcell.cells) <- tcell.cells$gene_id
tcell.cells <- tcell.cells[grepl(rownames(tcell.cells), pattern="ENS"), ]

tcell.meta <- read.table("~/Dropbox/Tcell/Nils/tcell-meta.tsv",
                         h=T, sep="\t", stringsAsFactors=F)
rownames(tcell.meta) <- tcell.meta$Sample
naive.tcells <- tcell.meta$Sample[tcell.meta$Condition == 'naive']
active.tcells <- tcell.meta$Sample[tcell.meta$Condition == 'active']

tcell.meta.active <- tcell.meta[intersect(tcell.meta$Sample, active.tcells), ]
tcell.meta <- tcell.meta[intersect(tcell.meta$Sample, naive.tcells), ]
tcells.active.cells <- tcell.cells[, intersect(colnames(tcell.cells), active.tcells)]
tcell.cells <- tcell.cells[, intersect(colnames(tcell.cells), naive.tcells)]

# split into young/old and B6/CAST for comparisons
young.b6 <- tcell.meta$Sample[tcell.meta$Age == "Young" & tcell.meta$Strain == "B6"]
old.b6 <- tcell.meta$Sample[tcell.meta$Age == "Old" & tcell.meta$Strain == "B6"]
young.cast <- tcell.meta$Sample[tcell.meta$Age == "Young" & tcell.meta$Strain == "CAST"]
old.cast <- tcell.meta$Sample[tcell.meta$Age == "Old" & tcell.meta$Strain == "CAST"]

# split into naive/active young and old
young.active <- tcell.meta.active$Sample[tcell.meta.active$Age == "Young"]
old.active <- tcell.meta.active$Sample[tcell.meta.active$Age == "Old"]


#############################
## Young activated B6 mice ##
#############################
active.young.b6 <- tcells.active.cells[, intersect(colnames(tcells.active.cells), young.active)]
active.means.young.b6 <- rowMeans(active.young.b6)
active.means.young.b6[active.means.young.b6 < 0] <- 0
active.vars.young.b6 <- apply(active.young.b6,
                             1, var)

active.median.young.b6 <- apply(active.young.b6,
                               1, median)
active.mad.young.b6 <- apply(active.young.b6,
                            1, FUN=function(M) mad(M, constant=1))


# create gene expression groups based on average expression over cells
active.young.exprs.groups <- as.factor(cut_number(active.means.young.b6, n=4))
active.young.gene.summary <- as.data.frame(cbind(active.means.young.b6, active.vars.young.b6, 
                                             active.median.young.b6, active.mad.young.b6,
                                             active.young.exprs.groups))
colnames(active.young.gene.summary) <- c("Mean", "Var", "Group")
active.young.gene.summary$CV2 <- active.young.gene.summary$Var/(active.young.gene.summary$Mean ** 2)
active.young.gene.summary$CV2[is.infinite(active.young.gene.summary$CV2)] <- 0
active.young.gene.summary$GENE <- rownames(active.young.gene.summary)

# set any gene expression values < 0 to 0
###################
## Young B6 mice ##
###################
tcell.young.b6 <- tcell.cells[, intersect(colnames(tcell.cells), intersect(young.b6, naive.tcells))]
tcell.means.young.b6 <- rowMeans(tcell.young.b6)
tcell.means.young.b6[tcell.means.young.b6 < 0] <- 0
tcell.vars.young.b6 <- apply(tcell.young.b6,
                   1, var)

tcell.median.young.b6 <- apply(tcell.young.b6,
                               1, median)
tcell.mad.young.b6 <- apply(tcell.young.b6,
                               1, FUN=function(M) mad(M, constant=1))


# create gene expression groups based on average expression over cells
b6.young.exprs.groups <- as.factor(cut_number(tcell.means.young.b6, n=4))
b6.young.gene.summary <- as.data.frame(cbind(tcell.means.young.b6, tcell.vars.young.b6, 
                                             tcell.median.young.b6, tcell.mad.young.b6,
                                             b6.young.exprs.groups))
colnames(b6.young.gene.summary) <- c("Mean", "Var", "Group")
b6.young.gene.summary$CV2 <- b6.young.gene.summary$Var/(b6.young.gene.summary$Mean ** 2)
b6.young.gene.summary$CV2[is.infinite(b6.young.gene.summary$CV2)] <- 0
b6.young.gene.summary$GENE <- rownames(b6.young.gene.summary)

#################
## Old B6 mice ##
#################
tcell.old.b6 <- tcell.cells[, intersect(colnames(tcell.cells), intersect(old.b6, naive.tcells))]
tcell.means.old.b6 <- rowMeans(tcell.old.b6)
tcell.means.old.b6[tcell.means.old.b6 < 0] <- 0
tcell.vars.old.b6 <- apply(tcell.old.b6,
                             1, var)
tcell.median.old.b6 <- apply(tcell.old.b6,
                               1, median)
tcell.mad.old.b6 <- apply(tcell.old.b6,
                            1, FUN=function(M) mad(M, constant=1))


# create gene expression groups based on average expression over cells
b6.old.exprs.groups <- as.factor(cut_number(tcell.means.old.b6, n=4))
b6.old.gene.summary <- as.data.frame(cbind(tcell.means.old.b6, tcell.vars.old.b6, 
                                           tcell.median.old.b6, tcell.mad.old.b6, 
                                           b6.old.exprs.groups))
colnames(b6.old.gene.summary) <- c("Mean", "Var", "Group")
b6.old.gene.summary$CV2 <- b6.old.gene.summary$Var/(b6.old.gene.summary$Mean ** 2)
b6.old.gene.summary$CV2[is.infinite(b6.old.gene.summary$CV2)] <- 0
b6.old.gene.summary$GENE <- rownames(b6.old.gene.summary)

###########################
## Old activated B6 mice ##
###########################
active.old.b6 <- tcells.active.cells[, intersect(colnames(tcells.active.cells), old.active)]
active.means.old.b6 <- rowMeans(active.old.b6)
active.means.old.b6[active.means.old.b6 < 0] <- 0
active.vars.old.b6 <- apply(active.old.b6,
                              1, var)

active.median.old.b6 <- apply(active.old.b6,
                                1, median)
active.mad.old.b6 <- apply(active.old.b6,
                             1, FUN=function(M) mad(M, constant=1))


# create gene expression groups based on average expression over cells
active.old.exprs.groups <- as.factor(cut_number(active.means.old.b6, n=4))
active.old.gene.summary <- as.data.frame(cbind(active.means.old.b6, active.vars.old.b6, 
                                             active.median.old.b6, active.mad.old.b6,
                                             active.old.exprs.groups))
colnames(active.old.gene.summary) <- c("Mean", "Var", "Group")
active.old.gene.summary$CV2 <- active.old.gene.summary$Var/(active.old.gene.summary$Mean ** 2)
active.old.gene.summary$CV2[is.infinite(active.old.gene.summary$CV2)] <- 0
active.old.gene.summary$GENE <- rownames(active.old.gene.summary)

#####################
## Young CAST mice ##
#####################
tcell.young.cast <- tcell.cells[, intersect(colnames(tcell.cells), intersect(young.cast, naive.tcells))]
tcell.means.young.cast <- rowMeans(tcell.young.cast)
tcell.means.young.cast[tcell.means.young.cast < 0] <- 0
tcell.vars.young.cast <- apply(tcell.young.cast,
                             1, var)
tcell.median.young.cast <- apply(tcell.young.cast,
                                 1, median)
tcell.mad.young.cast <- apply(tcell.young.cast,
                              1, FUN=function(M) mad(M, constant=1))


# create gene expression groups based on average expression over cells
cast.young.exprs.groups <- as.factor(cut_number(tcell.means.young.cast, n=4))
cast.young.gene.summary <- as.data.frame(cbind(tcell.means.young.cast, tcell.vars.young.cast, 
                                               tcell.median.young.cast, tcell.mad.young.cast,
                                               cast.young.exprs.groups))
colnames(cast.young.gene.summary) <- c("Mean", "Var", "Group")
cast.young.gene.summary$CV2 <- cast.young.gene.summary$Var/(cast.young.gene.summary$Mean ** 2)
cast.young.gene.summary$CV2[is.infinite(cast.young.gene.summary$CV2)] <- 0
cast.young.gene.summary$GENE <- rownames(cast.young.gene.summary)

###################
## Old CAST mice ##
###################
tcell.old.cast <- tcell.cells[, intersect(colnames(tcell.cells), intersect(old.cast, naive.tcells))]
tcell.means.old.cast <- rowMeans(tcell.old.cast)
tcell.means.old.cast[tcell.means.old.cast < 0] <- 0
tcell.vars.old.cast <- apply(tcell.old.cast,
                           1, var)
tcell.median.old.cast <- apply(tcell.old.cast,
                               1, median)
tcell.mad.old.cast <- apply(tcell.old.cast,
                            1, FUN=function(M) mad(M, constant=1))

# create gene expression groups based on average expression over cells
cast.old.exprs.groups <- as.factor(cut_number(tcell.means.old.cast, n=4))
cast.old.gene.summary <- as.data.frame(cbind(tcell.means.old.cast, tcell.vars.old.cast,
                                             tcell.median.old.cast, tcell.mad.old.cast,
                                             cast.old.exprs.groups))
colnames(cast.old.gene.summary) <- c("Mean", "Var", "Group")
cast.old.gene.summary$CV2 <- cast.old.gene.summary$Var/(cast.old.gene.summary$Mean ** 2)
cast.old.gene.summary$CV2[is.infinite(cast.old.gene.summary$CV2)] <- 0
cast.old.gene.summary$GENE <- rownames(cast.old.gene.summary)

summary.df <- as.data.frame(cbind(tcell.means.old.b6, tcell.means.old.cast, tcell.means.young.b6,
                                  tcell.means.young.cast, tcell.vars.old.b6, tcell.vars.old.cast,
                                  tcell.vars.young.b6, tcell.vars.young.cast,
                                  b6.old.gene.summary$CV2, cast.old.gene.summary$CV2,
                                  b6.young.gene.summary$CV2, cast.young.gene.summary$CV2))

colnames(summary.df) <- c("B6.Old.Mean", "CAST.Old.Mean", "B6.Young.Mean", "CAST.Young.Mean",
                          "B6.Old.Var", "CAST.Old.Var", "B6.Young.Var", "CAST.Young.Var",
                          "B6.Old.CV2", "CAST.Old.CV2", "B6.Young.CV2", "CAST.Young.CV2")
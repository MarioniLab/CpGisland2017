## Investigating genomic factors that influence gene expression noise
library(ggplot2)
library(reshape2)
library(Rtsne)
library(glmnet)
library(biomaRt)

dc.cells <- read.table("~/Dropbox/Dendritic_cells/GSE94820/bloodDC_QCtpm.tsv",
                       sep="\t", h=T, stringsAsFactors=F)
rownames(dc.cells) <- dc.cells$gene_id
dc.cells <- dc.cells[, 1:(dim(dc.cells)[2]-1)]
dc.cells <- dc.cells[grepl(rownames(dc.cells), pattern="ENS"), ]

dc.meta <- read.table("~/Dropbox/Dendritic_cells/GSE94820/bloodDC_meta.tsv",
                      h=T, stringsAsFactors=F, sep="\t")

# select a single community of cells dominated by a single cell type
# community 6 is a mix of many cell types, probably just poor quality cells
dc.cd141 <- dc.cells[, colnames(dc.cells) %in% dc.meta$Sample[dc.meta$CellType == "CD141" & 
                                                              dc.meta$Community == 8]]

dc.means <- rowMeans(dc.cd141[, 1:(dim(dc.cd141)[2]-1)])
dc.vars <- apply(dc.cd141[, 1:(dim(dc.cd141)[2]-1)],
                 1, var)
dc.median <- apply(dc.cd141[, 1:(dim(dc.cd141)[2]-1)],
                   1, median)
dc.mads <- apply(dc.cd141[, 1:(dim(dc.cd141)[2]-1)],
                 1, FUN=function(M) (mad(M, constant=1)))

# create gene expression groups based on average expression over cells
dc.exprs.groups <- as.factor(cut_number(dc.means, n=10))
cd141.gene.summary <- as.data.frame(cbind(dc.means, dc.vars, dc.median, dc.mads,
                                        dc.exprs.groups))
colnames(cd141.gene.summary) <- c("Mean", "Var", "Median", "MAD", "Group")
cd141.gene.summary$CV2 <- cd141.gene.summary$Var/(cd141.gene.summary$Mean** 2)
cd141.gene.summary$CV2[is.na(cd141.gene.summary$CV2)] <- 0
cd141.gene.summary$GENE <- rownames(cd141.gene.summary)
cd141.gene.summary <- cd141.gene.summary[(!cd141.gene.summary$Mean <= 0.05), ]

# get the variance mean of the counts on the linear scale
dc.count_var <- apply(dc.cd141[, 1:(dim(dc.cd141)[2]-1)],
                      1, FUN=function(Q) var(2**Q))
dc.count_mean <- apply(dc.cd141[, 1:(dim(dc.cd141)[2]-1)],
                       1, FUN=function(Q) mean(2**Q))

# remove genes with a log2(var) < 0?
cd141.gene.summary$CountVar <- log2(dc.count_var[!(dc.means <= 0.05)])
cd141.gene.summary$CountMean <- log2(dc.count_mean[!(dc.means <= 0.05)])

# estimate the over dispersion paramer, cd141, using support vector regression
set.seed(42)
dc.svm <- svm(CountVar ~ CountMean, cd141.gene.summary)
cd141.gene.summary$Alpha <- residuals(dc.svm)


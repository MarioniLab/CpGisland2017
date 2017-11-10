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
dc.axl <- dc.cells[, colnames(dc.cells) %in% dc.meta$Sample[dc.meta$CellType == "AXLSIGLEC6" & 
                                                               dc.meta$Community == 5]]

dc.means <- rowMeans(dc.axl[, 1:(dim(dc.axl)[2]-1)])
dc.vars <- apply(dc.axl[, 1:(dim(dc.axl)[2]-1)],
                   1, var)
dc.median <- apply(dc.axl[, 1:(dim(dc.axl)[2]-1)],
                     1, median)
dc.mads <- apply(dc.axl[, 1:(dim(dc.axl)[2]-1)],
                   1, FUN=function(M) (mad(M, constant=1)))

# create gene expression groups based on average expression over cells
dc.exprs.groups <- as.factor(cut_number(dc.means, n=10))
axl.gene.summary <- as.data.frame(cbind(dc.means, dc.vars, dc.median, dc.mads,
                                         dc.exprs.groups))
colnames(axl.gene.summary) <- c("Mean", "Var", "Median", "MAD", "Group")
axl.gene.summary$CV2 <- axl.gene.summary$Var/(axl.gene.summary$Mean** 2)
axl.gene.summary$CV2[is.na(axl.gene.summary$CV2)] <- 0
axl.gene.summary$GENE <- rownames(axl.gene.summary)
axl.gene.summary <- axl.gene.summary[(!axl.gene.summary$Mean <= 0.05), ]

# get the variance mean of the counts on the linear scale
dc.count_var <- apply(dc.axl[, 1:(dim(dc.axl)[2]-1)],
                        1, FUN=function(Q) var(2**Q))
dc.count_mean <- apply(dc.axl[, 1:(dim(dc.axl)[2]-1)],
                         1, FUN=function(Q) mean(2**Q))

# remove genes with a log2(var) < 0?
axl.gene.summary$CountVar <- log2(dc.count_var[!(dc.means <= 0.05)])
axl.gene.summary$CountMean <- log2(dc.count_mean[!(dc.means <= 0.05)])

# estimate the over dispersion paramer, axl, using support vector regression
set.seed(42)
dc.svm <- svm(CountVar ~ CountMean, axl.gene.summary)
axl.gene.summary$Alpha <- residuals(dc.svm)


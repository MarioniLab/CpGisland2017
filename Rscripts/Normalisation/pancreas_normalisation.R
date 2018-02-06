# pancreas normalization
library(statmod)
library(ggplot2)
library(scran)
library(limSolve)
library(gplots)
library(RColorBrewer)
library(dynamicTreeCut)
library(flashClust)
library(WGCNA)
library(Rtsne)
library(cluster)
library(stringr)
library(biomaRt)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

pancreas.df <- read.table("~/Dropbox/pancreas/GSE77980_mouse_islets_rpkm.txt",
                       h=T, sep="\t", stringsAsFactors=F)
ensembl <- useEnsembl(biomart='ensembl', dataset='mmusculus_gene_ensembl', GRCh=37)
gene_symbol <- getBM(attributes=c('ensembl_gene_id', 'entrezgene'),
                     filters='entrezgene', mart=ensembl,
                     values=pancreas.df$gene.id)
pancreas.genes <- merge(pancreas.df, gene_symbol,
                        by.x='gene.id', by.y='entrezgene')
pancreas <- pancreas.genes[!duplicated(pancreas.genes$ensembl_gene_id), ]
rownames(pancreas) <- pancreas$ensembl_gene_id

pancreas.meta <- read.table("~/Dropbox/pancreas/E-GEOD-77980.sdrf.txt",
                            h=T, sep="\t", stringsAsFactors=F)

rownames(pancreas.meta) <- unlist(lapply(strsplit(gsub(pancreas.meta$Comment..Sample_title.,
                                                       pattern="Islet cell ", replacement=""),
                                                  fixed=T, split=" "),
                                         FUN=function(C) paste(C, collapse="_")))

# sce normalisation
# remove cells and genes with all 0's
gene_sparsity <- (apply(pancreas == 0, MARGIN = 1, sum)/dim(pancreas)[2])
keep_genes <- gene_sparsity < 0.9
dim(pancreas[keep_genes, ])
pancreas.nz <- pancreas[keep_genes, ]

cell_sparsity <- apply(pancreas.nz == 0, MARGIN = 2, sum)/dim(pancreas.nz)[1]
keep_cells <- cell_sparsity < 0.9
dim(pancreas.nz[, keep_cells])
pancreas.nz <- pancreas.nz[, keep_cells]
pancreas.nz <- pancreas.nz[, -c(1, dim(pancreas.nz)[2])]

sce <- newSCESet(countData = pancreas.nz)
sce <- calculateQCMetrics(sce)

clusters <- quickCluster(sce, get.spikes=FALSE, min.size=50)
sce <- computeSumFactors(sce, sizes=c(5, 10, 15, 25), positive=T,
                         assay='counts', clusters=clusters)
summary(sizeFactors((sce)))
sce <- normalize(sce)

sf.norm <- as.data.frame(exprs(sce))
sf.norm$gene_id <- rownames(sf.norm)
write.table(sf.norm, "~/Dropbox/pancreas/pancreas_SFnorm.tsv",
            sep="\t", quote=F, row.names=F)

# calculate highly variable genes from a model fit
means <- rowMeans(exprs(sce), na.rm = T)
vars <- apply(exprs(sce), 1, var, na.rm=T)
cv2 <- vars/(means^2)

smoothScatter(means, cv2)

# find the minimum mean prior to fitting
minMeanForFit <- unname(quantile(means[which(cv2 > 0.2)], 0.8))

# select genes with mean value greater than min value for fitting
useForFit <- 1/means == 0
dim(cbind(a0 = 1, a1tilde=1/means[!useForFit]))

# fit with a gamma-distributed GLM
fit <- glmgam.fit(cbind(a0 = 1, a1tilde=1/means[!useForFit]), cv2[!useForFit])

# calculate % variance explained by the model fit
resid.var <- var(fitted.values(fit) - cv2[!useForFit])
total.var <- var(cv2[!useForFit])
resid.var/total.var

# get fitted values and mean-dispersion dependence line
a0 <- unname(fit$coefficients["a0"])
a1 <- unname(fit$coefficients["a1tilde"])

xg <- seq(0, max(means[means != Inf]), length.out=100000)
vfit <- (a1/xg) + a0
lines(xg, vfit, col='black', lwd=3)

# add confidence intervals
d.f <- ncol(sf.norm) - 1
lines(xg, vfit * qchisq(0.975, d.f)/d.f, lty=2, col='black')
lines(xg, vfit * qchisq(0.025, d.f)/d.f, lty=2, col='black')

# rank genes by the significance of their deviation from the fit
# to call HVGs

a.fit <- (a1/means) + a0
varFitRatio <- vars/(a.fit * means^2)
varOrder <- order(varFitRatio, decreasing=T)

oed <- exprs(sce)[varOrder, ]
smoothScatter(means, cv2)
lines(xg, vfit, col="black", lwd=3 )
lines(xg, vfit * qchisq(0.975, d.f)/d.f, lty=2, col="black")
lines(xg, vfit * qchisq(0.025, d.f)/d.f,lty=2,col="black")

# display the 100 most highly variable genes
points(means[varOrder[1:100]], cv2[varOrder[1:100]], col=2)

pvals <- pchisq(varFitRatio * d.f, d.f, lower.tail = F)
pvals[is.na(pvals)] <- 1.0
adj.pvals <- p.adjust(pvals, method='fdr')
HVG <- adj.pvals <= 1e-2
table(HVG)

hvg.df <- na.omit(data.frame(cbind(names(HVG)[HVG], HVG[HVG])))
write.table(rownames(hvg.df),
            file = "~/Dropbox/pancreas/pancreas-HVG.tsv",
            sep="\t", row.names=F, quote=F)

pancreas.hvg <- exprs(sce)[HVG, ]
pancreas.hvg[is.na(pancreas.hvg)] <- 0

# cluster cells on the basis of their HVGs
pancreas.cor <- cor(pancreas.hvg, method='spearman')
pancreas.cor[is.na(pancreas.cor)] <- 0
spearman.dist <- as.dist(1 - pancreas.cor)
spearman.tree <- flashClust(spearman.dist, method='average')
plot(spearman.tree)
spearman.cut <- cutreeDynamicTree(dendro=spearman.tree, minModuleSize = 30, 
                                  deepSplit=F)
spearman.cols <- labels2colors(spearman.cut)
pancreas.cluster <- data.frame(cbind(colnames(pancreas.cor), spearman.cols))
colnames(pancreas.cluster) <- c("Sample", "Cluster")

# get low dimensional embedding of genes by tSNE
pancreas.hvg <- pancreas.hvg[, !duplicated(t(pancreas.hvg))]
pancreas.tsne <- Rtsne(t(pancreas.hvg), perplexity=30)
pancreas.map <- data.frame(pancreas.tsne$Y)
colnames(pancreas.map) <- c("Dim1", "Dim2")
pancreas.map$Sample <- colnames(pancreas.hvg)

pancreas.uber <- merge(pancreas.map, pancreas.cluster,
                     by='Sample')
pancreas_t <- ggplot(pancreas.uber, aes(x=Dim1, y=Dim2,
                                        colour=Cluster)) + 
  geom_point(size=1) + theme_mike() +
  scale_colour_Publication()
pancreas_t

ggsave(pancreas_t,
       filename="~/Dropbox/pancreas/plots.dir/Mouse-HVG-tSNE.png",
       height=7.5, width=7.5, dpi=90)

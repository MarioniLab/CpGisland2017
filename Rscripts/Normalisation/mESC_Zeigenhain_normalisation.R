# HVGs a la Brennecke et al
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
source("/Users/michaelm/Desktop/Single_cell_datasets/R_sessions/theme_mike.R")

mesc <- "/Users/michaelm/Desktop/Single_cell_datasets/GSE75790_ziegenhain_complete_data.txt"
mesc.df <- read.table(mesc, sep="\t", h=T, stringsAsFactors=F)

mesc.meta <- read.table("/Users/michaelm/Desktop/Single_cell_datasets/data/mESC-design.tsv",
                        sep="\t", stringsAsFactors=F, h=T)

# remove cells and genes with all 0's
gene_sparsity <- (apply(mesc.df == 0, MARGIN = 1, sum)/dim(mesc.df)[2])
keep_genes <- gene_sparsity < 0.05
dim(mesc.df[keep_genes, ])
mesc.nz <- mesc.df[keep_genes, ]
spikes <- grepl(rownames(mesc.nz), pattern='gERCC')

cell_sparsity <- apply(mesc.nz == 0, MARGIN = 2, sum)/dim(mesc.nz)[1]
keep_cells <- cell_sparsity < 0.1
dim(mesc.nz[, keep_cells])
mesc.nz <- mesc.nz[, keep_cells]

sce <- newSCESet(countData = mesc.nz)
sce <- calculateQCMetrics(sce, feature_controls=list(Spikes=spikes))
isSpike(sce) <- 'Spikes'

clusters <- quickCluster(sce, get.spikes=TRUE, min.size=40)
sce <- computeSumFactors(sce, sizes=c(5, 10, 15, 20), positive=T,
                         assay='counts', clusters=clusters)
summary(sizeFactors((sce)))

mesc.sfs <- data.frame(sizeFactors(sce))
colnames(mesc.sfs) <- c('SizeFactor')
mesc.sfs$Sample <- rownames(mesc.sfs)

summary(sizeFactors((sce)))
sce <- normalize(sce)

sf.norm <- t(t(counts(sce))/sizeFactors(sce))

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
adj.pvals <- p.adjust(pvals, method='fdr')
HVG <- adj.pvals <= 1e-2
table(HVG)

hvg.df <- data.frame(cbind(names(HVG)[HVG], HVG[HVG]))
write.table(rownames(hvg.df),
            file = "/Users/michaelm/Desktop/Single_cell_datasets/data/mESC-HVG.tsv",
            sep="\t", row.names=F, quote=F)

mesc.hvg <- mesc.nz[HVG, ]

# cluster cells on the basis of their HVGs
mesc.cor = cor(mesc.hvg, method='spearman', use='pairwise.complete.obs')
spearman.dist <- as.dist(1 - mesc.cor)
spearman.tree <- flashClust(spearman.dist, method='average')
plot(spearman.tree)
spearman.cut <- cutreeDynamicTree(dendro=spearman.tree, minModuleSize = 30, 
                                  deepSplit=F)
spearman.cols <- labels2colors(spearman.cut)
mesc.cluster <- data.frame(cbind(colnames(mesc.cor), spearman.cols))
colnames(mesc.cluster) <- c("Sample", "Cluster")

# get low dimensional embedding of genes by tSNE
mesc.tsne <- Rtsne(t(mesc.hvg), perplexity=30)
mesc.map <- data.frame(mesc.tsne$Y)
colnames(mesc.map) <- c("Dim1", "Dim2")
mesc.map$Sample <- colnames(mesc.hvg)

# merge metadata and clustering info
mesc.cluster.meta <- merge(mesc.meta, mesc.cluster,
                               by='Sample')
mesc.uber <- merge(mesc.map, mesc.cluster.meta,
                       by='Sample')

# plot!
mesc_t <- ggplot(mesc.uber, aes(x=Dim1, y=Dim2, colour=Protocol)) + 
  geom_point(size=1) + theme_mike() +
  scale_colour_Publication()

ggsave(mesc_t, filename='/Users/michaelm/Desktop/Single_cell_datasets/plots.dir/mESC/RtSNE-HVGs.png',
       height=6, width=6, dpi=300)
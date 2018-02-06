#### human embryonic stem cell data
library(statmod)
library(scater)
library(scran)
library(ggplot2)
library(limSolve)
library(biomaRt)
library(flashClust)
library(WGCNA)
library(Rtsne)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

old.sce <- readRDS("~/Dropbox/hESC/sce_all.rds")

# each batch has a different ERCC concentration, which is screwing up the
# mean variance trend!
# just select batch2, which is larger
# only select cells in G1 phase
count.mat <- counts(old.sce)[, (pData(old.sce)$batch == 2)]

hesc.meta <- as.data.frame(cbind(colnames(count.mat),
                                 as.character(old.sce$phenotype)[(pData(old.sce)$batch == 2)]))
colnames(hesc.meta) <- c("Sample", "Condition")
rownames(hesc.meta) <- hesc.meta$Sample

gene_sparsity <- (apply(count.mat == 0, MARGIN = 1, sum)/dim(count.mat)[2])
keep_genes <- gene_sparsity < 0.99
dim(count.mat[keep_genes, ])
hesc.nz <- count.mat[keep_genes, ]

# remove a cell with very low counts, order of magnitude lower than all others
cell_sparsity <- apply(hesc.nz == 0, MARGIN = 2, sum)/dim(hesc.nz)[1]
keep_cells <- cell_sparsity < 0.80
dim(hesc.nz[, keep_cells])
hesc.nz <- hesc.nz[, keep_cells]

hesc.meta <- hesc.meta[intersect(rownames(hesc.meta), colnames(hesc.nz)), ]
naive.cells <- hesc.meta$Sample[hesc.meta$Condition == "primed"]

# select ERCC spike-ins
spikes <- grepl(rownames(hesc.nz), pattern="ERCC-")
pd <- new("AnnotatedDataFrame", data=hesc.meta[intersect(colnames(hesc.nz), naive.cells), ])
sce <- newSCESet(countData=hesc.nz[, intersect(colnames(hesc.nz), naive.cells)],
                 phenoData=pd)
sce <- calculateQCMetrics(sce, feature_controls=list(Spikes=spikes))
isSpike(sce) <- "Spikes"
plotPhenoData(sce, aes(x=Condition,
                       y=log10_counts_endogenous_features, colour=log10_total_counts))

clusters <- quickCluster(sce, min.size=25, get.spikes=TRUE)
sce <- computeSumFactors(sce, sizes=c(10, 20, 30, 40), positive=T,
                         assay='counts',
                         get.spikes=TRUE)
summary(sizeFactors((sce)))
hist(sizeFactors((sce)), breaks=100)
sce <- normalize(sce)

sf.norm <- as.data.frame(exprs(sce))
write.table(hesc.meta,
            file="~/Dropbox/hESC/hESC-meta.tsv",
            sep="\t", quote=F, row.names=F)

sf.norm$gene_id <- rownames(sf.norm)

ensembl <- useEnsembl(biomart='ensembl', dataset='hsapiens_gene_ensembl')
gene_symbol <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                     filters='external_gene_name', mart=ensembl,
                     values=rownames(sf.norm))
sf.merge <- merge(sf.norm, gene_symbol, by.x='gene_id',
                  by.y='external_gene_name')

sf.merge <- sf.merge[, -1]
write.table(sf.merge,
            file="~/Dropbox/hESC/hESC_norm.tsv",
            quote=F, row.names=F, sep="\t")

# calculate highly variable genes from a model fit
exprs.nz <- exprs(sce)
means <- rowMeans(exprs.nz, na.rm = T)
vars <- apply(exprs.nz, 1, var, na.rm=T)
cv2 <- vars/(means^2)

smoothScatter(means, cv2)

# find the minimum mean prior to fitting
minMeanForFit <- unname(quantile(means[which(cv2 > 0.2)], 0.8))

# select genes with mean value greater than min value for fitting
useForFit <- 1/means < 0

# fit with a gamma-distributed GLM
fit <- glmgam.fit(cbind(a0 = 1, a1tilde=1/means[!useForFit]), 
                  cv2[!useForFit])

# calculate % variance explained by the model fit
resid.var <- var(fitted.values(fit) - cv2[!useForFit])
total.var <- var(cv2[!useForFit])

# get fitted values and mean-dispersion dependence line
a0 <- unname(fit$coefficients["a0"])
a1 <- unname(fit$coefficients["a1tilde"])

xg <- seq(0, max(means[means != Inf]), length.out=100000)
vfit <- (a1/xg) + a0

smoothScatter(x=means, y=cv2,
              xlab="Mean log Expression", ylab=expression("CV"^2))
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

oed <- exprs.nz[varOrder, ]
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
            file = "~/Dropbox/hESC/hESC-HVG.tsv",
            sep="\t", row.names=F, quote=F)

hesc.hvg <- exprs.nz[HVG, ]
hesc.hvg[is.na(hesc.hvg)] <- 0

# cluster cells on the basis of their HVGs
hesc.cor <- cor(hesc.hvg, method='spearman')
hesc.cor[is.na(hesc.cor)] <- 0
spearman.dist <- as.dist(1 - hesc.cor)
spearman.tree <- flashClust(spearman.dist, method='average')
plot(spearman.tree)
spearman.cut <- cutreeDynamicTree(dendro=spearman.tree, minModuleSize = 30, 
                                  deepSplit=F)
spearman.cols <- labels2colors(spearman.cut)
hesc.cluster <- data.frame(cbind(colnames(hesc.cor), spearman.cols))
colnames(hesc.cluster) <- c("Sample", "Cluster")

# get low dimensional embedding of genes by tSNE
hesc.hvg <- hesc.hvg[, !duplicated(t(hesc.hvg))]
hesc.tsne <- Rtsne(t(hesc.hvg), perplexity=20)
hesc.map <- data.frame(hesc.tsne$Y)
colnames(hesc.map) <- c("Dim1", "Dim2")
hesc.map$Sample <- colnames(hesc.hvg)
hesc.merge <- merge(hesc.map, hesc.meta, by='Sample')
hesc.uber <- merge(hesc.merge, hesc.cluster,
                 by='Sample')

hesc_t <- ggplot(hesc.uber,
               aes(x=Dim1, y=Dim2, colour=Condition)) + 
  geom_point(size=1) + theme_mike() +
  scale_colour_Publication()
hesc_t

ggsave(hesc_t, 
       filename="~/Dropbox/hESC/plots.dir/hESC-HVG-tsne.png",
       height=6.5, width=6.5, dpi=300)

## old and young BL6/CAST F1 naive and activated CD4 T cells from Eling et al
# tcell normalization
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

tcells <- read.table("~/Dropbox/Tcell/Nils/QC_B6_CAST_all.txt",
                     sep="\t", h=T, stringsAsFactors=F)
stim <- unlist(lapply(strsplit(colnames(tcells), split="_", fixed=T),
                      FUN=function(X) paste0(X[2])))
plate <- unlist(lapply(strsplit(colnames(tcells), split="_", fixed=T),
                       FUN=function(X) paste0(X[1])))

tcell.meta <- as.data.frame(cbind(colnames(tcells), stim, plate))
colnames(tcell.meta) <- c("Sample", "Condition", "Plate")
tcell.meta$Age <- ""
tcell.meta[tcell.meta$Plate %in% c("SS51", "SS52", "SS19", "SS25"), ]$Age <- "Young"
tcell.meta[tcell.meta$Plate %in% c("SS63", "SS64", "SS53", "SS54"), ]$Age <- "Old"

tcell.meta$Strain <- ""
tcell.meta[tcell.meta$Plate %in% c("SS51", "SS52", "SS63", "SS64"), ]$Strain <- "B6"
tcell.meta[tcell.meta$Plate %in% c("SS19", "SS25", "SS53", "SS54"), ]$Strain <- "CAST"
rownames(tcell.meta) <- tcell.meta$Sample

write.table(tcell.meta,
            file="~/Dropbox/Tcell/Nils/tcell-meta.tsv",
            sep="\t", row.names=F, quote=F)
# sce normalisation
# remove cells and genes with all 0's
gene_sparsity <- (apply(tcells == 0, MARGIN = 1, sum)/dim(tcells)[2])
keep_genes <- gene_sparsity < 0.95
dim(tcells[keep_genes, ])
tcell.nz <- tcells[keep_genes, ]

cell_sparsity <- apply(tcell.nz == 0, MARGIN = 2, sum)/dim(tcell.nz)[1]
keep_cells <- cell_sparsity < 0.9
dim(tcell.nz[, keep_cells])
tcell.nz <- tcell.nz[, keep_cells]
tcell.nz <- tcell.nz[, -c(1, dim(tcell.nz)[2])]

tcell.spikes <- grepl(rownames(tcell.nz), pattern="ERCC")
pd <- new("AnnotatedDataFrame", tcell.meta[intersect(rownames(tcell.meta), colnames(tcell.nz)), ])
sce <- newSCESet(countData = tcell.nz, phenoData=pd)
sce <- calculateQCMetrics(sce, feature_controls=list(Spikes=tcell.spikes))

plotPhenoData(sce, aes(x=Plate, y=log10_counts_endogenous_features, colour=log10_total_counts))

# plate SS54 has an overall lower yield than the other plates
plotPhenoData(sce, aes(x=log10_total_counts,
                       y=pct_exprs_top_100_features,
                       colour=Plate))

clusters <- quickCluster(sce, get.spikes=TRUE, min.size=150)
sce <- computeSumFactors(sce, sizes=c(30, 40, 60, 75), positive=T,
                         assay='counts', clusters=clusters, get.spikes=TRUE)
summary(sizeFactors((sce)))
hist(sizeFactors(sce), 100)
sce <- normalize(sce)

sf.norm <- as.data.frame(exprs(sce))
sf.norm$gene_id <- rownames(sf.norm)
write.table(sf.norm, "~/Dropbox/Tcell/Nils/tcell_SFnorm.tsv",
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
            file = "~/Dropbox/Tcell/Nils/tcell-HVG.tsv",
            sep="\t", row.names=F, quote=F)

tcell.hvg <- exprs(sce)[HVG, ]
tcell.hvg[is.na(tcell.hvg)] <- 0

# cluster cells on the basis of their HVGs
tcell.cor <- cor(tcell.hvg, method='spearman')
tcell.cor[is.na(tcell.cor)] <- 0
spearman.dist <- as.dist(1 - tcell.cor)
spearman.tree <- flashClust(spearman.dist, method='average')
plot(spearman.tree)
spearman.cut <- cutreeDynamicTree(dendro=spearman.tree, minModuleSize = 30, 
                                  deepSplit=F)
spearman.cols <- labels2colors(spearman.cut)
tcell.cluster <- data.frame(cbind(colnames(tcell.cor), spearman.cols))
colnames(tcell.cluster) <- c("Sample", "Cluster")

# get low dimensional embedding of genes by tSNE
tcell.hvg <- tcell.hvg[, !duplicated(t(tcell.hvg))]
tcell.tsne <- Rtsne(t(tcell.hvg), perplexity=30)
tcell.map <- data.frame(tcell.tsne$Y)
colnames(tcell.map) <- c("Dim1", "Dim2")
tcell.map$Sample <- colnames(tcell.hvg)

tcell.merge <- merge(tcell.map, tcell.meta,
                     by='Sample')
tcell.uber <- merge(tcell.merge, tcell.cluster,
                       by='Sample')
tcell_t <- ggplot(tcell.uber, aes(x=Dim1, y=Dim2,
                                        colour=Condition)) + 
  geom_point(size=1) + theme_mike() +
  scale_colour_Publication()
tcell_t

ggsave(tcell_t,
       filename="~/Dropbox/Tcell/Nils/plots.dir/Mouse-HVG-tSNE.png",
       height=7.5, width=7.5, dpi=90)


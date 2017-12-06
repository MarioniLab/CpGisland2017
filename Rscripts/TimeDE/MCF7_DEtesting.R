#######################################################
## GSE78167 - ER stimulated MCF7 breast cancer cells ##
#######################################################
library(DESeq2)
library(ggplot2)
library(scales)
library(stringr)
library(RColorBrewer)
library(reshape2)
library(biomaRt)
library(ggrepel)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")
source("~/Dropbox/R_sessions/Noise/human_genomic_noise_features.R")

gse78167 <- read.table("~/CI_filesystem/mnt/nas2-data/jmlab/group_folders/morgan02/Noise_genetics/Timeseries/GSE78167_counts.tsv.gz",
                       sep="\t", h=T, stringsAsFactors=F)

gse78167.meta <- read.table("~/CI_filesystem/mnt/nas2-data/jmlab/group_folders/morgan02/Noise_genetics/Timeseries/GSE78167_meta.tsv.gz",
                            sep="\t", h=T, stringsAsFactors=F)

# get gene IDs as second component of <>|<>
# these are entrez gene IDs it looks like
genes <- unlist(lapply(strsplit(gse78167$gene_id, split="|", fixed=T),
                       FUN=function(Q) paste0(Q[2])))
rownames(gse78167) <- genes

# remove column 1:3 <- gene symbols, transcript ID and length
counts <- gse78167[, -c(1:2)]

samples <- colnames(counts)
n.samples <- length(samples)

# remove sparse genes, mean > 5 reads
nonZero <- rowMeans(counts) > 5
counts.nz <- counts[nonZero, ]

# remove .int suffix from gene id
ensembl <- useEnsembl(biomart='ensembl', dataset='hsapiens_gene_ensembl')

gene_symbol <- getBM(attributes=c('external_gene_name', 'entrezgene', 'ensembl_gene_id'),
                     filters='entrezgene', values=genes,
                     mart=ensembl)

# size factor normalisation
dds <- DESeqDataSetFromMatrix(countData=floor(counts.nz), colData=gse78167.meta,
                              design= ~ Replicate + Timepoint)

# run DESeq2
dds <- DESeq(dds)

## 0hr vs 1hr
de <- results(dds, contrast=c("Timepoint", "1hrs", "0hrs"))
res01 <- as.data.frame(lfcShrink(dds, res=de, contrast=c("Timepoint", "1hrs", "0hrs")))

res01$padj[is.na(res01$padj)] <- 1.0
res01$Sig <- as.factor(as.numeric(res01$padj <= 0.01))

res01$gene_id <- rownames(res01)
res01.merge <- merge(res01,  gene_symbol,
                     by.x='gene_id', by.y='entrezgene', all.x=TRUE)

## 1hr vs 2hr
de <- results(dds, contrast=c("Timepoint", "2hrs", "1hrs"))
res12 <- as.data.frame(lfcShrink(dds, res=de, contrast=c("Timepoint", "2hrs", "1hrs")))

res12$padj[is.na(res12$padj)] <- 1.0
res12$Sig <- as.factor(as.numeric(res12$padj <= 0.01))

res12$gene_id <- rownames(res12)
res12.merge <- merge(res12,  gene_symbol,
                     by.x='gene_id', by.y='entrezgene', all.x=TRUE)

## 2hr vs 3hr
de <- results(dds, contrast=c("Timepoint", "3hrs", "2hrs"))
res23 <- as.data.frame(lfcShrink(dds, res=de, contrast=c("Timepoint", "3hrs", "2hrs")))

res23$padj[is.na(res23$padj)] <- 1.0
res23$Sig <- as.factor(as.numeric(res23$padj <= 0.01))

res23$gene_id <- rownames(res23)
res23.merge <- merge(res23,  gene_symbol,
                     by.x='gene_id', by.y='entrezgene', all.x=TRUE)

## 3hr vs 4hr
de <- results(dds, contrast=c("Timepoint", "4hrs", "3hrs"))
res34 <- as.data.frame(lfcShrink(dds, res=de, contrast=c("Timepoint", "4hrs", "3hrs")))

res34$padj[is.na(res34$padj)] <- 1.0
res34$Sig <- as.factor(as.numeric(res34$padj <= 0.01))

res34$gene_id <- rownames(res34)
res34.merge <- merge(res34,  gene_symbol,
                     by.x='gene_id', by.y='entrezgene', all.x=TRUE)

## 4hr vs 5hr
de <- results(dds, contrast=c("Timepoint", "5hrs", "4hrs"))
res45 <- as.data.frame(lfcShrink(dds, res=de, contrast=c("Timepoint", "4hrs", "5hrs")))

res45$padj[is.na(res45$padj)] <- 1.0
res45$Sig <- as.factor(as.numeric(res45$padj <= 0.01))

res45$gene_id <- rownames(res45)
res45.merge <- merge(res45,  gene_symbol,
                     by.x='gene_id', by.y='entrezgene', all.x=TRUE)

## 5hr vs 6hr
de <- results(dds, contrast=c("Timepoint", "5hrs", "6hrs"))
res56 <- as.data.frame(lfcShrink(dds, res=de, contrast=c("Timepoint", "5hrs", "6hrs")))

res56$padj[is.na(res56$padj)] <- 1.0
res56$Sig <- as.factor(as.numeric(res56$padj <= 0.01))

res56$gene_id <- rownames(res56)
res56.merge <- merge(res56,  gene_symbol,
                     by.x='gene_id', by.y='entrezgene', all.x=TRUE)


## 6hr vs 8hr
de <- results(dds, contrast=c("Timepoint", "6hrs", "8hrs"))
res68 <- as.data.frame(lfcShrink(dds, res=de, contrast=c("Timepoint", "6hrs", "8hrs")))

res68$padj[is.na(res68$padj)] <- 1.0
res68$Sig <- as.factor(as.numeric(res68$padj <= 0.01))

res68$gene_id <- rownames(res68)
res68.merge <- merge(res68,  gene_symbol,
                     by.x='gene_id', by.y='entrezgene', all.x=TRUE)


## 8hr vs 12hr
de <- results(dds, contrast=c("Timepoint", "8hrs", "12hrs"))
res812 <- as.data.frame(lfcShrink(dds, res=de, contrast=c("Timepoint", "8hrs", "12hrs")))

res812$padj[is.na(res812$padj)] <- 1.0
res812$Sig <- as.factor(as.numeric(res812$padj <= 0.01))

res812$gene_id <- rownames(res812)
res812.merge <- merge(res812,  gene_symbol,
                      by.x='gene_id', by.y='entrezgene', all.x=TRUE)


# match up to CpG islands
# 0vs1
ER.1.up <- res01.merge$ensembl_gene_id[(res01.merge$log2FoldChange > 0) & (res01.merge$Sig == 1)]
ER.1.up.size <- human.genomic.features$CGI_SIZE.kb[human.genomic.features$GENE %in% ER.1.up]
ER.1.up_df <- cbind.data.frame(human.genomic.features$GENE[human.genomic.features$GENE %in% ER.1.up], ER.1.up.size)
colnames(ER.1.up_df) <- c("GeneID", "CGI_SIZE.kb")
ER.1.up_df$Time <- "0v1"
ER.1.up_df$Direction <- "Up"


ER.1.down <- res01.merge$ensembl_gene_id[(res01.merge$log2FoldChange < 0) & (res01.merge$Sig == 1)]
ER.1.down.size <- human.genomic.features$CGI_SIZE.kb[human.genomic.features$GENE %in% ER.1.down]
ER.1.down_df <- cbind.data.frame(human.genomic.features$GENE[human.genomic.features$GENE %in% ER.1.down], ER.1.down.size)
colnames(ER.1.down_df) <- c("GeneID", "CGI_SIZE.kb")
ER.1.down_df$Time <- "0v1"
ER.1.down_df$Direction <- "Down"


# 1vs2
ER.2.up <- res12.merge$ensembl_gene_id[(res12.merge$log2FoldChange > 0) & (res12.merge$Sig == 1)]
ER.2.up.size <- human.genomic.features$CGI_SIZE.kb[human.genomic.features$GENE %in% ER.2.up]
ER.2.up_df <- cbind.data.frame(human.genomic.features$GENE[human.genomic.features$GENE %in% ER.2.up], ER.2.up.size)
colnames(ER.2.up_df) <- c("GeneID", "CGI_SIZE.kb")
ER.2.up_df$Time <- "1v2"
ER.2.up_df$Direction <- "Up"


ER.2.down <- res12.merge$ensembl_gene_id[(res12.merge$log2FoldChange < 0) & (res12.merge$Sig == 1)]
ER.2.down.size <- human.genomic.features$CGI_SIZE.kb[human.genomic.features$GENE %in% ER.2.down]
ER.2.down_df <- cbind.data.frame(human.genomic.features$GENE[human.genomic.features$GENE %in% ER.2.down], ER.2.down.size)
colnames(ER.2.down_df) <- c("GeneID", "CGI_SIZE.kb")
ER.2.down_df$Time <- "1v2"
ER.2.down_df$Direction <- "Down"


# 2vs3
ER.3.up <- res23.merge$ensembl_gene_id[(res23.merge$log2FoldChange > 0) & (res23.merge$Sig == 1)]
ER.3.up.size <- human.genomic.features$CGI_SIZE.kb[human.genomic.features$GENE %in% ER.3.up]
ER.3.up_df <- cbind.data.frame(human.genomic.features$GENE[human.genomic.features$GENE %in% ER.3.up], ER.3.up.size)
colnames(ER.3.up_df) <- c("GeneID", "CGI_SIZE.kb")
ER.3.up_df$Time <- "2v3"
ER.3.up_df$Direction <- "Up"


ER.3.down <- res23.merge$ensembl_gene_id[(res23.merge$log2FoldChange < 0) & (res23.merge$Sig == 1)]
ER.3.down.size <- human.genomic.features$CGI_SIZE.kb[human.genomic.features$GENE %in% ER.3.down]
ER.3.down_df <- cbind.data.frame(human.genomic.features$GENE[human.genomic.features$GENE %in% ER.3.down], ER.3.down.size)
colnames(ER.3.down_df) <- c("GeneID", "CGI_SIZE.kb")
ER.3.down_df$Time <- "2v3"
ER.3.down_df$Direction <- "Down"


# 3vs4
ER.4.up <- res34.merge$ensembl_gene_id[(res34.merge$log2FoldChange > 0) & (res34.merge$Sig == 1)]
ER.4.up.size <- human.genomic.features$CGI_SIZE.kb[human.genomic.features$GENE %in% ER.4.up]
ER.4.up_df <- cbind.data.frame(human.genomic.features$GENE[human.genomic.features$GENE %in% ER.4.up], ER.4.up.size)
colnames(ER.4.up_df) <- c("GeneID", "CGI_SIZE.kb")
ER.4.up_df$Time <- "3v4"
ER.4.up_df$Direction <- "Up"


ER.4.down <- res34.merge$ensembl_gene_id[(res34.merge$log2FoldChange < 0) & (res34.merge$Sig == 1)]
ER.4.down.size <- human.genomic.features$CGI_SIZE.kb[human.genomic.features$GENE %in% ER.4.down]
ER.4.down_df <- cbind.data.frame(human.genomic.features$GENE[human.genomic.features$GENE %in% ER.4.down], ER.4.down.size)
colnames(ER.4.down_df) <- c("GeneID", "CGI_SIZE.kb")
ER.4.down_df$Time <- "3v4"
ER.4.down_df$Direction <- "Down"


# 4vs5
ER.5.up <- res45.merge$ensembl_gene_id[(res45.merge$log2FoldChange > 0) & (res45.merge$Sig == 1)]
ER.5.up.size <- human.genomic.features$CGI_SIZE.kb[human.genomic.features$GENE %in% ER.5.up]
ER.5.up_df <- cbind.data.frame(human.genomic.features$GENE[human.genomic.features$GENE %in% ER.5.up], ER.5.up.size)
colnames(ER.5.up_df) <- c("GeneID", "CGI_SIZE.kb")
ER.5.up_df$Time <- "4v5"
ER.5.up_df$Direction <- "Up"


ER.5.down <- res45.merge$ensembl_gene_id[(res45.merge$log2FoldChange < 0) & (res45.merge$Sig == 1)]
ER.5.down.size <- human.genomic.features$CGI_SIZE.kb[human.genomic.features$GENE %in% ER.5.down]
ER.5.down_df <- cbind.data.frame(human.genomic.features$GENE[human.genomic.features$GENE %in% ER.5.down], ER.5.down.size)
colnames(ER.5.down_df) <- c("GeneID", "CGI_SIZE.kb")
ER.5.down_df$Time <- "4v5"
ER.5.down_df$Direction <- "Down"


# 5vs6
ER.6.up <- res56.merge$ensembl_gene_id[(res56.merge$log2FoldChange > 0) & (res56.merge$Sig == 1)]
ER.6.up.size <- human.genomic.features$CGI_SIZE.kb[human.genomic.features$GENE %in% ER.6.up]
ER.6.up_df <- cbind.data.frame(human.genomic.features$GENE[human.genomic.features$GENE %in% ER.6.up], ER.6.up.size)
colnames(ER.6.up_df) <- c("GeneID", "CGI_SIZE.kb")
ER.6.up_df$Time <- "5v6"
ER.6.up_df$Direction <- "Up"


ER.6.down <- res56.merge$ensembl_gene_id[(res56.merge$log2FoldChange < 0) & (res56.merge$Sig == 1)]
ER.6.down.size <- human.genomic.features$CGI_SIZE.kb[human.genomic.features$GENE %in% ER.6.down]
ER.6.down_df <- cbind.data.frame(human.genomic.features$GENE[human.genomic.features$GENE %in% ER.6.down], ER.6.down.size)
colnames(ER.6.down_df) <- c("GeneID", "CGI_SIZE.kb")
ER.6.down_df$Time <- "5v6"
ER.6.down_df$Direction <- "Down"


# 6vs8
ER.8.up <- res68.merge$ensembl_gene_id[(res68.merge$log2FoldChange > 0) & (res68.merge$Sig == 1)]
ER.8.up.size <- human.genomic.features$CGI_SIZE.kb[human.genomic.features$GENE %in% ER.8.up]
ER.8.up_df <- cbind.data.frame(human.genomic.features$GENE[human.genomic.features$GENE %in% ER.8.up], ER.8.up.size)
colnames(ER.8.up_df) <- c("GeneID", "CGI_SIZE.kb")
ER.8.up_df$Time <- "6v8"
ER.8.up_df$Direction <- "Up"


ER.8.down <- res68.merge$ensembl_gene_id[(res68.merge$log2FoldChange < 0) & (res68.merge$Sig == 1)]
ER.8.down.size <- human.genomic.features$CGI_SIZE.kb[human.genomic.features$GENE %in% ER.8.down]
ER.8.down_df <- cbind.data.frame(human.genomic.features$GENE[human.genomic.features$GENE %in% ER.8.down], ER.8.down.size)
colnames(ER.8.down_df) <- c("GeneID", "CGI_SIZE.kb")
ER.8.down_df$Time <- "6v8"
ER.8.down_df$Direction <- "Down"


# 8vs12
ER.12.up <- res812.merge$ensembl_gene_id[(res812.merge$log2FoldChange > 0) & (res812.merge$Sig == 1)]
ER.12.up.size <- human.genomic.features$CGI_SIZE.kb[human.genomic.features$GENE %in% ER.12.up]
ER.12.up_df <- cbind.data.frame(human.genomic.features$GENE[human.genomic.features$GENE %in% ER.12.up], ER.12.up.size)
colnames(ER.12.up_df) <- c("GeneID", "CGI_SIZE.kb")
ER.12.up_df$Time <- "8v12"
ER.12.up_df$Direction <- "Up"


ER.12.down <- res812.merge$ensembl_gene_id[(res812.merge$log2FoldChange < 0) & (res812.merge$Sig == 1)]
ER.12.down.size <- human.genomic.features$CGI_SIZE.kb[human.genomic.features$GENE %in% ER.12.down]
ER.12.down_df <- cbind.data.frame(human.genomic.features$GENE[human.genomic.features$GENE %in% ER.12.down], ER.12.down.size)
colnames(ER.12.down_df) <- c("GeneID", "CGI_SIZE.kb")
ER.12.down_df$Time <- "8v12"
ER.12.down_df$Direction <- "Down"

size.df <- do.call(rbind.data.frame, list(ER.1.down_df, ER.1.up_df, ER.2.down_df, ER.2.up_df, ER.3.down_df, ER.3.up_df,
                                          ER.4.down_df, ER.4.up_df, ER.5.down_df, ER.5.up_df, ER.6.down_df, ER.6.up_df,
                                          ER.8.down_df, ER.8.up_df, ER.12.down_df, ER.12.up_df)) 
plot.df <- size.df[size.df$CGI_SIZE.kb != 0 & size.df$Direction == "Up",]

er.size.dens <- ggplot(plot.df[plot.df$Time %in% c("0v1", "1v2"), ],
                       aes(x=CGI_SIZE.kb, colour=Time)) +
  geom_density(size=2, alpha=0.5) + theme_mike() +
  scale_colour_manual(values=c("#D86200", "#D8AA00")) +
  theme(axis.text=element_text(size=16)) +
  guides(colour=FALSE) +
  labs(x="CpG island Size (kb)", y="Density") +
  geom_segment(aes(y=0, yend=1.5,
                   x=median(plot.df$CGI_SIZE.kb[plot.df$Time == "0v1"]),
                   xend=median(plot.df$CGI_SIZE.kb[plot.df$Time == "0v1"])),
               linetype="dashed", colour="#D86200", size=2) +
  geom_segment(aes(y=0, yend=1.5,
                   x=median(plot.df$CGI_SIZE.kb[plot.df$Time == "1v2"]),
                   xend=median(plot.df$CGI_SIZE.kb[plot.df$Time == "1v2"])),
               linetype="dashed", colour="#D8AA00", size=2) +
  scale_x_continuous(limits=c(0, 3), oob=censor)
er.size.dens

ggsave(er.size.dens,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/MCF7_ER-CGIsize_density.png",
       height=4.25, width=4.75, dpi=300)

er.size.box <- ggplot(plot.df[plot.df$Time %in% c("0v1", "1v2"), ],
                       aes(x=Time, y=CGI_SIZE.kb, fill=Time)) +
  geom_boxplot() +
  theme_mike() +
  scale_fill_manual(values=c("#D86200", "#D8AA00")) +
  theme(axis.text=element_text(size=16), axis.title.x=element_blank()) +
  guides(fill=FALSE) +
  labs(y="CpG island Size (kb)", x="Time Comparison") +
  scale_y_continuous(limits=c(0, 3), oob=censor)
er.size.box

ggsave(er.size.box,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/MCF7_ER-CGIsize_boxplot.png",
       height=4.25, width=3.25, dpi=300)

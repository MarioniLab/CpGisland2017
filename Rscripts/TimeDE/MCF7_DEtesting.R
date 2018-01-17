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

# ## 3hr vs 4hr
# de <- results(dds, contrast=c("Timepoint", "4hrs", "3hrs"))
# res34 <- as.data.frame(lfcShrink(dds, res=de, contrast=c("Timepoint", "4hrs", "3hrs")))
# 
# res34$padj[is.na(res34$padj)] <- 1.0
# res34$Sig <- as.factor(as.numeric(res34$padj <= 0.01))
# 
# res34$gene_id <- rownames(res34)
# res34.merge <- merge(res34,  gene_symbol,
#                      by.x='gene_id', by.y='entrezgene', all.x=TRUE)
# 
# ## 4hr vs 5hr
# de <- results(dds, contrast=c("Timepoint", "5hrs", "4hrs"))
# res45 <- as.data.frame(lfcShrink(dds, res=de, contrast=c("Timepoint", "4hrs", "5hrs")))
# 
# res45$padj[is.na(res45$padj)] <- 1.0
# res45$Sig <- as.factor(as.numeric(res45$padj <= 0.01))
# 
# res45$gene_id <- rownames(res45)
# res45.merge <- merge(res45,  gene_symbol,
#                      by.x='gene_id', by.y='entrezgene', all.x=TRUE)
# 
# ## 5hr vs 6hr
# de <- results(dds, contrast=c("Timepoint", "5hrs", "6hrs"))
# res56 <- as.data.frame(lfcShrink(dds, res=de, contrast=c("Timepoint", "5hrs", "6hrs")))
# 
# res56$padj[is.na(res56$padj)] <- 1.0
# res56$Sig <- as.factor(as.numeric(res56$padj <= 0.01))
# 
# res56$gene_id <- rownames(res56)
# res56.merge <- merge(res56,  gene_symbol,
#                      by.x='gene_id', by.y='entrezgene', all.x=TRUE)
# 
# 
# ## 6hr vs 8hr
# de <- results(dds, contrast=c("Timepoint", "6hrs", "8hrs"))
# res68 <- as.data.frame(lfcShrink(dds, res=de, contrast=c("Timepoint", "6hrs", "8hrs")))
# 
# res68$padj[is.na(res68$padj)] <- 1.0
# res68$Sig <- as.factor(as.numeric(res68$padj <= 0.01))
# 
# res68$gene_id <- rownames(res68)
# res68.merge <- merge(res68,  gene_symbol,
#                      by.x='gene_id', by.y='entrezgene', all.x=TRUE)
# 
# 
# ## 8hr vs 12hr
# de <- results(dds, contrast=c("Timepoint", "8hrs", "12hrs"))
# res812 <- as.data.frame(lfcShrink(dds, res=de, contrast=c("Timepoint", "8hrs", "12hrs")))
# 
# res812$padj[is.na(res812$padj)] <- 1.0
# res812$Sig <- as.factor(as.numeric(res812$padj <= 0.01))
# 
# res812$gene_id <- rownames(res812)
# res812.merge <- merge(res812,  gene_symbol,
#                       by.x='gene_id', by.y='entrezgene', all.x=TRUE)

################################################################################################
#### use a binomial test to find differences in the CpG island size ranks between timepoints ###
################################################################################################
ER.0v1.merge <- merge(res01.merge, human.genomic.features, by.x='ensembl_gene_id', by.y='GENE')
ER.1v2.merge <- merge(res12.merge, human.genomic.features, by.x='ensembl_gene_id', by.y='GENE')

# select only the up-regulated genes
ER.0v1.merge.cgi <- ER.0v1.merge[ER.0v1.merge$CGI_SIZE.kb != 0 &
                                   !is.na(ER.0v1.merge$CGI_SIZE.kb) &
                                   ER.0v1.merge$log2FoldChange > 0, ]
ER.1v2.merge.cgi <- ER.1v2.merge[ER.1v2.merge$CGI_SIZE.kb != 0 & 
                                   !is.na(ER.1v2.merge$CGI_SIZE.kb) &
                                   ER.1v2.merge$log2FoldChange > 0, ]

# order based on t-statistic
ER.res01.size_rank <- ER.0v1.merge.cgi$CGI_SIZE.kb[order(ER.0v1.merge.cgi$stat, decreasing=TRUE)]
ER.res26.size_rank <- ER.1v2.merge.cgi$CGI_SIZE.kb[order(ER.1v2.merge.cgi$stat, decreasing=TRUE)]

ER.sign.size_rank <- as.numeric(ER.res01.size_rank < ER.res26.size_rank)
binom.test(sum(ER.sign.size_rank), n=length(ER.sign.size_rank),
           alternative="greater")


# get CGI sizes for up-regulated genes at each time point
# top 250 genes from each only
ER.up.0v1 <- ER.0v1.merge.cgi$CGI_SIZE.kb
ER.up.0v1 <- ER.up.0v1[order(ER.0v1.merge.cgi$stat, decreasing=TRUE)][1:250]

ER.up.1v2 <- ER.1v2.merge.cgi$CGI_SIZE.kb
ER.up.1v2 <- ER.up.1v2[order(ER.1v2.merge.cgi$stat, decreasing=TRUE)][1:250]

ER.size_list <- list("0v1"=cbind(ER.up.0v1, rep("0v1", length(ER.up.0v1))),
                           "1v2"=cbind(ER.up.1v2, rep("1v2", length(ER.up.1v2))))

ER.size.df <- data.frame(do.call(rbind, ER.size_list))
colnames(ER.size.df) <- c("CGI_SIZE.kb", "Comparison")
ER.size.df$CGI_SIZE.kb <- as.numeric(as.character(ER.size.df$CGI_SIZE.kb))
ER.size.df$Comparison <- factor(ER.size.df$Comparison,
                                      levels=c("0v1", "1v2"),
                                      labels=c("0v1", "1v2"))

#################
### plotting ####
#################
er.size.dens <- ggplot(ER.size.df[ER.size.df$Comparison %in% c("0v1", "1v2"), ],
                       aes(x=CGI_SIZE.kb, colour=Comparison)) +
  geom_density(size=2, alpha=0.5) + theme_mike() +
  scale_colour_manual(values=c("#D86200", "#D8AA00")) +
  theme(axis.text=element_text(size=16)) +
  guides(colour=FALSE) +
  labs(x="CpG island Size (kb)", y="Density") +
  geom_segment(aes(y=0, yend=1.5,
                   x=median(ER.size.df$CGI_SIZE.kb[ER.size.df$Comparison == "0v1"]),
                   xend=median(ER.size.df$CGI_SIZE.kb[ER.size.df$Comparison == "0v1"])),
               linetype="dashed", colour="#D86200", size=2) +
  geom_segment(aes(y=0, yend=1.5,
                   x=median(ER.size.df$CGI_SIZE.kb[ER.size.df$Comparison == "1v2"]),
                   xend=median(ER.size.df$CGI_SIZE.kb[ER.size.df$Comparison == "1v2"])),
               linetype="dashed", colour="#D8AA00", size=2) +
  scale_x_continuous(limits=c(0, 3), oob=censor)
er.size.dens

ggsave(er.size.dens,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/MCF7_ER-CGIsize_density.png",
       height=4.25, width=4.75, dpi=300)

er.size.box <- ggplot(ER.size.df[ER.size.df$Comparison %in% c("0v1", "1v2"), ],
                      aes(x=Comparison, y=CGI_SIZE.kb, fill=Comparison)) +
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


# plot heatmap of ranked CpG island size
# to illustrate how I test for an enrichment of short CpG island genes


#############
## ERalpha ##
#############
ER.diff.df <- do.call(cbind.data.frame,
                       list("SizeDiff"=ER.sign.size_rank,
                            "Size"=ER.res26.size_rank))
ER.diff.df$Comparison <- "Test"

ER.null.df <- do.call(cbind.data.frame,
                       list("SizeDiff"=rbinom(n=length(ER.res26.size_rank),
                                              size=1, prob=0.5),
                            "Size"=ER.res01.size_rank))
ER.null.df$Comparison <- "Null"

ER.binom.df <- do.call(rbind.data.frame,
                        list("0v1"=ER.diff.df,
                             "null"=ER.null.df))

ER.binom.df$Comparison <- factor(ER.binom.df$Comparison,
                                  levels=c("Null", "Test"),
                                  labels=c("Null", "Test"))

ER.0v1.heat <- ggplot(ER.binom.df,
                       aes(x=Size, y=Comparison)) +
  geom_tile(aes(fill=SizeDiff)) +
  scale_fill_gradient(low="grey", high="darkred") +
  scale_x_continuous(limits=c(0, 3), oob=censor) +
  theme_mike()  +
  theme(panel.grid=element_blank()) +
  labs(x="CpG island size (kb)", y="Comparison") +
  guides(fill=FALSE)
ER.0v1.heat


ggsave(ER.0v1.heat,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/MCF7_ER-binom_heat.png",
       width=8.25, height=2.25, dpi=300)


fc.01.df <- do.call(cbind.data.frame,
                    list("FC"=ER.0v1.merge.cgi$log2FoldChange,
                         "Size"=ER.0v1.merge.cgi$CGI_SIZE.kb))
fc.01.df$Comparison <- "0v1"

fc.121.df <- do.call(cbind.data.frame,
                    list("FC"=ER.1v2.merge.cgi$log2FoldChange,
                         "Size"=ER.1v2.merge.cgi$CGI_SIZE.kb))
fc.121.df$Comparison <- "1v2"

fc.df <- do.call(rbind.data.frame,
                 list("0v1"=fc.01.df,
                      "1v2"=fc.121.df))

ggplot(fc.df,
       aes(y=FC, x=Comparison)) +
  geom_tile(aes(colour=FC), size=1) +
  scale_colour_gradient(low="lightyellow", high="darkred") +
  theme_mike()  +
  theme(panel.grid=element_blank()) +
  labs(y=expression(paste("Log"[2], " fold change")),
       x="Comparison") +
  guides(fill=FALSE)





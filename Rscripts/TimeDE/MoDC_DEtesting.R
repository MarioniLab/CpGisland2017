#########################################################################################
### GSE84865 - Ebola virus glycoprotein challenge of monocyte-derived Dendritic Cells ###
#########################################################################################
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
human.genomic.features$CGI_SIZE.kb <- human.genomic.features$CGI_SIZE/1000

# these look like expected counts from Kallisto or Salmon.  No details in the linked MS however.
gse84865 <- read.table("~/Dropbox/Noise_genomics/Timeseries/GSE84865_genes_expression_expected_count.tsv",
                       sep="\t", h=TRUE, stringsAsFactors=FALSE)

# gene IDs are gene symbols
ensembl <- useEnsembl(biomart='ensembl', dataset='hsapiens_gene_ensembl')

gene_symbol <- getBM(attributes=c('external_gene_name', 'entrezgene', 'ensembl_gene_id'),
                     filters='external_gene_name', values=gse84865$gene,
                     mart=ensembl)

gse84865.merge <- merge(gse84865, gene_symbol, by.x='gene', by.y='external_gene_name')
gse84865.merge <- gse84865.merge[!duplicated(gse84865.merge$ensembl_gene_id), ]
gse84865.counts <- gse84865.merge[, colnames(gse84865.merge)[!grepl(colnames(gse84865.merge),
                                                                    pattern="(gene|transcript)")]]

gse84865.counts <- as.data.frame(apply(gse84865.counts, 2, FUN=function(X) floor(as.numeric(X))))
rownames(gse84865.counts) <- gse84865.merge$ensembl_gene_id

# get metadata info from colnames
colnames(gse84865.counts) <- gsub(colnames(gse84865.counts), pattern="([A-Z])_([A-Z])", replacement="\\1.\\2")
colnames(gse84865.counts) <- gsub(colnames(gse84865.counts), pattern="([0-9])_([0-9])", replacement="\\1_WT_\\2")

noms <- strsplit(colnames(gse84865.counts), split="_", fixed=TRUE)
donor <- unlist(lapply(noms, FUN=function(X) paste0(X[1])))
genotype <- unlist(lapply(noms, FUN=function(X) paste0(X[2])))
timepoint <- unlist(lapply(noms, FUN=function(X) paste0(X[3])))

gse84865.meta <- do.call(cbind.data.frame,
                         list("Sample"=colnames(gse84865.counts),
                              "Donor"=donor,
                              "Genotype"=genotype,
                              "Timepoint"=timepoint))
gse84865.meta$Sample <- as.character(gse84865.meta$Sample)
rownames(gse84865.meta) <- gse84865.meta$Sample

# remove sparse genes, mean > 5 reads
nonZero <- rowMeans(gse84865.counts) > 5
gse84865.nz <- gse84865.counts[nonZero, ]

# size factor normalisation
# keep WT samples only
#wt.samps <- gse84865.meta$Sample[gse84865.meta$Genotype == "WT"]

dds <- DESeqDataSetFromMatrix(countData=floor(gse84865.nz),
                              colData=gse84865.meta,
                              design= ~ Donor + Genotype + Timepoint)

# run DESeq2
dds <- DESeq(dds)

## 0hr vs 1hr
de <- results(dds, contrast=c("Timepoint", "1hr", "0hr"))
res01 <- as.data.frame(lfcShrink(dds, res=de, contrast=c("Timepoint", "1hr", "0hr")))

res01$padj[is.na(res01$padj)] <- 1.0
res01$Sig <- as.factor(as.numeric(res01$padj <= 0.1))

res01$gene_id <- rownames(res01)
res01.merge <- merge(res01,  human.genomic.features,
                     by.x='gene_id', by.y='GENE', all.x=TRUE)

## 1hr vs 2hr
de <- results(dds, contrast=c("Timepoint", "2hr", "1hr"))
res12 <- as.data.frame(lfcShrink(dds, res=de, contrast=c("Timepoint", "2hr", "1hr")))

res12$padj[is.na(res12$padj)] <- 1.0
res12$Sig <- as.factor(as.numeric(res12$padj <= 0.1))

res12$gene_id <- rownames(res12)
res12.merge <- merge(res12,  human.genomic.features,
                     by.x='gene_id', by.y='GENE', all.x=TRUE)


################################################################################################
#### use a binomial test to find differences in the CpG island size ranks between timepoints ###
################################################################################################
# select only the up-regulated genes
modc.0v1.merge.cgi <- res01.merge[res01.merge$CGI_SIZE.kb != 0 &
                                    !is.na(res01.merge$CGI_SIZE.kb) &
                                    res01.merge$log2FoldChange > 0, ]
modc.1v2.merge.cgi <- res12.merge[res12.merge$CGI_SIZE.kb != 0 & 
                                    !is.na(res12.merge$CGI_SIZE.kb) &
                                    res12.merge$log2FoldChange > 0, ]

# order based on t-statistic
modc.res01.size_rank <- modc.0v1.merge.cgi$CGI_SIZE.kb[order(modc.0v1.merge.cgi$stat, decreasing=TRUE)]
modc.res26.size_rank <- modc.1v2.merge.cgi$CGI_SIZE.kb[order(modc.1v2.merge.cgi$stat, decreasing=TRUE)]

modc.sign.size_rank <- as.numeric(modc.res01.size_rank < modc.res26.size_rank)
binom.test(sum(modc.sign.size_rank), n=length(modc.sign.size_rank),
           alternative="greater")

# get CGI sizes for up-regulated genes at each time point
# top 250 genes from each only
modc.up.0v1 <- modc.0v1.merge.cgi$CGI_SIZE.kb
modc.up.0v1 <- modc.up.0v1[order(modc.0v1.merge.cgi$stat, decreasing=TRUE)][1:250]

modc.up.1v2 <- modc.1v2.merge.cgi$CGI_SIZE.kb
modc.up.1v2 <- modc.up.1v2[order(modc.1v2.merge.cgi$stat, decreasing=TRUE)][1:250]

modc.size_list <- list("0v1"=cbind(modc.up.0v1, rep("0v1", length(modc.up.0v1))),
                       "1v2"=cbind(modc.up.1v2, rep("1v2", length(modc.up.1v2))))

modc.size.df <- data.frame(do.call(rbind, modc.size_list))
colnames(modc.size.df) <- c("CGI_SIZE.kb", "Comparison")
modc.size.df$CGI_SIZE.kb <- as.numeric(as.character(modc.size.df$CGI_SIZE.kb))
modc.size.df$Comparison <- factor(modc.size.df$Comparison,
                                  levels=c("0v1", "1v2"),
                                  labels=c("0v1", "1v2"))

#################
### plotting ####
#################
modc.size.dens <- ggplot(modc.size.df[modc.size.df$Comparison %in% c("0v1", "1v2"), ],
                         aes(x=CGI_SIZE.kb, colour=Comparison)) +
  geom_density(size=2, alpha=0.5) + theme_mike() +
  scale_colour_manual(values=c("#D86200", "#D8AA00")) +
  theme(axis.text=element_text(size=16)) +
  guides(colour=FALSE) +
  labs(x="CpG island Size (kb)", y="Density") +
  geom_segment(aes(y=0, yend=1.5,
                   x=median(modc.size.df$CGI_SIZE.kb[modc.size.df$Comparison == "0v1"]),
                   xend=median(modc.size.df$CGI_SIZE.kb[modc.size.df$Comparison == "0v1"])),
               linetype="dashed", colour="#D86200", size=2) +
  geom_segment(aes(y=0, yend=1.5,
                   x=median(modc.size.df$CGI_SIZE.kb[modc.size.df$Comparison == "1v2"]),
                   xend=median(modc.size.df$CGI_SIZE.kb[modc.size.df$Comparison == "1v2"])),
               linetype="dashed", colour="#D8AA00", size=2) +
  scale_x_continuous(limits=c(0, 3), oob=censor)

ggsave(modc.size.dens,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/MODC_EBOLA-CGIsize_density.png",
       height=4.25, width=4.75, dpi=300)


modc.size.box <- ggplot(modc.size.df[modc.size.df$Comparison %in% c("0v1", "1v2"), ],
                        aes(x=Comparison, y=CGI_SIZE.kb, fill=Comparison)) +
  geom_boxplot() +
  theme_mike() +
  scale_fill_manual(values=c("#D86200", "#D8AA00")) +
  theme(axis.text=element_text(size=16), axis.title.x=element_blank()) +
  guides(fill=FALSE) +
  labs(y="CpG island Size (kb)", x="Time Comparison") +
  scale_y_continuous(limits=c(0, 3), oob=censor)

ggsave(modc.size.box,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/MODC_EBOLA-CGIsize_boxplot.png",
       height=4.25, width=3.25, dpi=300)








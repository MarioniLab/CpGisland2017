## read in genomic feature information
## Investigating genomic factors that influence gene expression noise
library(ggplot2)
library(reshape2)
library(Rtsne)
library(glmnet)
library(biomaRt)
source("~/Dropbox/R_sessions/GGMike/palette_256.R")
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

# need the summarised table, not this one
gene.feature <- read.table("~/Dropbox/ENSEMBL/mm10/mm10_exon-stats.tsv",
                           sep="\t", h=T, stringsAsFactors=F)
phast.feature <- read.table("~/Dropbox/UCSC/mm10/mm10_phastCons60way-summary.tsv",
                            sep="\t", h=T, stringsAsFactors=F)
colnames(phast.feature) <- c("GENE", "MED_PHAST", "NALIGN_PHAST", "SUM_PHAST")

gc.feature <- read.table("~/Dropbox/ENSEMBL/mm10/mm10_ensembl86_transcripts_promoter-gc.tsv",
                         h=F, sep="\t", stringsAsFactors=F)
colnames(gc.feature) <- c("GENE", "GC")

cpg.feature <- read.table("~/Dropbox/UCSC/mm10/mm10_CpG-Stats.tsv",
                          h=T, sep="\t", stringsAsFactors=F)
cpg.feature <- cpg.feature[!duplicated(cpg.feature$GENE), c("GENE", "cpg_GCNUM",
                                                            "cpg_Overlap", "cpg_RATIO",
                                                            "N_CpG")]
cpg.feature$N_CpG[cpg.feature$N_CpG >= 1] <- 1

cpg.bed <- read.table("~/Dropbox/ENSEMBL/mm10/mm10_ensembl86_transcripts_promoter-CGI_intersect.bed.gz",
                      h=F, sep="\t", stringsAsFactors=F)
colnames(cpg.bed) <- c("CHR", "START", "END", "GENE", "SCORE", "STRAND", "CGI_CHR", "CGI_START", "CGI_END", "CGI_NAME", "CGI_cpgnum", "CGI_STRAND", "CGI_OVERLAP")
cpg.bed$CGI_SIZE <- abs(cpg.bed$CGI_END - cpg.bed$CGI_START)

cpg.feature <- merge(cpg.feature, cpg.bed[, c("GENE", "CGI_SIZE")], by="GENE", all.x=TRUE)

mesc.sp1 <- read.table("~/Dropbox/ENSEMBL/mm10/mm10_ensembl86_transcripts_promoter-SP1_JASPAR.bed.gz",
                       h=F, sep="\t", stringsAsFactors=F)
colnames(mesc.sp1) <- c("CHR", "START", "END", "GENE", "SCORE", "STRAND", "SP1")

mesc.tbp <- read.table("~/Dropbox/ENSEMBL/mm10/mm10_ensembl86_transcripts_promoter-TBP_JASPAR.bed.gz",
                       h=F, sep="\t", stringsAsFactors=F)
colnames(mesc.tbp) <- c("CHR", "START", "END", "GENE", "SCORE", "STRAND", "TBP")

mesc.nmi <- read.table("~/Dropbox/ENSEMBL/mm10/mm10_ensembl86_transcripts_promoter-mESC_NMI.bed.gz",
                       sep="\t", h=F, stringsAsFactors=F)
colnames(mesc.nmi) <- c("CHR", "START", "END", "GENE", "SCORE", "STRAND", "NMI")

genomic.features <- merge(gc.feature, cpg.feature, by='GENE')
genomic.features <- merge(genomic.features, mesc.sp1[, c("GENE", "SP1")], by='GENE')
genomic.features <- merge(genomic.features, mesc.tbp[, c("GENE", "TBP")], by='GENE')
genomic.features <- merge(genomic.features, phast.feature, by='GENE')
genomic.features <- merge(genomic.features, gene.feature, by='GENE')
genomic.features <- merge(genomic.features, mesc.nmi[, c("GENE", "NMI")], by='GENE')
genomic.features <- merge(genomic.features, chrom.cast, by='GENE')

# drop duplicates
mouse.genomic.features <- genomic.features[!duplicated(genomic.features$GENE), ]

tra.df <- read.table("~/Dropbox/fantom5/Tau_TRA_level3.tsv",
                     sep="\t", h=T, stringsAsFactors=F)
mouse.consistent.genes <- tra.df$ensembl.gene[tra.df$tau <= 0.4]
mouse.tra.genes <- tra.df$ensembl.gene[tra.df$tau >= 0.8]

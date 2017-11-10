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
gene.feature <- read.table("~/Dropbox/ENSEMBL/hg38/hg38_exon-stats.tsv",
                           sep="\t", h=T, stringsAsFactors=F)

gc.feature <- read.table("~/Dropbox/ENSEMBL/hg38/hg38_ensembl86_transcripts_promoter-gc.tsv",
                         h=F, sep="\t", stringsAsFactors=F)
colnames(gc.feature) <- c("GENE", "GC")

cpg.size <- read.table("~/Dropbox/ENSEMBL/hg19/hg19_ensembl86_transcripts_promoter-CpG_size.bed.gz",
                       h=F, stringsAsFactors=F, sep="\t")
colnames(cpg.size) <- c("CGI_CHR", "CGI_START", "CGI_END", "CGI_ID", "CGI_SIZE", "CGI_STRAND",
                        "CHR", "START", "END", "GENE", "SCORE", "STRAND", "OVERLAP")

# there are considerable inconsistencies between hg19 CGIs and hg38 ones
cpg.feature <- read.table("~/Dropbox/UCSC/hg38/hg38_CpG-Stats.tsv",
                          h=T, sep="\t", stringsAsFactors=F)
cpg.feature <- cpg.feature[!duplicated(cpg.feature$GENE), c("GENE", "cpg_GCNUM",
                                                            "cpg_Overlap", "cpg_RATIO",
                                                            "N_CpG")]
cpg.feature$N_CpG[cpg.feature$N_CpG >= 1] <- 1

hesc.sp1 <- read.table("~/Dropbox/ENSEMBL/hg38/hg38_ensembl86_transcripts_promoter-SP1_JASPAR.bed.gz",
                       h=F, sep="\t", stringsAsFactors=F)
colnames(hesc.sp1) <- c("CHR", "START", "END", "GENE", "SCORE", "STRAND", "SP1")

hesc.tbp <- read.table("~/Dropbox/ENSEMBL/hg38/hg38_ensembl86_transcripts_promoter-TBP_JASPAR.bed.gz",
                       h=F, sep="\t", stringsAsFactors=F)
colnames(hesc.tbp) <- c("CHR", "START", "END", "GENE", "SCORE", "STRAND", "TBP")

genomic.features <- merge(gc.feature, cpg.feature, by='GENE')
genomic.features <- merge(genomic.features, hesc.sp1[, c("GENE", "SP1")], by='GENE')
genomic.features <- merge(genomic.features, hesc.tbp[, c("GENE", "TBP")], by='GENE')
genomic.features <- merge(genomic.features, gene.feature, by='GENE')
genomic.features <- merge(genomic.features, cpg.size[, c("GENE", "CGI_SIZE", "CGI_ID")], by='GENE', all.x=TRUE)

genomic.features$N_CpG <- as.numeric(!is.na(genomic.features$CGI_ID))

# drop duplicates
human.genomic.features <- genomic.features[!duplicated(genomic.features$GENE), ]
human.genomic.features$CGI_ID[is.na(human.genomic.features$CGI_ID)] <- "NONE"
human.genomic.features$CGI_SIZE[is.na(human.genomic.features$CGI_SIZE)] <- 0

human.genomic.features$CGI_SIZE.kb <- human.genomic.features$CGI_SIZE/1000

# split NMIs into size bins on quartiles
human.genomic.features$CGI_SIZE.group <- as.character(cut(human.genomic.features$CGI_SIZE.kb,
                                                          breaks=quantile(human.genomic.features$CGI_SIZE.kb[human.genomic.features$CGI_SIZE.kb != 0],
                                                                          probs=c(0, 0.25, 0.5, 0.75, 1.0))))
human.genomic.features$CGI_SIZE.group[(human.genomic.features$CGI_SIZE.kb <= 0.201)] <- "(0,0.201]"
human.genomic.features$CGI_SIZE.group[is.na(human.genomic.features$CGI_SIZE.group) |human.genomic.features$CGI_SIZE.kb == 0] <- "Absent"

# merge the V.Short and Short CGIs
human.genomic.features$CGI_SIZE.group[human.genomic.features$CGI_SIZE.group == "(0,0.201]"] <- "(0.201,0.494]"
human.genomic.features$CGI_SIZE.group <- factor(human.genomic.features$CGI_SIZE.group,
                                                labels=c("Absent", "Short", "Short.Mid", "Long.Mid", "Long"),
                                                levels=c("Absent", "(0.201,0.494]", "(0.494,0.787]", "(0.787,1.2]", "(1.2,40.1]"))

tra.df <- read.table("~/Dropbox/fantom5/hg19/Human-Tau_TRA_level3.tsv",
                     sep="\t", h=T, stringsAsFactors=F)
human.consistent.genes <- tra.df$ensembl.gene[tra.df$tau <= 0.4]
human.tra.genes <- tra.df$ensembl.gene[tra.df$tau >= 0.8]
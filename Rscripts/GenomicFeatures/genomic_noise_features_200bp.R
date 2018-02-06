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

cpg.bed <- read.table("~/Dropbox/ENSEMBL/mm10/promoters/promoter_200/mm10_ensembl86_transcripts_promoter200-CpGisland.bed.gz",
                      h=F, sep="\t", stringsAsFactors=F)
colnames(cpg.bed) <- c("CHR", "START", "END", "GENE", "SCORE", "STRAND", "CGI_CHR", "CGI_START", "CGI_END", "CGI_NAME", "CGI_LEN", "CGI_cpgnum", "CGI_GCNUM",
                       "CGI_perCpG", "CGI_perGC", "CGI_RATIO", "OVERLAP")
cpg.bed$CGI_SIZE <- abs(cpg.bed$CGI_END - cpg.bed$CGI_START)

cpg.feature <- merge(cpg.feature, cpg.bed[, c("GENE", "CGI_SIZE")], by="GENE", all.x=TRUE)

mesc.sp1 <- read.table("~/Dropbox/ENSEMBL/mm10/promoters/promoter_200/mm10_ensembl86_transcripts_promoter200-SP1.bed.gz",
                       h=F, sep="\t", stringsAsFactors=F)
colnames(mesc.sp1) <- c("CHR", "START", "END", "GENE", "SCORE", "STRAND", "SP1_CHR", "SP1_START", "SP1_END", "MOTIF", "SIGNAL", "SP1_STRAND", "OVERLAP")

# summarise the number of SP1 motifs per gene
mesc.sp1.count <- do.call(cbind.data.frame,
                          list("GENE"=unique(mesc.sp1$GENE),
                               "SP1"=as.numeric(by(mesc.sp1$SP1_CHR, INDICES=mesc.sp1$GENE, FUN=function(X) sum(X != '.')))))

mesc.tbp <- read.table("~/Dropbox/ENSEMBL/mm10/promoters/promoter_200/mm10_ensembl86_transcripts_promoter200-TBP.bed.gz",
                       h=F, sep="\t", stringsAsFactors=F)
colnames(mesc.tbp) <- c("CHR", "START", "END", "GENE", "SCORE", "STRAND", "TBP_CHR", "TBP_START", "TBP_END", "MOTIF", "SIGNAL", "TBP_STRAND", "OVERLAP")
# summarise the number of SP1 motifs per gene
mesc.tbp.count <- do.call(cbind.data.frame,
                          list("GENE"=unique(mesc.tbp$GENE),
                               "TBP"=as.numeric(by(mesc.tbp$TBP_CHR, INDICES=mesc.tbp$GENE, FUN=function(X) sum(X != '.')))))

# bring in additional info from EPDnew
epd.motifs <- read.table("~/Dropbox/ENSEMBL/EPDnew/mm9_promoter_motifs.txt",
                         h=TRUE, sep="\t", stringsAsFactors=FALSE, comment.char='')
colnames(epd.motifs) <- c("EPD", "TATA.box", "Inr", "CCAAT.box", "GC.box")

epd.ensembl <- read.table("~/Dropbox/ENSEMBL/EPDnew/mm9_promoter_ensembl.txt",
                          h=FALSE, sep="\t", stringsAsFactors=FALSE)
colnames(epd.ensembl) <- c("EPD", "GENE")

epd.merge <- merge(epd.ensembl, epd.motifs, by='EPD', all=TRUE)

genomic.features <- merge(gc.feature, cpg.feature, by='GENE')
genomic.features <- merge(genomic.features, mesc.sp1.count, by='GENE')
genomic.features <- merge(genomic.features, mesc.tbp.count, by='GENE')
genomic.features <- merge(genomic.features, phast.feature, by='GENE')
genomic.features <- merge(genomic.features, gene.feature, by='GENE')
#genomic.features <- merge(genomic.features, mesc.nmi[, c("GENE", "NMI")], by='GENE')
genomic.features <- merge(genomic.features, epd.merge, by='GENE')
genomic.features <- merge(genomic.features, chrom.cast, by='GENE')

# drop duplicates
mouse.genomic.features <- genomic.features[!duplicated(genomic.features$GENE), ]

tra.df <- read.table("~/Dropbox/fantom5/Tau_TRA_level3.tsv",
                     sep="\t", h=T, stringsAsFactors=F)
mouse.consistent.genes <- tra.df$ensembl.gene[tra.df$tau <= 0.4]
mouse.tra.genes <- tra.df$ensembl.gene[tra.df$tau >= 0.8]

## Investigating genomic factors that influence gene expression noise
library(ggplot2)
library(reshape2)
library(Rtsne)
library(glmnet)
source("~/Dropbox/R_sessions/GGMike/palette_256.R")
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

liver.cells <- read.table("~/Dropbox/Liver/liver_SFnorm.tsv",
                          sep="\t", h=T, stringsAsFactors=F)
rownames(liver.cells) <- liver.cells$gene_id
liver.cells <- liver.cells[grepl(rownames(liver.cells), pattern="ENS"), ]

genomic.features <- read.table("~/Dropbox/ENCODE/Liver/Liver_ENCODE-counts_matrix.tsv",
                               h=T, sep="\t", stringsAsFactors=F)
rownames(genomic.features) <- genomic.features$GENE

# need the summarised table, not this one
gene.feature <- read.table("~/Dropbox/ENSEMBL/mm10_exon-stats.tsv",
                           sep="\t", h=T, stringsAsFactors=F)

phast.feature <- read.table("~/Dropbox/UCSC/mm10_phastCons60way-summary.tsv",
                            sep="\t", h=T, stringsAsFactors=F)
colnames(phast.feature) <- c("GENE", "MED_PHAST", "NALIGN_PHAST", "SUM_PHAST")

gc.feature <- read.table("~/Dropbox/ENSEMBL/mm10_ensembl86_transcripts_promoter-gc.tsv",
                         h=F, sep="\t", stringsAsFactors=F)
colnames(gc.feature) <- c("GENE", "GC")

cpg.feature <- read.table("~/Dropbox/UCSC/mm10_CpG-Stats.tsv",
                          h=T, sep="\t", stringsAsFactors=F)
cpg.feature <- cpg.feature[!duplicated(cpg.feature$GENE), c("GENE", "cpg_GCNUM",
                                                            "cpg_Overlap", "cpg_RATIO",
                                                            "N_CpG")]
tfbs.feature <- read.table("~/Dropbox/ENSEMBL/mm10_ensembl86-liver_RegTFBS.tsv.gz",
                           sep="\t", h=F, stringsAsFactors=F)
colnames(tfbs.feature) <- c("CHR", "START", "END", "GENE", "SCORE", "STRAND", "TFBS_COUNT")
tfbs.feature <- tfbs.feature[, c("GENE", "TFBS_COUNT")]

open.counts <- read.table("~/Dropbox/ENSEMBL/mm10_ensembl86-liver_RegOpenChrom.tsv.gz",
                           sep="\t", h=F, stringsAsFactors=F)
colnames(open.counts) <- c("CHR", "START", "END", "GENE", "SCORE", "STRAND", "OC_COUNT")

open.overlap <- read.table("~/Dropbox/ENSEMBL/mm10_ensembl86-liver_RegOpen_overlap.bed.gz",
                           sep="\t", h=F, stringsAsFactors=F)
colnames(open.overlap) <- c("CHR", "START", "END", "GENE", "SCORE", "STRAND",
                            "ov_CHR", "ov_START", "ov_END", "ov_NAME", "ov_SCORE",
                            "ov_STRAND", "OVERLAP")
open.feature <- merge(open.counts, open.overlap,
                      by=c("CHR", "START", "END", "GENE", "SCORE", "STRAND"))
open.feature <- open.feature[, c("GENE", "OC_COUNT", "OVERLAP")]

enhancer.feature <- read.table("~/Dropbox/ENSEMBL/mm10_ensembl86_liver-Enhancer_stats.tsv",
                               h=T, sep="\t", stringsAsFactors=F)

genomic.features <- merge(genomic.features, gc.feature, by='GENE')
genomic.features <- merge(genomic.features, cpg.feature, by='GENE')
genomic.features <- merge(genomic.features, phast.feature, by='GENE')
genomic.features <- merge(genomic.features, gene.feature, by='GENE')
genomic.features <- merge(genomic.features, tfbs.feature, by='GENE')
genomic.features <- merge(genomic.features, open.feature, by='GENE')
genomic.features <- merge(genomic.features, enhancer.feature, by='GENE')

# drop duplicates
genomic.features <- genomic.features[!duplicated(genomic.features$GENE), ]

tra.df <- read.table("~/Dropbox/fantom5/Tau_TRA_level3.tsv",
                     sep="\t", h=T, stringsAsFactors=F)
consistent.genes <- tra.df$ensembl.gene[tra.df$tau <= 0.4]
tra.genes <- tra.df$ensembl.gene[tra.df$tau >= 0.8]
liver.means <- rowMeans(liver.cells[, 1:(dim(liver.cells)[2]-1)])
liver.vars <- apply(liver.cells[, 1:(dim(liver.cells)[2]-1)],
                    1, var)
# create gene expression groups based on average expression over cells
exprs.groups <- as.factor(cut_number(liver.means, n=10))
gene.summary <- as.data.frame(cbind(liver.means, liver.vars, exprs.groups))
colnames(gene.summary) <- c("Mean", "Var", "Group")
gene.summary$CV2 <- (sqrt(gene.summary$Var)/gene.summary$Mean) ** 2
gene.summary$GENE <- rownames(gene.summary)
gene.summary$GeneGroup <- "Misc"
gene.summary$GeneGroup[gene.summary$GENE %in% consistent.genes] <- "Consistent"
gene.summary$GeneGroup[gene.summary$GENE %in% tra.genes] <- "TissueSpec"
gene.summary$GeneGroup <- factor(gene.summary$GeneGroup,
                                 labels=c(1, 2, 3),
                                 levels=c("TissueSpec", "Misc", "Consistent"))

liver.match <- merge(gene.summary, genomic.features,
                     by='GENE')
rownames(liver.match) <- liver.match$GENE
# liver.consist <- na.omit(liver.match[consistent.genes, ])
# liver.consist <- na.omit(liver.match[tra.genes ,])
liver.consist <- liver.match

# linear fit features vs. CV^2
# don't plot features with all zeros
colSums(as.matrix(apply(liver.consist[, c(4, 6:(dim(liver.consist)[2]))],2,as.numeric)))

feature.lmfit <- lm.fit(x=as.matrix(liver.consist[, c(4, 6:(dim(liver.consist)[2]))]),
                        y=liver.consist[, 5])

# plot the features against the CV^2, which ones look like they may have some
# influence on expression noise?
ggplot(liver.consist, aes(x=as.factor(CTCF_0), y=CV2)) + 
  geom_boxplot() + geom_jitter(alpha=0.2) + 
  theme_mike()

ggplot(liver.consist, aes(x=as.factor(ACTIVE_ENHANCER_COUNT), y=CV2)) + 
  geom_boxplot() + geom_jitter(alpha=0.2) + 
  theme_mike()

ggplot(liver.consist, aes(x=log10(ACTIVE_ENHANCER_AVDIST), y=CV2)) + 
  geom_point() +
  theme_mike()

ggplot(liver.consist, aes(x=ACTIVE_ENHANCER_AVSIZE, y=CV2)) + 
  geom_point() +
  theme_mike()

ggplot(liver.consist, aes(x=ACTIVE_ENHANCER_COUNT, y=CV2)) + 
  geom_point() +
  theme_mike()

ggplot(liver.consist, aes(x=as.factor(TFBS_COUNT), y=CV2)) + 
  geom_boxplot() + geom_jitter(alpha=0.2) + 
  theme_mike()

ggplot(liver.consist, aes(x=GeneGroup, y=CV2)) + 
  geom_boxplot() + geom_jitter(alpha=0.2) + 
  theme_mike()

ggplot(liver.consist, aes(x=as.factor(OC_COUNT), y=CV2)) + 
  geom_boxplot() + geom_jitter(alpha=0.2) + 
  theme_mike()

ggplot(liver.consist, aes(x=log10(OVERLAP), y=CV2)) + 
  geom_point() +
  theme_mike()

ggplot(liver.consist, aes(x=log10(EXON_COUNT), y=CV2)) + 
  geom_point() +
  theme_mike()

ggplot(liver.consist, aes(x=log10(EXON_AVLENGTH), y=CV2)) + 
  geom_point() +
  theme_mike()

ggplot(liver.consist, aes(x=log10(EXON_TOTLENGTH), y=CV2)) + 
  geom_point() +
  theme_mike()

ggplot(liver.consist, aes(x=GC, y=CV2)) + 
  geom_point() +
  theme_mike()

ggplot(liver.consist, aes(x=MED_PHAST, y=CV2)) + 
  geom_point() +
  theme_mike()

ggplot(liver.consist, aes(x=NALIGN_PHAST, y=CV2)) + 
  geom_point() +
  theme_mike()

ggplot(liver.consist, aes(x=SUM_PHAST, y=CV2)) + 
  geom_point() +
  theme_mike()

ggplot(liver.consist, aes(x=cpg_GCNUM, y=CV2)) + 
  geom_point() +
  theme_mike()

ggplot(liver.consist, aes(x=as.factor(CTCF_8), y=CV2)) + 
  geom_boxplot() + geom_jitter(alpha=0.2) + 
  theme_mike()

ggplot(liver.consist, aes(x=as.factor(H3K36me3_8), y=CV2)) + 
  geom_boxplot() + geom_jitter(alpha=0.2) + 
  theme_mike()

ggplot(liver.consist, aes(x=as.factor(N_CpG), y=CV2)) + 
  geom_boxplot() + geom_jitter(alpha=0.2) + 
  theme_mike()

ggplot(liver.consist, aes(x=as.factor(Group), y=CV2)) + 
  geom_boxplot() + geom_jitter(alpha=0.2) + 
  theme_mike()

genomic.vars <- paste(c("Group", "GeneGroup", colnames(genomic.features)[2:dim(genomic.features)[2]]),
                      collapse=" + ")
glm.form <- as.formula(paste("CV2", genomic.vars, sep=" ~ "))

feature.glm <- glm(glm.form, data=liver.consist, family=gaussian('identity'))
summary(feature.glm)

# use LASSO to select meaningful features, set alpha=1 to make
# it a full LASSO rather than elastic net
lasso <- glmnet(x=as.matrix(liver.consist[, c(4, 6:(dim(liver.consist)[2]))]),
                y=liver.consist[, 5],
                family="gaussian", alpha=1)
plot(lasso, xvar='dev', label=TRUE)
coef(lasso, s=0.2)

cv.lasso <- cv.glmnet(x=as.matrix(apply(liver.consist[, c(4, 6:(dim(liver.consist)[2]))],2,as.numeric)),
                      y=liver.consist[, 5],
                      family="gaussian", alpha=1,
                      type.measure = "mse", nfolds = 20)
coef(cv.lasso, s=cv.lasso$lambda.1se)
plot(cv.lasso) # the only meaningful feature is mean expression group

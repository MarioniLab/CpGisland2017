## Investigating genomic factors that influence gene expression noise
library(ggplot2)
library(reshape2)
library(Rtsne)
library(glmnet)
source("~/Dropbox/R_sessions/GGMike/palette_256.R")
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

tcell.cells <- read.table("~/Dropbox/Tcell/Tcell_SFnorm.tsv",
                          sep="\t", h=T, stringsAsFactors=F)
rownames(tcell.cells) <- tcell.cells$gene_id
tcell.cells <- tcell.cells[grepl(rownames(tcell.cells), pattern="ENS"), ]

# # need the summarised table, not this one
# gene.feature <- read.table("~/Dropbox/ENSEMBL/mm10_exon-stats.tsv",
#                            sep="\t", h=T, stringsAsFactors=F)
# phast.feature <- read.table("~/Dropbox/UCSC/mm10_phastCons60way-summary.tsv",
#                             sep="\t", h=T, stringsAsFactors=F)
# colnames(phast.feature) <- c("GENE", "MED_PHAST", "NALIGN_PHAST", "SUM_PHAST")
# 
# gc.feature <- read.table("~/Dropbox/ENSEMBL/mm10_ensembl86_transcripts_promoter-gc.tsv",
#                          h=F, sep="\t", stringsAsFactors=F)
# colnames(gc.feature) <- c("GENE", "GC")
# 
# cpg.feature <- read.table("~/Dropbox/UCSC/mm10_CpG-Stats.tsv",
#                           h=T, sep="\t", stringsAsFactors=F)
# cpg.feature <- cpg.feature[!duplicated(cpg.feature$GENE), c("GENE", "cpg_GCNUM",
#                                                             "cpg_Overlap", "cpg_RATIO",
#                                                             "N_CpG")]
# 
# genomic.features <- merge(gc.feature, cpg.feature, by='GENE')
# genomic.features <- merge(genomic.features, phast.feature, by='GENE')
# genomic.features <- merge(genomic.features, gene.feature, by='GENE')
# 
# # drop duplicates
# genomic.features <- genomic.features[!duplicated(genomic.features$GENE), ]
# 
# tra.df <- read.table("~/Dropbox/fantom5/Tau_TRA_level3.tsv",
#                      sep="\t", h=T, stringsAsFactors=F)
# consistent.genes <- tra.df$ensembl.gene[tra.df$tau <= 0.4]
# tra.genes <- tra.df$ensembl.gene[tra.df$tau >= 0.8]
tcell.means <- rowMeans(tcell.cells[, 1:(dim(tcell.cells)[2]-1)])
tcell.vars <- apply(tcell.cells[, 1:(dim(tcell.cells)[2]-1)],
                    1, var)
tcell.median <- apply(tcell.cells[, 1:(dim(tcell.cells)[2]-1)],
                      1, median)
tcell.mad <- apply(tcell.cells[, 1:(dim(tcell.cells)[2]-1)],
                      1, mad)

# create gene expression groups based on average expression over cells
tcell.exprs.groups <- as.factor(cut_number(tcell.means, n=10))
tcell.gene.summary <- as.data.frame(cbind(tcell.means, tcell.vars, tcell.median, tcell.mad,
                                          tcell.exprs.groups))
colnames(tcell.gene.summary) <- c("Mean", "Var", "Median", "MAD", "Group")
tcell.gene.summary$CV2 <- tcell.gene.summary$Var/(tcell.gene.summary$Mean ** 2)
tcell.gene.summary$GENE <- rownames(tcell.gene.summary)
# gene.summary$GeneGroup <- "Misc"
# gene.summary$GeneGroup[gene.summary$GENE %in% consistent.genes] <- "Consistent"
# gene.summary$GeneGroup[gene.summary$GENE %in% tra.genes] <- "TissueSpec"
# gene.summary$GeneGroup <- factor(gene.summary$GeneGroup,
#                                  labels=c(1, 2, 3),
#                                  levels=c("TissueSpec", "Misc", "Consistent"))
# 
# tcell.match <- merge(gene.summary, genomic.features,
#                      by='GENE')
# rownames(tcell.match) <- tcell.match$GENE
# # tcell.consist <- na.omit(tcell.match[consistent.genes, ])
# # tcell.consist <- na.omit(tcell.match[tra.genes ,])
# tcell.consist <- tcell.match
# 
# # linear fit features vs. CV^2
# # don't plot features with all zeros
# colSums(as.matrix(apply(tcell.consist[, c(4, 6:(dim(tcell.consist)[2]))],2,as.numeric)))
# 
# feature.lmfit <- lm.fit(x=as.matrix(tcell.consist[, c(4, 6:(dim(tcell.consist)[2]))]),
#                         y=tcell.consist[, 5])
# 
# # plot the features against the CV^2, which ones look like they may have some
# # influence on expression noise?
# p_geneGroup <- ggplot(tcell.consist, aes(x=GeneGroup, y=CV2, fill=GeneGroup)) + 
#   geom_boxplot() + geom_jitter(alpha=0.2) + 
#   theme_mike() + scale_colour_Publication()
# ggsave(p_geneGroup, filename="~/Dropbox/Noise_genomics/Tcell-GeneGroup_vs_CV2.png",
#        width=4, height=4*1.618, dpi=90)
# 
# p_exonLen <- ggplot(tcell.consist, aes(x=log10(EXON_TOTLENGTH), y=CV2)) + 
#   geom_point() +
#   theme_mike()
# ggsave(p_exonLen, filename="~/Dropbox/Noise_genomics/Tcell-ExonLength_vs_CV2.png",
#        width=4, height=4*1.618, dpi=90)
# 
# p_avexon <- ggplot(tcell.consist, aes(x=log10(EXON_AVLENGTH), y=CV2)) + 
#   geom_point() +
#   theme_mike()
# ggsave(p_avexon, filename="~/Dropbox/Noise_genomics/Tcell-AvExonLength_vs_CV2.png",
#        width=4, height=4*1.618, dpi=90)
# 
# p_sumphast <- ggplot(tcell.consist, aes(x=SUM_PHAST, y=CV2, colour=GeneGroup)) + 
#   geom_point(alpha=0.2) + scale_colour_Publication() +
#   theme_mike()
# ggsave(p_sumphast, filename="~/Dropbox/Noise_genomics/Tcell-SumPhastCons_vs_CV2.png",
#        height=4, width=4*1.618, dpi=90)
# 
# p_align <- ggplot(tcell.consist, aes(x=NALIGN_PHAST, y=CV2)) + 
#   geom_point() +
#   theme_mike()
# ggsave(p_align, filename="~/Dropbox/Noise_genomics/Tcell-AlignedBases_vs_CV2.png",
#        height=4, width=4*1.618, dpi=90)
# 
# p_gc <- ggplot(tcell.consist, aes(x=GC, y=CV2)) + 
#   geom_point() +
#   theme_mike()
# ggsave(p_gc, filename="~/Dropbox/Noise_genomics/Tcell-GC_vs_CV2.png",
#        height=4, width=4*1.618, dpi=90)
# 
# p_group <- ggplot(tcell.consist, aes(x=as.factor(Group), y=CV2)) + 
#   geom_boxplot() + geom_jitter(alpha=0.2) + 
#   theme_mike()
# ggsave(p_group, filename="~/Dropbox/Noise_genomics/Tcell-ExprsGroup_vs_CV2.png",
#        height=5, width=5*1.618, dpi=90)
# 
# genomic.vars <- paste(c("Group", "GeneGroup", colnames(genomic.features)[2:dim(genomic.features)[2]]),
#                       collapse=" + ")
# glm.form <- as.formula(paste("CV2", genomic.vars, sep=" ~ "))
# 
# feature.glm <- glm(glm.form, data=tcell.consist, family=gaussian('identity'))
# summary(feature.glm)
# write.table(summary(feature.glm)$coefficients,
#             file="~/Dropbox/Noise_genomics/Tcell-linearFit_CV2.tsv",
#             sep="\t", quote=F)
# 
# # use LASSO to select meaningful features, set alpha=1 to make
# # it a full LASSO rather than elastic net
# lasso <- glmnet(x=as.matrix(tcell.consist[, c(4, 6:(dim(tcell.consist)[2]))]),
#                 y=tcell.consist[, 5],
#                 family="gaussian", alpha=1)
# plot(lasso, xvar='dev', label=TRUE)
# coef(lasso, s=0.2)
# 
# cv.lasso <- cv.glmnet(x=as.matrix(apply(tcell.consist[, c(4, 6:(dim(tcell.consist)[2]))],2,
#                                         as.numeric)),
#                       y=tcell.consist[, 5],
#                       family="gaussian", alpha=1,
#                       type.measure = "mse", nfolds = 10)
# write.table(as.matrix(coef(cv.lasso, s=cv.lasso$lambda.1se)),
#             file="~/Dropbox/Noise_genomics/Tcell-LASSO_CV2.tsv",
#             sep="\t", quote=F)
# plot(cv.lasso) # the only meaningful feature is mean expression group

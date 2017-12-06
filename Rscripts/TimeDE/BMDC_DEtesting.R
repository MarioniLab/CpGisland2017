# BMDC time-series
library(lsa)
library(mgcv)
library(MASS)
library(org.Mm.eg.db)
library(goseq)
library(stringr)
library(scales)
library(ggrepel)
library(cowplot)
library(biomaRt)
library(RColorBrewer)
library(limma)
library(VennDiagram)

source("~/Dropbox/R_sessions/Noise/mESC_chromHMM.R")
source("~/Dropbox/R_sessions/Noise/genomic_noise_features.R")
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")
source("~/Dropbox/R_sessions/SingleCellFunctions/single_cell_functions.R")

ensembl <- useEnsembl(biomart='ensembl', dataset='mmusculus_gene_ensembl', GRCh=37)

bmdc.cells <- read.table("~/Dropbox/Dendritic_cells/All_BMDC/dendritic_SFnorm.tsv", h=T, 
                         sep="\t", stringsAsFactors=F)
rownames(bmdc.cells) <- bmdc.cells$gene_id

bmdc.meta <- read.table("~/Dropbox/Dendritic_cells/All_BMDC/dendritic-metadata.tsv", h=T,
                        sep="\t", stringsAsFactors=F)

gene_symbol <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name'),
                     filters='ensembl_gene_id', mart=ensembl,
                     values=unique(bmdc.cells$gene_id))

############################################
## DE testing between pairs of timepoints ##
############################################
# use limma for DE testing
# run each stimulation simulataneously

lps.cells <- bmdc.meta$Sample[bmdc.meta$Stim %in% c("Unstim", "LPS")]
pam.cells <- bmdc.meta$Sample[bmdc.meta$Stim %in% c("Unstim", "PAM")]
pic.cells <- bmdc.meta$Sample[bmdc.meta$Stim %in% c("Unstim", "PIC")]

bmdc.meta$Timepoint <- factor(bmdc.meta$Time,
                              labels=c("0h", "1h", "2h", "4h", "6h"),
                              levels=c("0h", "1h", "2h", "4h", "6h"))

bmdc.lps.design <- model.matrix(~ 0 + Time, data=bmdc.meta[bmdc.meta$Sample %in% lps.cells, ])
bmdc.pam.design <- model.matrix(~ 0 + Time, data=bmdc.meta[bmdc.meta$Sample %in% pam.cells, ])
bmdc.pic.design <- model.matrix(~ 0 + Time, data=bmdc.meta[bmdc.meta$Sample %in% pic.cells, ])

bmdc.lps.exprs <- bmdc.cells[, colnames(bmdc.cells) %in% lps.cells]
bmdc.pam.exprs <- bmdc.cells[, colnames(bmdc.cells) %in% pam.cells]
bmdc.pic.exprs <- bmdc.cells[, colnames(bmdc.cells) %in% pic.cells]

###############
### Test LPS ##
###############
bmdc.lps.fit <- lmFit(bmdc.lps.exprs, bmdc.lps.design)
bmdc.lps.constrast <- makeContrasts(contrasts=c("Time1h-Time0h",
                                                       "Time2h-Time1h",
                                                       "Time4h-Time2h",
                                                       "Time6h-Time4h"),
                                           levels=bmdc.lps.design)
bmdc.lps.fit <- contrasts.fit(bmdc.lps.fit, contrasts=bmdc.lps.constrast)
bmdc.lps.fit <- eBayes(bmdc.lps.fit)
bmdc.lps.sum.res <- summary(decideTests(bmdc.lps.fit))

############################
### DE Unstimulated vs 1h ##
############################

bmdc.lps.de.res.0v1 <- topTable(bmdc.lps.fit, coef=1, n=Inf, sort="p", p=1.0)
bmdc.lps.de.res.0v1$Sig <- 0
bmdc.lps.de.res.0v1$Sig[bmdc.lps.de.res.0v1$adj.P.Val <= 0.01] <- 1
bmdc.lps.de.res.0v1$Sig <- as.factor(bmdc.lps.de.res.0v1$Sig)

bmdc.lps.de.res.0v1$Diff <- 0
bmdc.lps.de.res.0v1$Diff[bmdc.lps.de.res.0v1$logFC < 0 & bmdc.lps.de.res.0v1$Sig == 1] <- -1
bmdc.lps.de.res.0v1$Diff[bmdc.lps.de.res.0v1$logFC > 0 & bmdc.lps.de.res.0v1$Sig == 1] <- 1
bmdc.lps.de.res.0v1$Diff <- as.factor(bmdc.lps.de.res.0v1$Diff)

bmdc.lps.de.res.0v1$gene_id <- rownames(bmdc.lps.de.res.0v1)

##################
### DE 1h vs 2h ##
##################

bmdc.lps.de.res.1v2 <- topTable(bmdc.lps.fit, coef=2, n=Inf, sort="p", p=1.0)
bmdc.lps.de.res.1v2$Sig <- 0
bmdc.lps.de.res.1v2$Sig[bmdc.lps.de.res.1v2$adj.P.Val <= 0.01] <- 1
bmdc.lps.de.res.1v2$Sig <- as.factor(bmdc.lps.de.res.1v2$Sig)

bmdc.lps.de.res.1v2$Diff <- 0
bmdc.lps.de.res.1v2$Diff[bmdc.lps.de.res.1v2$logFC < 0 & bmdc.lps.de.res.1v2$Sig == 1] <- -1
bmdc.lps.de.res.1v2$Diff[bmdc.lps.de.res.1v2$logFC > 0 & bmdc.lps.de.res.1v2$Sig == 1] <- 1
bmdc.lps.de.res.1v2$Diff <- as.factor(bmdc.lps.de.res.1v2$Diff)

bmdc.lps.de.res.1v2$gene_id <- rownames(bmdc.lps.de.res.1v2)

##################
### DE 2h vs 4h ##
##################

bmdc.lps.de.res.2v4 <- topTable(bmdc.lps.fit, coef=3, n=Inf, sort="p", p=1.0)
bmdc.lps.de.res.2v4$Sig <- 0
bmdc.lps.de.res.2v4$Sig[bmdc.lps.de.res.2v4$adj.P.Val <= 0.01] <- 1
bmdc.lps.de.res.2v4$Sig <- as.factor(bmdc.lps.de.res.2v4$Sig)

bmdc.lps.de.res.2v4$Diff <- 0
bmdc.lps.de.res.2v4$Diff[bmdc.lps.de.res.2v4$logFC < 0 & bmdc.lps.de.res.2v4$Sig == 1] <- -1
bmdc.lps.de.res.2v4$Diff[bmdc.lps.de.res.2v4$logFC > 0 & bmdc.lps.de.res.2v4$Sig == 1] <- 1
bmdc.lps.de.res.2v4$Diff <- as.factor(bmdc.lps.de.res.2v4$Diff)

bmdc.lps.de.res.2v4$gene_id <- rownames(bmdc.lps.de.res.2v4)

##################
### DE 4h vs 6h ##
##################

bmdc.lps.de.res.4v6 <- topTable(bmdc.lps.fit, coef=4, n=Inf, sort="p", p=1.0)
bmdc.lps.de.res.4v6$Sig <- 0
bmdc.lps.de.res.4v6$Sig[bmdc.lps.de.res.4v6$adj.P.Val <= 0.01] <- 1
bmdc.lps.de.res.4v6$Sig <- as.factor(bmdc.lps.de.res.4v6$Sig)

bmdc.lps.de.res.4v6$Diff <- 0
bmdc.lps.de.res.4v6$Diff[bmdc.lps.de.res.4v6$logFC < 0 & bmdc.lps.de.res.4v6$Sig == 1] <- -1
bmdc.lps.de.res.4v6$Diff[bmdc.lps.de.res.4v6$logFC > 0 & bmdc.lps.de.res.4v6$Sig == 1] <- 1
bmdc.lps.de.res.4v6$Diff <- as.factor(bmdc.lps.de.res.4v6$Diff)

bmdc.lps.de.res.4v6$gene_id <- rownames(bmdc.lps.de.res.4v6)

###############
### Test PAM ##
###############
bmdc.pam.fit <- lmFit(bmdc.pam.exprs, bmdc.pam.design)
bmdc.pam.constrast <- makeContrasts(contrasts=c("Time1h-Time0h",
                                                "Time2h-Time1h",
                                                "Time4h-Time2h",
                                                "Time6h-Time4h"),
                                    levels=bmdc.pam.design)
bmdc.pam.fit <- contrasts.fit(bmdc.pam.fit, contrasts=bmdc.pam.constrast)
bmdc.pam.fit <- eBayes(bmdc.pam.fit)
bmdc.pam.sum.res <- summary(decideTests(bmdc.pam.fit))

############################
### DE Unstimulated vs 1h ##
############################

bmdc.pam.de.res.0v1 <- topTable(bmdc.pam.fit, coef=1, n=Inf, sort="p", p=1.0)
bmdc.pam.de.res.0v1$Sig <- 0
bmdc.pam.de.res.0v1$Sig[bmdc.pam.de.res.0v1$adj.P.Val <= 0.01] <- 1
bmdc.pam.de.res.0v1$Sig <- as.factor(bmdc.pam.de.res.0v1$Sig)

bmdc.pam.de.res.0v1$Diff <- 0
bmdc.pam.de.res.0v1$Diff[bmdc.pam.de.res.0v1$logFC < 0 & bmdc.pam.de.res.0v1$Sig == 1] <- -1
bmdc.pam.de.res.0v1$Diff[bmdc.pam.de.res.0v1$logFC > 0 & bmdc.pam.de.res.0v1$Sig == 1] <- 1
bmdc.pam.de.res.0v1$Diff <- as.factor(bmdc.pam.de.res.0v1$Diff)

bmdc.pam.de.res.0v1$gene_id <- rownames(bmdc.pam.de.res.0v1)

##################
### DE 1h vs 2h ##
##################

bmdc.pam.de.res.1v2 <- topTable(bmdc.pam.fit, coef=2, n=Inf, sort="p", p=1.0)
bmdc.pam.de.res.1v2$Sig <- 0
bmdc.pam.de.res.1v2$Sig[bmdc.pam.de.res.1v2$adj.P.Val <= 0.01] <- 1
bmdc.pam.de.res.1v2$Sig <- as.factor(bmdc.pam.de.res.1v2$Sig)

bmdc.pam.de.res.1v2$Diff <- 0
bmdc.pam.de.res.1v2$Diff[bmdc.pam.de.res.1v2$logFC < 0 & bmdc.pam.de.res.1v2$Sig == 1] <- -1
bmdc.pam.de.res.1v2$Diff[bmdc.pam.de.res.1v2$logFC > 0 & bmdc.pam.de.res.1v2$Sig == 1] <- 1
bmdc.pam.de.res.1v2$Diff <- as.factor(bmdc.pam.de.res.1v2$Diff)

bmdc.pam.de.res.1v2$gene_id <- rownames(bmdc.pam.de.res.1v2)

##################
### DE 2h vs 4h ##
##################

bmdc.pam.de.res.2v4 <- topTable(bmdc.pam.fit, coef=3, n=Inf, sort="p", p=1.0)
bmdc.pam.de.res.2v4$Sig <- 0
bmdc.pam.de.res.2v4$Sig[bmdc.pam.de.res.2v4$adj.P.Val <= 0.01] <- 1
bmdc.pam.de.res.2v4$Sig <- as.factor(bmdc.pam.de.res.2v4$Sig)

bmdc.pam.de.res.2v4$Diff <- 0
bmdc.pam.de.res.2v4$Diff[bmdc.pam.de.res.2v4$logFC < 0 & bmdc.pam.de.res.2v4$Sig == 1] <- -1
bmdc.pam.de.res.2v4$Diff[bmdc.pam.de.res.2v4$logFC > 0 & bmdc.pam.de.res.2v4$Sig == 1] <- 1
bmdc.pam.de.res.2v4$Diff <- as.factor(bmdc.pam.de.res.2v4$Diff)

bmdc.pam.de.res.2v4$gene_id <- rownames(bmdc.pam.de.res.2v4)

##################
### DE 4h vs 6h ##
##################

bmdc.pam.de.res.4v6 <- limma::topTable(bmdc.pam.fit, coef=4, n=Inf, sort="p", p=1.0)
bmdc.pam.de.res.4v6$Sig <- 0
bmdc.pam.de.res.4v6$Sig[bmdc.pam.de.res.4v6$adj.P.Val <= 0.01] <- 1
bmdc.pam.de.res.4v6$Sig <- as.factor(bmdc.pam.de.res.4v6$Sig)

bmdc.pam.de.res.4v6$Diff <- 0
bmdc.pam.de.res.4v6$Diff[bmdc.pam.de.res.4v6$logFC < 0 & bmdc.pam.de.res.4v6$Sig == 1] <- -1
bmdc.pam.de.res.4v6$Diff[bmdc.pam.de.res.4v6$logFC > 0 & bmdc.pam.de.res.4v6$Sig == 1] <- 1
bmdc.pam.de.res.4v6$Diff <- as.factor(bmdc.pam.de.res.4v6$Diff)

bmdc.pam.de.res.4v6$gene_id <- rownames(bmdc.pam.de.res.4v6)

###############
### Test PIC ##
###############

bmdc.pic.fit <- lmFit(bmdc.pic.exprs, bmdc.pic.design)
bmdc.pic.constrast <- makeContrasts(contrasts=c("Time1h-Time0h",
                                                "Time2h-Time1h",
                                                "Time4h-Time2h",
                                                "Time6h-Time4h"),
                                    levels=bmdc.pic.design)
bmdc.pic.fit <- contrasts.fit(bmdc.pic.fit, contrasts=bmdc.pic.constrast)
bmdc.pic.fit <- eBayes(bmdc.pic.fit)
bmdc.pic.sum.res <- summary(decideTests(bmdc.pic.fit))

############################
### DE Unstimulated vs 1h ##
############################

bmdc.pic.de.res.0v1 <- topTable(bmdc.pic.fit, coef=1, n=Inf, sort="p", p=1.0)
bmdc.pic.de.res.0v1$Sig <- 0
bmdc.pic.de.res.0v1$Sig[bmdc.pic.de.res.0v1$adj.P.Val <= 0.01] <- 1
bmdc.pic.de.res.0v1$Sig <- as.factor(bmdc.pic.de.res.0v1$Sig)

bmdc.pic.de.res.0v1$Diff <- 0
bmdc.pic.de.res.0v1$Diff[bmdc.pic.de.res.0v1$logFC < 0 & bmdc.pic.de.res.0v1$Sig == 1] <- -1
bmdc.pic.de.res.0v1$Diff[bmdc.pic.de.res.0v1$logFC > 0 & bmdc.pic.de.res.0v1$Sig == 1] <- 1
bmdc.pic.de.res.0v1$Diff <- as.factor(bmdc.pic.de.res.0v1$Diff)

bmdc.pic.de.res.0v1$gene_id <- rownames(bmdc.pic.de.res.0v1)

##################
### DE 1h vs 2h ##
##################

bmdc.pic.de.res.1v2 <- topTable(bmdc.pic.fit, coef=2, n=Inf, sort="p", p=1.0)
bmdc.pic.de.res.1v2$Sig <- 0
bmdc.pic.de.res.1v2$Sig[bmdc.pic.de.res.1v2$adj.P.Val <= 0.01] <- 1
bmdc.pic.de.res.1v2$Sig <- as.factor(bmdc.pic.de.res.1v2$Sig)

bmdc.pic.de.res.1v2$Diff <- 0
bmdc.pic.de.res.1v2$Diff[bmdc.pic.de.res.1v2$logFC < 0 & bmdc.pic.de.res.1v2$Sig == 1] <- -1
bmdc.pic.de.res.1v2$Diff[bmdc.pic.de.res.1v2$logFC > 0 & bmdc.pic.de.res.1v2$Sig == 1] <- 1
bmdc.pic.de.res.1v2$Diff <- as.factor(bmdc.pic.de.res.1v2$Diff)

bmdc.pic.de.res.1v2$gene_id <- rownames(bmdc.pic.de.res.1v2)

##################
### DE 2h vs 4h ##
##################

bmdc.pic.de.res.2v4 <- topTable(bmdc.pic.fit, coef=3, n=Inf, sort="p", p=1.0)
bmdc.pic.de.res.2v4$Sig <- 0
bmdc.pic.de.res.2v4$Sig[bmdc.pic.de.res.2v4$adj.P.Val <= 0.01] <- 1
bmdc.pic.de.res.2v4$Sig <- as.factor(bmdc.pic.de.res.2v4$Sig)

bmdc.pic.de.res.2v4$Diff <- 0
bmdc.pic.de.res.2v4$Diff[bmdc.pic.de.res.2v4$logFC < 0 & bmdc.pic.de.res.2v4$Sig == 1] <- -1
bmdc.pic.de.res.2v4$Diff[bmdc.pic.de.res.2v4$logFC > 0 & bmdc.pic.de.res.2v4$Sig == 1] <- 1
bmdc.pic.de.res.2v4$Diff <- as.factor(bmdc.pic.de.res.2v4$Diff)

bmdc.pic.de.res.2v4$gene_id <- rownames(bmdc.pic.de.res.2v4)

##################
### DE 4h vs 6h ##
##################

bmdc.pic.de.res.4v6 <- limma::topTable(bmdc.pic.fit, coef=4, n=Inf, sort="p", p=1.0)
bmdc.pic.de.res.4v6$Sig <- 0
bmdc.pic.de.res.4v6$Sig[bmdc.pic.de.res.4v6$adj.P.Val <= 0.01] <- 1
bmdc.pic.de.res.4v6$Sig <- as.factor(bmdc.pic.de.res.4v6$Sig)

bmdc.pic.de.res.4v6$Diff <- 0
bmdc.pic.de.res.4v6$Diff[bmdc.pic.de.res.4v6$logFC < 0 & bmdc.pic.de.res.4v6$Sig == 1] <- -1
bmdc.pic.de.res.4v6$Diff[bmdc.pic.de.res.4v6$logFC > 0 & bmdc.pic.de.res.4v6$Sig == 1] <- 1
bmdc.pic.de.res.4v6$Diff <- as.factor(bmdc.pic.de.res.4v6$Diff)

bmdc.pic.de.res.4v6$gene_id <- rownames(bmdc.pic.de.res.4v6)

################################################################################################
################################################################################################
###################################################################################
## Merge DE results for each comparison and stimulation with CpG island features ##
###################################################################################
## LPS ##
#########

colnames(bmdc.lps.de.res.0v1) <- c(paste("t0_t1", colnames(bmdc.lps.de.res.0v1)[1:8],
                                         sep="."), "gene_id")

colnames(bmdc.lps.de.res.1v2) <- c(paste("t1_t2", colnames(bmdc.lps.de.res.1v2)[1:8],
                                         sep="."), "gene_id")

colnames(bmdc.lps.de.res.2v4) <- c(paste("t2_t4", colnames(bmdc.lps.de.res.2v4)[1:8],
                                         sep="."), "gene_id")

colnames(bmdc.lps.de.res.4v6) <- c(paste("t4_t6", colnames(bmdc.lps.de.res.4v6)[1:8],
                                         sep="."), "gene_id")


bmdc.lps.de_list <- list("t0_t1"=bmdc.lps.de.res.0v1,
                         "t1_t2h"=bmdc.lps.de.res.1v2,
                         "t2_t4h"=bmdc.lps.de.res.2v4,
                         "t4_t6h"=bmdc.lps.de.res.4v6)
bmdc.lps_de.merge <- Reduce(x=bmdc.lps.de_list,
                            f=function(x, y) merge(x, y, by=c('gene_id')))

bmdc.lps.genomic <- merge(genomic.features, bmdc.lps_de.merge, by.x='GENE', by.y='gene_id')
bmdc.lps.genomic$CGI_SIZE.kb <- bmdc.lps.genomic$CGI_SIZE/1000

## PAM ##
#########
colnames(bmdc.pam.de.res.0v1) <- c(paste("t0_t1", colnames(bmdc.pam.de.res.0v1)[1:8],
                                         sep="."), "gene_id")

colnames(bmdc.pam.de.res.1v2) <- c(paste("t1_t2", colnames(bmdc.pam.de.res.1v2)[1:8],
                                         sep="."), "gene_id")

colnames(bmdc.pam.de.res.2v4) <- c(paste("t2_t4", colnames(bmdc.pam.de.res.2v4)[1:8],
                                         sep="."), "gene_id")

colnames(bmdc.pam.de.res.4v6) <- c(paste("t4_t6", colnames(bmdc.pam.de.res.4v6)[1:8],
                                         sep="."), "gene_id")

bmdc.pam.de_list <- list("t0_t1"=bmdc.pam.de.res.0v1,
                         "t1_t2h"=bmdc.pam.de.res.1v2,
                         "t2_t4h"=bmdc.pam.de.res.2v4,
                         "t4_t6h"=bmdc.pam.de.res.4v6)
bmdc.pam_de.merge <- Reduce(x=bmdc.pam.de_list,
                            f=function(x, y) merge(x, y, by=c('gene_id')))

bmdc.pam.genomic <- merge(genomic.features, bmdc.pam_de.merge, by.x='GENE', by.y='gene_id')
bmdc.pam.genomic$CGI_SIZE.kb <- bmdc.pam.genomic$CGI_SIZE/1000

## PIC ##
#########
colnames(bmdc.pic.de.res.0v1) <- c(paste("t0_t1", colnames(bmdc.pic.de.res.0v1)[1:8],
                                                      sep="."), "gene_id")

colnames(bmdc.pic.de.res.1v2) <- c(paste("t1_t2", colnames(bmdc.pic.de.res.1v2)[1:8],
                                                      sep="."), "gene_id")

colnames(bmdc.pic.de.res.2v4) <- c(paste("t2_t4", colnames(bmdc.pic.de.res.2v4)[1:8],
                                                      sep="."), "gene_id")

colnames(bmdc.pic.de.res.4v6) <- c(paste("t4_t6", colnames(bmdc.pic.de.res.4v6)[1:8],
                                                      sep="."), "gene_id")


bmdc.pic.de_list <- list("t0_t1"=bmdc.pic.de.res.0v1,
                         "t1_t2h"=bmdc.pic.de.res.1v2,
                         "t2_t4h"=bmdc.pic.de.res.2v4,
                         "t4_t6h"=bmdc.pic.de.res.4v6)
bmdc.pic_de.merge <- Reduce(x=bmdc.pic.de_list,
                            f=function(x, y) merge(x, y, by=c('gene_id')))

bmdc.pic.genomic <- merge(genomic.features, bmdc.pic_de.merge, by.x='GENE', by.y='gene_id')
bmdc.pic.genomic$CGI_SIZE.kb <- bmdc.pic.genomic$CGI_SIZE/1000

# within each condition, compare the timepoint specific CGI size distributions of DE genes
# get CGI sizes for up-regulated genes at each time point
bmdc.lps.up.0v1 <- bmdc.lps.genomic$CGI_SIZE.kb[bmdc.lps.genomic$t0_t1.Diff == 1 & bmdc.lps.genomic$CGI_SIZE.kb > 0]
bmdc.lps.up.1v2 <- bmdc.lps.genomic$CGI_SIZE.kb[bmdc.lps.genomic$t1_t2.Diff == 1 & bmdc.lps.genomic$CGI_SIZE.kb > 0]
bmdc.lps.up.2v4 <- bmdc.lps.genomic$CGI_SIZE.kb[bmdc.lps.genomic$t2_t4.Diff == 1 & bmdc.lps.genomic$CGI_SIZE.kb > 0]
bmdc.lps.up.4v6 <- bmdc.lps.genomic$CGI_SIZE.kb[bmdc.lps.genomic$t4_t6.Diff == 1 & bmdc.lps.genomic$CGI_SIZE.kb > 0]

# pam
bmdc.pam.up.0v1 <- bmdc.pam.genomic$CGI_SIZE.kb[bmdc.pam.genomic$t0_t1.Diff == 1 & bmdc.pam.genomic$CGI_SIZE.kb > 0]
bmdc.pam.up.1v2 <- bmdc.pam.genomic$CGI_SIZE.kb[bmdc.pam.genomic$t1_t2.Diff == 1 & bmdc.pam.genomic$CGI_SIZE.kb > 0]
bmdc.pam.up.2v4 <- bmdc.pam.genomic$CGI_SIZE.kb[bmdc.pam.genomic$t2_t4.Diff == 1 & bmdc.pam.genomic$CGI_SIZE.kb > 0]
bmdc.pam.up.4v6 <- bmdc.pam.genomic$CGI_SIZE.kb[bmdc.pam.genomic$t4_t6.Diff == 1 & bmdc.pam.genomic$CGI_SIZE.kb > 0]

# pic
bmdc.pic.up.0v1 <- bmdc.pic.genomic$CGI_SIZE.kb[bmdc.pic.genomic$t0_t1.Diff == 1 & bmdc.pic.genomic$CGI_SIZE.kb > 0]
bmdc.pic.up.1v2 <- bmdc.pic.genomic$CGI_SIZE.kb[bmdc.pic.genomic$t1_t2.Diff == 1 & bmdc.pic.genomic$CGI_SIZE.kb > 0]
bmdc.pic.up.2v4 <- bmdc.pic.genomic$CGI_SIZE.kb[bmdc.pic.genomic$t2_t4.Diff == 1 & bmdc.pic.genomic$CGI_SIZE.kb > 0]
bmdc.pic.up.4v6 <- bmdc.pic.genomic$CGI_SIZE.kb[bmdc.pic.genomic$t4_t6.Diff == 1 & bmdc.pic.genomic$CGI_SIZE.kb > 0]

bmdc.lps.size_list <- list("0v1"=cbind(bmdc.lps.up.0v1, rep("0v1", length(bmdc.lps.up.0v1))),
                           "1v2"=cbind(bmdc.lps.up.1v2, rep("1v2", length(bmdc.lps.up.1v2))),
                           "2v4"=cbind(bmdc.lps.up.2v4, rep("2v4", length(bmdc.lps.up.2v4))),
                           "4v6"=cbind(bmdc.lps.up.4v6, rep("4v6", length(bmdc.lps.up.4v6))))

bmdc.lps.size.df <- data.frame(do.call(rbind, bmdc.lps.size_list))
colnames(bmdc.lps.size.df) <- c("CGI_SIZE.kb", "Comparison")
bmdc.lps.size.df$CGI_SIZE.kb <- as.numeric(as.character(bmdc.lps.size.df$CGI_SIZE.kb))
bmdc.lps.size.df$Comparison <- factor(bmdc.lps.size.df$Comparison,
                                      levels=c("0v1", "1v2", "2v4", "4v6"),
                                      labels=c("0v1", "1v2", "2v4", "4v6"))


bmdc.pam.size_list <- list("0v1"=cbind(bmdc.pam.up.0v1, rep("0v1", length(bmdc.pam.up.0v1))),
                           "1v2"=cbind(bmdc.pam.up.1v2, rep("1v2", length(bmdc.pam.up.1v2))),
                           "2v4"=cbind(bmdc.pam.up.2v4, rep("2v4", length(bmdc.pam.up.2v4))),
                           "4v6"=cbind(bmdc.pam.up.4v6, rep("4v6", length(bmdc.pam.up.4v6))))

bmdc.pam.size.df <- data.frame(do.call(rbind, bmdc.pam.size_list))
colnames(bmdc.pam.size.df) <- c("CGI_SIZE.kb", "Comparison")
bmdc.pam.size.df$CGI_SIZE.kb <- as.numeric(as.character(bmdc.pam.size.df$CGI_SIZE.kb))
bmdc.pam.size.df$Comparison <- factor(bmdc.pam.size.df$Comparison,
                                      levels=c("0v1", "1v2", "2v4", "4v6"),
                                      labels=c("0v1", "1v2", "2v4", "4v6"))


bmdc.pic.size_list <- list("0v1"=cbind(bmdc.pic.up.0v1, rep("0v1", length(bmdc.pic.up.0v1))),
                           "1v2"=cbind(bmdc.pic.up.1v2, rep("1v2", length(bmdc.pic.up.1v2))),
                           "2v4"=cbind(bmdc.pic.up.2v4, rep("2v4", length(bmdc.pic.up.2v4))),
                           "4v6"=cbind(bmdc.pic.up.4v6, rep("4v6", length(bmdc.pic.up.4v6))))

bmdc.pic.size.df <- data.frame(do.call(rbind, bmdc.pic.size_list))
colnames(bmdc.pic.size.df) <- c("CGI_SIZE.kb", "Comparison")
bmdc.pic.size.df$CGI_SIZE.kb <- as.numeric(as.character(bmdc.pic.size.df$CGI_SIZE.kb))
bmdc.pic.size.df$Comparison <- factor(bmdc.pic.size.df$Comparison,
                                      levels=c("0v1", "1v2", "2v4", "4v6"),
                                      labels=c("0v1", "1v2", "2v4", "4v6"))

#################
### plotting ####
#################
## LPS ##
#########
bmdc.lps.size.box <- ggplot(bmdc.lps.size.df[bmdc.lps.size.df$Comparison %in% c("0v1", "1v2"), ],
                            aes(x=Comparison, y=CGI_SIZE.kb, fill=Comparison)) +
  geom_boxplot() +
  theme_mike() +
  scale_fill_manual(values=c("#D86200", "#D8AA00")) +
  theme(axis.text=element_text(size=16), axis.title.x=element_blank()) +
  guides(fill=FALSE) +
  labs(y="CpG island Size (kb)", x="Time Comparison") +
  scale_y_continuous(limits=c(0, 3), oob=censor)

ggsave(bmdc.lps.size.box,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/BMDC_LPS-CGIsize_boxplot.png",
       height=4.25, width=3.25, dpi=300)

bmdc.lps.size.dens <- ggplot(bmdc.lps.size.df[bmdc.lps.size.df$Comparison %in% c("0v1", "1v2"), ],
                              aes(x=CGI_SIZE.kb, colour=Comparison)) +
  geom_density(size=2, alpha=0.5) + theme_mike() +
  scale_colour_manual(values=c("#D86200", "#D8AA00")) +
  theme(axis.text=element_text(size=16)) +
  guides(colour=FALSE) +
  labs(x="CpG island Size (kb)", y="Density") +
  geom_segment(aes(y=0, yend=1.5,
                   x=median(bmdc.lps.size.df$CGI_SIZE.kb[bmdc.lps.size.df$Comparison == "0v1"]),
                   xend=median(bmdc.lps.size.df$CGI_SIZE.kb[bmdc.lps.size.df$Comparison == "0v1"])),
               linetype="dashed", colour="#D86200", size=2) +
  geom_segment(aes(y=0, yend=1.5,
                   x=median(bmdc.lps.size.df$CGI_SIZE.kb[bmdc.lps.size.df$Comparison == "1v2"]),
                   xend=median(bmdc.lps.size.df$CGI_SIZE.kb[bmdc.lps.size.df$Comparison == "1v2"])),
               linetype="dashed", colour="#D8AA00", size=2) +
  scale_x_continuous(limits=c(0, 3), oob=censor)

ggsave(bmdc.lps.size.dens,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/BMDC_LPS-CGIsize_density.png",
       height=4.25, width=4.75, dpi=300)

## PAM ##
#########

bmdc.pam.size.box <- ggplot(bmdc.pam.size.df[bmdc.pam.size.df$Comparison %in% c("0v1", "1v2"), ],
                            aes(x=Comparison, y=CGI_SIZE.kb, fill=Comparison)) +
  geom_boxplot() +
  theme_mike() +
  scale_fill_manual(values=c("#D86200", "#D8AA00")) +
  theme(axis.text=element_text(size=16), axis.title.x=element_blank()) +
  guides(fill=FALSE) +
  labs(y="CpG island Size (kb)", x="Time Comparison") +
  scale_y_continuous(limits=c(0, 3), oob=censor)

ggsave(bmdc.pam.size.box,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/BMDC_PAM-CGIsize_boxplot.png",
       height=4.25, width=3.25, dpi=300)

bmdc.pam.size.dens <- ggplot(bmdc.pam.size.df[bmdc.pam.size.df$Comparison %in% c("0v1", "1v2"), ],
                             aes(x=CGI_SIZE.kb, colour=Comparison)) +
  geom_density(size=2, alpha=0.5) + theme_mike() +
  scale_colour_manual(values=c("#D86200", "#D8AA00")) +
  theme(axis.text=element_text(size=16)) +
  guides(colour=FALSE) +
  labs(x="CpG island Size (kb)", y="Density") +
  geom_segment(aes(y=0, yend=1.5,
                   x=median(bmdc.pam.size.df$CGI_SIZE.kb[bmdc.pam.size.df$Comparison == "0v1"]),
                   xend=median(bmdc.pam.size.df$CGI_SIZE.kb[bmdc.pam.size.df$Comparison == "0v1"])),
               linetype="dashed", colour="#D86200", size=2) +
  geom_segment(aes(y=0, yend=1.5,
                   x=median(bmdc.lps.size.df$CGI_SIZE.kb[bmdc.pam.size.df$Comparison == "1v2"]),
                   xend=median(bmdc.lps.size.df$CGI_SIZE.kb[bmdc.pam.size.df$Comparison == "1v2"])),
               linetype="dashed", colour="#D8AA00", size=2) +
  scale_x_continuous(limits=c(0, 3), oob=censor)

ggsave(bmdc.pam.size.dens,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/BMDC_PAM-CGIsize_density.png",
       height=4.25, width=4.75, dpi=300)

## PIC ##
#########
bmdc.pic.size.box <- ggplot(bmdc.pic.size.df[bmdc.pic.size.df$Comparison %in% c("0v1", "1v2"), ],
                            aes(x=Comparison, y=CGI_SIZE.kb, fill=Comparison)) +
  geom_boxplot() +
  theme_mike() +
  scale_fill_manual(values=c("#D86200", "#D8AA00")) +
  theme(axis.text=element_text(size=16), axis.title.x=element_blank()) +
  guides(fill=FALSE) +
  labs(y="CpG island Size (kb)", x="Time Comparison") +
  scale_y_continuous(limits=c(0, 3), oob=censor)

ggsave(bmdc.pic.size.box,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/BMDC_PIC-CGIsize_boxplot.png",
       height=4.25, width=3.25, dpi=300)
bmdc.pic.size.dens <- ggplot(bmdc.pic.size.df[bmdc.pic.size.df$Comparison %in% c("0v1", "1v2"), ],
                             aes(x=CGI_SIZE.kb, colour=Comparison)) +
  geom_density(size=2, alpha=0.5) + theme_mike() +
  scale_colour_manual(values=c("#D86200", "#D8AA00")) +
  theme(axis.text=element_text(size=16)) +
  guides(colour=FALSE) +
  labs(x="CpG island Size (kb)", y="Density") +
  geom_segment(aes(y=0, yend=1.5,
                   x=median(bmdc.pic.size.df$CGI_SIZE.kb[bmdc.pic.size.df$Comparison == "0v1"]),
                   xend=median(bmdc.pic.size.df$CGI_SIZE.kb[bmdc.pic.size.df$Comparison == "0v1"])),
               linetype="dashed", colour="#D86200", size=2) +
  geom_segment(aes(y=0, yend=1.5,
                   x=median(bmdc.pic.size.df$CGI_SIZE.kb[bmdc.pic.size.df$Comparison == "1v2"]),
                   xend=median(bmdc.pic.size.df$CGI_SIZE.kb[bmdc.pic.size.df$Comparison == "1v2"])),
               linetype="dashed", colour="#D8AA00", size=2) +
  scale_x_continuous(limits=c(0, 3), oob=censor)

ggsave(bmdc.pic.size.dens,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/BMDC_PIC-CGIsize_density.png",
       height=4.25, width=4.75, dpi=300)

# # how much overlap is there in the CGI up-regulated genes between conditions - 0 vs 1
# bmdc.lps.0v1.genes <- bmdc.lps.genomic$GENE[bmdc.lps.genomic$t0_t1.Diff == 1 & bmdc.lps.genomic$CGI_SIZE.kb > 0]
# bmdc.pam.0v1.genes <- bmdc.pam.genomic$GENE[bmdc.pam.genomic$t0_t1.Diff == 1 & bmdc.pam.genomic$CGI_SIZE.kb > 0]
# bmdc.pic.0v1.genes <- bmdc.pic.genomic$GENE[bmdc.pic.genomic$t0_t1.Diff == 1 & bmdc.pic.genomic$CGI_SIZE.kb > 0]
# 
# grid.newpage()
# bmdc.venn.0v1 <- draw.triple.venn(area1=length(bmdc.lps.0v1.genes),
#                                   area2=length(bmdc.pam.0v1.genes),
#                                   area3=length(bmdc.pic.0v1.genes),
#                                   n12=length(intersect(bmdc.lps.0v1.genes, bmdc.pam.0v1.genes)),
#                                   n13=length(intersect(bmdc.lps.0v1.genes, bmdc.pic.0v1.genes)),
#                                   n23=length(intersect(bmdc.pam.0v1.genes, bmdc.pic.0v1.genes)),
#                                   n123=length(intersect(bmdc.lps.0v1.genes,
#                                                         intersect(bmdc.pam.0v1.genes, bmdc.pic.0v1.genes))),
#                                   category=c("LPS", "PAM", "PIC"))
# grid.draw(bmdc.venn.0v1)
# 
# # how much overlap is there in the CGI up-regulated genes between conditions - 1 vs 2
# bmdc.lps.1v2.genes <- bmdc.lps.genomic$GENE[bmdc.lps.genomic$t1_t2.Diff == 1 & bmdc.lps.genomic$CGI_SIZE.kb > 0]
# bmdc.pam.1v2.genes <- bmdc.pam.genomic$GENE[bmdc.pam.genomic$t1_t2.Diff == 1 & bmdc.pam.genomic$CGI_SIZE.kb > 0]
# bmdc.pic.1v2.genes <- bmdc.pic.genomic$GENE[bmdc.pic.genomic$t1_t2.Diff == 1 & bmdc.pic.genomic$CGI_SIZE.kb > 0]
# 
# grid.newpage()
# bmdc.venn.1v2 <- draw.triple.venn(area1=length(bmdc.lps.1v2.genes),
#                                   area2=length(bmdc.pam.1v2.genes),
#                                   area3=length(bmdc.pic.1v2.genes),
#                                   n12=length(intersect(bmdc.lps.1v2.genes, bmdc.pam.1v2.genes)),
#                                   n13=length(intersect(bmdc.lps.1v2.genes, bmdc.pic.1v2.genes)),
#                                   n23=length(intersect(bmdc.pam.1v2.genes, bmdc.pic.1v2.genes)),
#                                   n123=length(intersect(bmdc.lps.1v2.genes,
#                                                         intersect(bmdc.pam.1v2.genes, bmdc.pic.1v2.genes))),
#                                   category=c("LPS", "PAM", "PIC"))
# grid.draw(bmdc.venn.1v2)




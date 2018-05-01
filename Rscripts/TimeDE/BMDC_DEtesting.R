# BMDC time-series
library(statmod)
library(lsa)
library(mgcv)
library(MASS)
library(robustbase)
library(stringr)
library(scales)
library(ggrepel)
library(cowplot)
library(biomaRt)
library(RColorBrewer)
library(limma)
library(VennDiagram)
library(mclust)

source("~/Dropbox/R_sessions/Noise/mESC_chromHMM.R")
source("~/Dropbox/R_sessions/Noise/genomic_noise_features.R")
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")
source("~/Dropbox/R_sessions/SingleCellFunctions/single_cell_functions.R")
mouse.genomic.features$CGI_SIZE.kb <- mouse.genomic.features$CGI_SIZE/1000
genomic.features$CGI_SIZE.kb <- genomic.features$CGI_SIZE/1000

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
                                                "Time6h-Time4h",
                                                "Time2h-Time0h"),
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


##################
### DE 0h vs 2h ##
##################

bmdc.lps.de.res.0v2 <- limma::topTable(bmdc.lps.fit, coef=5, n=Inf, sort="p", p=1.0)
bmdc.lps.de.res.0v2$Sig <- 0
bmdc.lps.de.res.0v2$Sig[bmdc.lps.de.res.0v2$adj.P.Val <= 0.01] <- 1
bmdc.lps.de.res.0v2$Sig <- as.factor(bmdc.lps.de.res.0v2$Sig)

bmdc.lps.de.res.0v2$Diff <- 0
bmdc.lps.de.res.0v2$Diff[bmdc.lps.de.res.0v2$logFC < 0 & bmdc.lps.de.res.0v2$Sig == 1] <- -1
bmdc.lps.de.res.0v2$Diff[bmdc.lps.de.res.0v2$logFC > 0 & bmdc.lps.de.res.0v2$Sig == 1] <- 1
bmdc.lps.de.res.0v2$Diff <- as.factor(bmdc.lps.de.res.0v2$Diff)

bmdc.lps.de.res.0v2$gene_id <- rownames(bmdc.lps.de.res.0v2)


###############
### Test PAM ##
###############
bmdc.pam.fit <- lmFit(bmdc.pam.exprs, bmdc.pam.design)
bmdc.pam.constrast <- makeContrasts(contrasts=c("Time1h-Time0h",
                                                "Time2h-Time1h",
                                                "Time4h-Time2h",
                                                "Time6h-Time4h",
                                                "Time2h-Time0h"),
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


##################
### DE 0h vs 2h ##
##################

bmdc.pam.de.res.0v2 <- limma::topTable(bmdc.pam.fit, coef=5, n=Inf, sort="p", p=1.0)
bmdc.pam.de.res.0v2$Sig <- 0
bmdc.pam.de.res.0v2$Sig[bmdc.pam.de.res.0v2$adj.P.Val <= 0.01] <- 1
bmdc.pam.de.res.0v2$Sig <- as.factor(bmdc.pam.de.res.0v2$Sig)

bmdc.pam.de.res.0v2$Diff <- 0
bmdc.pam.de.res.0v2$Diff[bmdc.pam.de.res.0v2$logFC < 0 & bmdc.pam.de.res.0v2$Sig == 1] <- -1
bmdc.pam.de.res.0v2$Diff[bmdc.pam.de.res.0v2$logFC > 0 & bmdc.pam.de.res.0v2$Sig == 1] <- 1
bmdc.pam.de.res.0v2$Diff <- as.factor(bmdc.pam.de.res.0v2$Diff)

bmdc.pam.de.res.0v2$gene_id <- rownames(bmdc.pam.de.res.0v2)

###############
### Test PIC ##
###############

bmdc.pic.fit <- lmFit(bmdc.pic.exprs, bmdc.pic.design)
bmdc.pic.constrast <- makeContrasts(contrasts=c("Time1h-Time0h",
                                                "Time2h-Time1h",
                                                "Time4h-Time2h",
                                                "Time6h-Time4h",
                                                "Time2h-Time0h"),
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

##################
### DE 0h vs 2h ##
##################

bmdc.pic.de.res.0v2 <- limma::topTable(bmdc.pic.fit, coef=5, n=Inf, sort="p", p=1.0)
bmdc.pic.de.res.0v2$Sig <- 0
bmdc.pic.de.res.0v2$Sig[bmdc.pic.de.res.0v2$adj.P.Val <= 0.01] <- 1
bmdc.pic.de.res.0v2$Sig <- as.factor(bmdc.pic.de.res.0v2$Sig)

bmdc.pic.de.res.0v2$Diff <- 0
bmdc.pic.de.res.0v2$Diff[bmdc.pic.de.res.0v2$logFC < 0 & bmdc.pic.de.res.0v2$Sig == 1] <- -1
bmdc.pic.de.res.0v2$Diff[bmdc.pic.de.res.0v2$logFC > 0 & bmdc.pic.de.res.0v2$Sig == 1] <- 1
bmdc.pic.de.res.0v2$Diff <- as.factor(bmdc.pic.de.res.0v2$Diff)

bmdc.pic.de.res.0v2$gene_id <- rownames(bmdc.pic.de.res.0v2)

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

colnames(bmdc.lps.de.res.0v2) <- c(paste("t0_t2", colnames(bmdc.lps.de.res.0v2)[1:8],
                                         sep="."), "gene_id")


bmdc.lps.de_list <- list("t0_t1"=bmdc.lps.de.res.0v1,
                         "t1_t2h"=bmdc.lps.de.res.1v2,
                         "t2_t4h"=bmdc.lps.de.res.2v4,
                         "t4_t6h"=bmdc.lps.de.res.4v6,
                         "t0_t2h"=bmdc.lps.de.res.0v2)
bmdc.lps_de.merge <- Reduce(x=bmdc.lps.de_list,
                            f=function(x, y) merge(x, y, by=c('gene_id')))

bmdc.lps.genomic <- merge(mouse.genomic.features, bmdc.lps_de.merge, by.x='GENE', by.y='gene_id')
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

colnames(bmdc.pam.de.res.0v2) <- c(paste("t0_t2", colnames(bmdc.pam.de.res.0v2)[1:8],
                                         sep="."), "gene_id")


bmdc.pam.de_list <- list("t0_t1"=bmdc.pam.de.res.0v1,
                         "t1_t2h"=bmdc.pam.de.res.1v2,
                         "t2_t4h"=bmdc.pam.de.res.2v4,
                         "t4_t6h"=bmdc.pam.de.res.4v6,
                         "t0_t2h"=bmdc.pam.de.res.0v2)
bmdc.pam_de.merge <- Reduce(x=bmdc.pam.de_list,
                            f=function(x, y) merge(x, y, by=c('gene_id')))

bmdc.pam.genomic <- merge(mouse.genomic.features, bmdc.pam_de.merge, by.x='GENE', by.y='gene_id')
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

colnames(bmdc.pic.de.res.0v2) <- c(paste("t0_t2", colnames(bmdc.pic.de.res.0v2)[1:8],
                                         sep="."), "gene_id")

bmdc.pic.de_list <- list("t0_t1"=bmdc.pic.de.res.0v1,
                         "t1_t2h"=bmdc.pic.de.res.1v2,
                         "t2_t4h"=bmdc.pic.de.res.2v4,
                         "t4_t6h"=bmdc.pic.de.res.4v6,
                         "t0_t2h"=bmdc.pic.de.res.0v2)
bmdc.pic_de.merge <- Reduce(x=bmdc.pic.de_list,
                            f=function(x, y) merge(x, y, by=c('gene_id')))

bmdc.pic.genomic <- merge(mouse.genomic.features, bmdc.pic_de.merge, by.x='GENE', by.y='gene_id')

################################################################################################
#### use a binomial test to find differences in the CpG island size ranks between timepoints ###
################################################################################################
###########
### LPS ###
###########
bmdc.lps.0v1.merge <- merge(bmdc.lps.de.res.0v1, mouse.genomic.features, by.x='gene_id', by.y='GENE')
bmdc.lps.1v2.merge <- merge(bmdc.lps.de.res.1v2, mouse.genomic.features, by.x='gene_id', by.y='GENE')
bmdc.lps.0v2.merge <- merge(bmdc.lps.de.res.0v2, mouse.genomic.features, by.x='gene_id', by.y='GENE')

# select only the up-regulated genes
bmdc.lps.0v1.merge.cgi <- bmdc.lps.0v1.merge[bmdc.lps.0v1.merge$CGI_SIZE.kb != 0 &
                                               !is.na(bmdc.lps.0v1.merge$CGI_SIZE.kb) &
                                               bmdc.lps.0v1.merge$t0_t1.logFC > 0, ]
bmdc.lps.1v2.merge.cgi <- bmdc.lps.1v2.merge[bmdc.lps.1v2.merge$CGI_SIZE.kb != 0 & 
                                               !is.na(bmdc.lps.1v2.merge$CGI_SIZE.kb) &
                                               bmdc.lps.1v2.merge$t1_t2.logFC > 0, ]
bmdc.lps.0v2.merge.cgi <- bmdc.lps.0v2.merge[bmdc.lps.0v2.merge$CGI_SIZE.kb != 0 & 
                                               !is.na(bmdc.lps.0v2.merge$CGI_SIZE.kb) &
                                               bmdc.lps.0v2.merge$t0_t2.logFC > 0, ]

# order based on t-statistic
bmdc.lps.res01.size_rank <- bmdc.lps.0v1.merge.cgi$CGI_SIZE.kb[order(bmdc.lps.0v1.merge.cgi$t0_t1.t, decreasing=TRUE)]
bmdc.lps.res26.size_rank <- bmdc.lps.1v2.merge.cgi$CGI_SIZE.kb[order(bmdc.lps.1v2.merge.cgi$t1_t2.t, decreasing=TRUE)]
bmdc.lps.res02.size_rank <- bmdc.lps.0v2.merge.cgi$CGI_SIZE.kb[order(bmdc.lps.0v2.merge.cgi$t0_t2.t, decreasing=TRUE)]

bmdc.lps.sign.size_rank <- as.numeric(bmdc.lps.res01.size_rank < 
                                        bmdc.lps.res26.size_rank)
binom.test(sum(bmdc.lps.sign.size_rank, na.rm=TRUE), n=length(bmdc.lps.sign.size_rank),
           alternative="greater")

# test 0v1 Vs 0v2
bmdc.lps.sign.size_rank.0v2 <- as.numeric(bmdc.lps.res01.size_rank < bmdc.lps.res02.size_rank)
binom.test(sum(bmdc.lps.sign.size_rank.0v2, na.rm=TRUE), n=length(bmdc.lps.sign.size_rank.0v2),
           alternative="greater")

# test 1v2 Vs 0v2
bmdc.lps.sign.size_rank.0v2 <- as.numeric(bmdc.lps.res02.size_rank < bmdc.lps.res26.size_rank)
binom.test(sum(bmdc.lps.sign.size_rank.0v2, na.rm=TRUE), n=length(bmdc.lps.sign.size_rank.0v2),
           alternative="greater")


###########
### PIC ###
###########

bmdc.pic.0v1.merge <- merge(bmdc.pic.de.res.0v1, mouse.genomic.features, by.x='gene_id', by.y='GENE')
bmdc.pic.1v2.merge <- merge(bmdc.pic.de.res.1v2, mouse.genomic.features, by.x='gene_id', by.y='GENE')

# select only the up-regulated genes
bmdc.pic.0v1.merge.cgi <- bmdc.pic.0v1.merge[bmdc.pic.0v1.merge$CGI_SIZE.kb != 0 &
                                               !is.na(bmdc.pic.0v1.merge$CGI_SIZE.kb) &
                                               bmdc.pic.0v1.merge$t0_t1.logFC > 0, ]
bmdc.pic.1v2.merge.cgi <- bmdc.pic.1v2.merge[bmdc.pic.1v2.merge$CGI_SIZE.kb != 0 & 
                                               !is.na(bmdc.pic.1v2.merge$CGI_SIZE.kb) &
                                               bmdc.pic.1v2.merge$t1_t2.logFC > 0, ]

# order based on t-statistic
bmdc.pic.res01.size_rank <- bmdc.pic.0v1.merge.cgi$CGI_SIZE.kb[order(bmdc.pic.0v1.merge.cgi$t0_t1.t, decreasing=TRUE)]
bmdc.pic.res26.size_rank <- bmdc.pic.1v2.merge.cgi$CGI_SIZE.kb[order(bmdc.pic.1v2.merge.cgi$t1_t2.t, decreasing=TRUE)]

bmdc.pic.sign.size_rank <- as.numeric(bmdc.pic.res01.size_rank < bmdc.pic.res26.size_rank)
binom.test(sum(bmdc.pic.sign.size_rank), n=length(bmdc.pic.sign.size_rank),
           alternative="greater")


###########
### PAM ###
###########

bmdc.pam.0v1.merge <- merge(bmdc.pam.de.res.0v1, mouse.genomic.features, by.x='gene_id', by.y='GENE')
bmdc.pam.1v2.merge <- merge(bmdc.pam.de.res.1v2, mouse.genomic.features, by.x='gene_id', by.y='GENE')

# select only the up-regulated genes
bmdc.pam.0v1.merge.cgi <- bmdc.pam.0v1.merge[bmdc.pam.0v1.merge$CGI_SIZE.kb != 0 &
                                               !is.na(bmdc.pam.0v1.merge$CGI_SIZE.kb) &
                                               bmdc.pam.0v1.merge$t0_t1.logFC > 0, ]
bmdc.pam.1v2.merge.cgi <- bmdc.pam.1v2.merge[bmdc.pam.1v2.merge$CGI_SIZE.kb != 0 & 
                                               !is.na(bmdc.pam.1v2.merge$CGI_SIZE.kb) &
                                               bmdc.pam.1v2.merge$t1_t2.logFC > 0, ]
# order based on t-statistic
bmdc.pam.res01.size_rank <- bmdc.pam.0v1.merge.cgi$CGI_SIZE.kb[order(bmdc.pam.0v1.merge.cgi$t0_t1.t, decreasing=TRUE)]
bmdc.pam.res26.size_rank <- bmdc.pam.1v2.merge.cgi$CGI_SIZE.kb[order(bmdc.pam.1v2.merge.cgi$t1_t2.t, decreasing=TRUE)]

bmdc.pam.sign.size_rank <- as.numeric(bmdc.pam.res01.size_rank < bmdc.pam.res26.size_rank)
binom.test(sum(bmdc.pam.sign.size_rank), n=length(bmdc.pam.sign.size_rank),
           alternative="greater")


# within each condition, compare the timepoint specific CGI size distributions of DE genes
# get CGI sizes for up-regulated genes at each time point
# top 250 genes from each only
bmdc.lps.up.0v1 <- bmdc.lps.0v1.merge.cgi$CGI_SIZE.kb
bmdc.lps.up.0v1 <- bmdc.lps.up.0v1[order(bmdc.lps.0v1.merge.cgi$t0_t1.t, decreasing=TRUE)][1:250]

bmdc.lps.up.1v2 <- bmdc.lps.1v2.merge.cgi$CGI_SIZE.kb
bmdc.lps.up.1v2 <- bmdc.lps.up.1v2[order(bmdc.lps.1v2.merge.cgi$t1_t2.t, decreasing=TRUE)][1:250]


bmdc.lps.size_list <- list("0v1"=cbind(bmdc.lps.up.0v1, rep("0v1", length(bmdc.lps.up.0v1))),
                           "1v2"=cbind(bmdc.lps.up.1v2, rep("1v2", length(bmdc.lps.up.1v2))))

bmdc.lps.size.df <- data.frame(do.call(rbind, bmdc.lps.size_list))
colnames(bmdc.lps.size.df) <- c("CGI_SIZE.kb", "Comparison")
bmdc.lps.size.df$CGI_SIZE.kb <- as.numeric(as.character(bmdc.lps.size.df$CGI_SIZE.kb))
bmdc.lps.size.df$Comparison <- factor(bmdc.lps.size.df$Comparison,
                                      levels=c("0v1", "1v2"),
                                      labels=c("0v1", "1v2"))

# pam
bmdc.pam.up.0v1 <- bmdc.pic.0v1.merge.cgi$CGI_SIZE.kb
bmdc.pam.up.0v1 <- bmdc.pam.up.0v1[order(bmdc.pic.0v1.merge.cgi$t0_t1.t, decreasing=TRUE)][1:250]

bmdc.pam.up.1v2 <- bmdc.pic.1v2.merge.cgi$CGI_SIZE.kb
bmdc.pam.up.1v2 <- bmdc.pam.up.1v2[order(bmdc.pic.1v2.merge.cgi$t1_t2.t, decreasing=TRUE)][1:250]

bmdc.pam.size_list <- list("0v1"=cbind(bmdc.pam.up.0v1, rep("0v1", length(bmdc.pam.up.0v1))),
                           "1v2"=cbind(bmdc.pam.up.1v2, rep("1v2", length(bmdc.pam.up.1v2))))

bmdc.pam.size.df <- data.frame(do.call(rbind, bmdc.pam.size_list))
colnames(bmdc.pam.size.df) <- c("CGI_SIZE.kb", "Comparison")
bmdc.pam.size.df$CGI_SIZE.kb <- as.numeric(as.character(bmdc.pam.size.df$CGI_SIZE.kb))
bmdc.pam.size.df$Comparison <- factor(bmdc.pam.size.df$Comparison,
                                      levels=c("0v1", "1v2"),
                                      labels=c("0v1", "1v2"))

# pic
bmdc.pic.up.0v1 <- bmdc.pam.0v1.merge.cgi$CGI_SIZE.kb
bmdc.pic.up.0v1 <- bmdc.pic.up.0v1[order(bmdc.pam.0v1.merge.cgi$t0_t1.t, decreasing=TRUE)][1:250]

bmdc.pic.up.1v2 <- bmdc.pam.1v2.merge.cgi$CGI_SIZE.kb
bmdc.pic.up.1v2 <- bmdc.pic.up.1v2[order(bmdc.pam.1v2.merge.cgi$t1_t2.t, decreasing=TRUE)][1:250]

bmdc.pic.size_list <- list("0v1"=cbind(bmdc.pic.up.0v1, rep("0v1", length(bmdc.pic.up.0v1))),
                           "1v2"=cbind(bmdc.pic.up.1v2, rep("1v2", length(bmdc.lps.up.1v2))))

bmdc.pic.size.df <- data.frame(do.call(rbind, bmdc.pic.size_list))
colnames(bmdc.pic.size.df) <- c("CGI_SIZE.kb", "Comparison")
bmdc.pic.size.df$CGI_SIZE.kb <- as.numeric(as.character(bmdc.pic.size.df$CGI_SIZE.kb))
bmdc.pic.size.df$Comparison <- factor(bmdc.pic.size.df$Comparison,
                                      levels=c("0v1", "1v2"),
                                      labels=c("0v1", "1v2"))

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

###########################################################################################################
# we should expect that the noisy genes in the unstimulated cells are also the ones that are up-regulated #
# post-stimulation
# split by stimulation
#########
## LPS ##
#########
# calculate the variability in expression at the 0h

lps.0h.cells <- na.omit(bmdc.meta$Sample[bmdc.meta$Timepoint == "0h"])
lps.0h.mean <- rowMeans(bmdc.lps.exprs[, lps.0h.cells])
lps.0h.var <- apply(bmdc.lps.exprs[, lps.0h.cells], 1, var)
lps.0h.cv2 <- lps.0h.var/(lps.0h.mean**2)

# calculate the residual CV^2
useForFit <- lps.0h.mean < 0.1

# fit with a gamma-distributed GLM
fit <- glmgam.fit(cbind(a0 = 1, a1tilde=1/lps.0h.mean[!useForFit]), 
                  lps.0h.cv2[!useForFit])

# the negative values don't make a whole lot of sense, they should start at 0
lps.0h.rCV2 <- abs( lps.0h.cv2[!useForFit] - fitted.values(fit))

# plot gene expression varibility Vs logFC
lps.0h.sum <- do.call(cbind.data.frame,
                      list("Mean"=lps.0h.mean[!useForFit],
                           "Var"=lps.0h.var[!useForFit],
                           "CV2"=lps.0h.cv2[!useForFit],
                           "rCV2"=lps.0h.rCV2,
                           "gene_id"=rownames(bmdc.lps.exprs)[!useForFit]))

# merge with genomic features
lps.0h.demerge <- merge(lps.0h.sum, bmdc.lps.de.res.0v1, by='gene_id')
lps.1h.demerge <- merge(lps.0h.demerge, bmdc.lps.de.res.1v2, by='gene_id')
lps.0h.genomic <- merge(lps.0h.demerge, mouse.genomic.features, by.x='gene_id', by.y='GENE')
lps.1h.genomic <- merge(lps.1h.demerge, mouse.genomic.features, by.x='gene_id', by.y='GENE')

# split CGIs into size bins on quartiles
lps.0h.genomic$CGI_SIZE.group <- as.character(cut(lps.0h.genomic$CGI_SIZE.kb,
                                                  breaks=c(0.5, 1.0, 1.5, 2.0)))
lps.0h.genomic$CGI_SIZE.group[(lps.0h.genomic$CGI_SIZE.kb <= 0.5)] <- "(0,0.5]"
lps.0h.genomic$CGI_SIZE.group[(lps.0h.genomic$CGI_SIZE.kb > 2)] <- "(2,4.585]"
lps.0h.genomic$CGI_SIZE.group[is.na(lps.0h.genomic$CGI_SIZE.group) | lps.0h.genomic$CGI_SIZE.kb == 0] <- "Absent"
lps.0h.genomic$CGI_SIZE.group <- factor(lps.0h.genomic$CGI_SIZE.group,
                                                labels=c("Absent", "(0,0.5]", "(0.5,1]", "(1,1.5]",  "(1.5,2]",
                                                         "(2,4.585]"),
                                                levels=c("Absent", "(0,0.5]", "(0.5,1]", "(1,1.5]",  "(1.5,2]",
                                                         "(2,4.585]"))

# split CGIs into size bins on quartiles
lps.1h.genomic$CGI_SIZE.group <- as.character(cut(lps.1h.genomic$CGI_SIZE.kb,
                                                  breaks=c(0.5, 1.0, 1.5, 2.0)))
lps.1h.genomic$CGI_SIZE.group[(lps.1h.genomic$CGI_SIZE.kb <= 0.5)] <- "(0,0.5]"
lps.1h.genomic$CGI_SIZE.group[(lps.1h.genomic$CGI_SIZE.kb > 2)] <- "(2,4.585]"
lps.1h.genomic$CGI_SIZE.group[is.na(lps.1h.genomic$CGI_SIZE.group) | lps.1h.genomic$CGI_SIZE.kb == 0] <- "Absent"
lps.1h.genomic$CGI_SIZE.group <- factor(lps.1h.genomic$CGI_SIZE.group,
                                        labels=c("Absent", "(0,0.5]", "(0.5,1]", "(1,1.5]",  "(1.5,2]",
                                                 "(2,4.585]"),
                                        levels=c("Absent", "(0,0.5]", "(0.5,1]", "(1,1.5]",  "(1.5,2]",
                                                 "(2,4.585]"))
# exclude the bad genes
# bmdc.bad_genes <- as.character(lps.0h.genomic$gene_id[lps.0h.genomic$Mclust %in% c(1, 8)])
# lps.0h.genomic <- lps.0h.genomic[!lps.0h.genomic$gene_id %in% bmdc.bad_genes,]

# remove the absent factor level before running rLM
lps.0h.genomic$CGI_SIZE.group <- factor(lps.0h.genomic$CGI_SIZE.group,
                                        labels=c("(0,0.5]", "(0.5,1]", "(1,1.5]",  "(1.5,2]", "(2,4.585]"),
                                        levels=c("(0,0.5]", "(0.5,1]", "(1,1.5]",  "(1.5,2]", "(2,4.585]"))

# remove the absent factor level before running rLM
lps.1h.genomic$CGI_SIZE.group <- factor(lps.1h.genomic$CGI_SIZE.group,
                                        labels=c("(0,0.5]", "(0.5,1]", "(1,1.5]",  "(1.5,2]", "(2,4.585]"),
                                        levels=c("(0,0.5]", "(0.5,1]", "(1,1.5]",  "(1.5,2]", "(2,4.585]"))

model.control <- lmrob.control(max.it=500, k.max=500, rel.tol=1e-7)
rlm.size <- lmrob(rCV2 ~ t0_t1.Diff + CGI_SIZE.group, control=model.control,
                  data=lps.0h.genomic[lps.0h.genomic$N_CpG == 1,])


lps.cv2.by.cgi <- ggplot(lps.0h.genomic[lps.0h.genomic$N_CpG == 1,],
       aes(y=CV2, x=CGI_SIZE.group, fill=t0_t1.Diff,
           colour=t0_t1.Diff)) +
  geom_jitter(alpha=0.5, 
              position=position_jitterdodge(jitter.width=0.5)) +
  stat_summary(fun.y=mean, colour="grey", geom="point", size=4,
               aes(group=t0_t1.Diff), position=position_dodge(width=0.75)) +
  #geom_boxplot(colour="black") + 
  theme_mike() +
  #scale_fill_Publication(values=c("darkblue", "yellow", "darkred")) +
  scale_colour_manual(values=c("darkblue", "yellow", "darkred")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
  labs(x="CpG island size interval (kb)",
       y=expression(paste("Unstimulated CV"^2))) +
  guides(colour=FALSE, fill=FALSE) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16))

ggsave(lps.cv2.by.cgi,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/BMDC_LPS-CGIsize-CV2-scatter.png",
       width=7.75, height=4.25, dpi=300)


lps.rcv2.by.cgi <- ggplot(lps.0h.genomic[lps.0h.genomic$N_CpG == 1,],
                         aes(y=rCV2, x=CGI_SIZE.group, fill=t0_t1.Diff,
                             colour=t0_t1.Diff, group=t0_t1.Diff)) +
  geom_jitter(alpha=0.5, 
              position=position_jitterdodge(jitter.width=0.5)) +
  #geom_boxplot(colour="black") + 
  theme_mike() +
  stat_summary(fun.y=mean, colour="grey", geom="point", size=4,
               aes(group=t0_t1.Diff), position=position_dodge(width=0.75)) + 
  scale_colour_manual(values=c("darkblue", "yellow", "darkred")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
  labs(x="CpG island size interval (kb)",
       y=expression(paste("Unstimulated Residual CV"^2))) +
  guides(colour=FALSE, fill=FALSE) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16))

ggsave(lps.rcv2.by.cgi,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/BMDC_LPS-CGIsize-residualCV2-scatter.png",
       width=7.75, height=4.25, dpi=300)

# show the baseline mean expression
lps.mean.by.cgi <- ggplot(lps.0h.genomic[lps.0h.genomic$N_CpG == 1,],
       aes(y=Mean, x=CGI_SIZE.group, fill=t0_t1.Diff,
           colour=t0_t1.Diff, group=t0_t1.Diff)) +
  geom_jitter(alpha=0.5, 
              position=position_jitterdodge(jitter.width=0.5)) +
  theme_mike() +
  stat_summary(fun.y=mean, colour="grey", geom="point", size=4,
               aes(group=t0_t1.Diff), position=position_dodge(width=0.75)) + 
  scale_colour_manual(values=c("darkblue", "yellow", "darkred")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
  labs(x="CpG island size interval (kb)",
       y=expression(paste("Mean Log"[2], " Normalized  Expression"))) +
  guides(colour=FALSE, fill=FALSE) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16))

ggsave(lps.mean.by.cgi,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/BMDC_LPS-CGIsize-Mean-scatter.png",
       width=7.75, height=5.25, dpi=300)


## what about the non-CpG island genes?  Is the same pattern seen there?
lps.cv2.by.non_cgi <- ggplot(lps.0h.genomic[lps.0h.genomic$N_CpG == 0,],
                          aes(y=rCV2, x=CGI_SIZE.group, fill=t0_t1.Diff,
                              colour=t0_t1.Diff, group=t0_t1.Diff)) +
  #geom_boxplot(colour="black", alpha=0.4) + 
  geom_jitter(alpha=0.5, 
              position=position_jitterdodge(jitter.width=0.5)) +
  theme_mike() +
  stat_summary(fun.y=mean, colour="grey", geom="point", size=4,
               aes(group=t0_t1.Diff), position=position_dodge(width=0.75)) + 
  scale_colour_manual(values=c("darkblue", "yellow", "darkred")) +
  scale_fill_manual(values=c("darkblue", "yellow", "darkred")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
  labs(x="Differential Expression",
       y=expression(paste("Unstimulated Residual CV"^2))) +
  guides(colour=FALSE, fill=FALSE) +
  theme(axis.text=element_text(size=16),
        axis.text.x=element_blank(),
        axis.title=element_text(size=16))

ggsave(lps.cv2.by.non_cgi,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/BMDC_LPS-nonCGI-residualCV2-scatter.png",
       width=7.75, height=5.25, dpi=300)



### plot the t1_t2.Diff as well to see if noisy non-CGI promtors respond later
# ggplot(lps.1h.genomic[lps.1h.genomic$N_CpG == 0,],
#        aes(y=rCV2, x=CGI_SIZE.group, fill=t1_t2.Diff,
#            colour=t1_t2.Diff, group=t1_t2.Diff)) +
#   #geom_boxplot(colour="black", alpha=0.4) + 
#   geom_jitter(alpha=0.5, 
#               position=position_jitterdodge(jitter.width=0.5)) +
#   theme_mike() +
#   stat_summary(fun.y=mean, colour="grey", geom="point", size=4,
#                aes(group=t0_t1.Diff), position=position_dodge(width=0.75)) + 
#   scale_colour_manual(values=c("darkblue", "yellow", "darkred")) +
#   scale_fill_manual(values=c("darkblue", "yellow", "darkred")) +
#   theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
#   labs(x="Differential Expression",
#        y=expression(paste("Unstimulated Residual CV"^2))) +
#   guides(colour=FALSE, fill=FALSE) +
#   theme(axis.text=element_text(size=16),
#         axis.text.x=element_blank(),
#         axis.title=element_text(size=16))




lps.cv2.by.cgi.cdf <- ggplot(lps.0h.genomic[lps.0h.genomic$N_CpG == 1,],
       aes(x=CV2,
           colour=t0_t1.Diff)) +
  stat_ecdf() +
  theme_mike() +
  #scale_fill_Publication(values=c("darkblue", "yellow", "darkred")) +
  scale_colour_manual(values=c("darkblue", "yellow", "darkred")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
  labs(x="CpG island size interval (kb)",
       y=expression(paste("Unstimulated CV"^2))) +
  facet_wrap(~CGI_SIZE.group,
             ncol=3) +
  guides(colour=FALSE, fill=FALSE)

ggsave(lps.cv2.by.cgi.cdf,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/BMDC_LPS-CGIsize-CV2-ecdf.png",
       width=7.75, height=4.25, dpi=300)

raw.x <- lps.0h.genomic$CV2[lps.0h.genomic$N_CpG == 1 & lps.0h.genomic$t0_t1.Diff == 1]
raw.y <- lps.0h.genomic$CV2[lps.0h.genomic$N_CpG == 1 & lps.0h.genomic$t0_t1.Diff == 0]
raw.down <- lps.0h.genomic$CV2[lps.0h.genomic$N_CpG == 1 & lps.0h.genomic$t0_t1.Diff == -1]
ks.test(raw.x, raw.y, alternative="less")
ks.test(raw.x, raw.down, alternative="less")$p.value

## set up a kolmogorov-smirnoff test between the undiff and up-diff for each CpG island size interval
# (0,0.5]
int1.x <- lps.0h.genomic$CV2[lps.0h.genomic$CGI_SIZE.group == "(0,0.5]" & lps.0h.genomic$t0_t1.Diff == 1]
int1.y <- lps.0h.genomic$CV2[lps.0h.genomic$CGI_SIZE.group == "(0,0.5]" & lps.0h.genomic$t0_t1.Diff == 0]
ks.test(int1.x, int1.y, alternative="less")

# (0.5,1]
int2.x <- lps.0h.genomic$CV2[lps.0h.genomic$CGI_SIZE.group == "(0.5,1]" & lps.0h.genomic$t0_t1.Diff == 1]
int2.y <- lps.0h.genomic$CV2[lps.0h.genomic$CGI_SIZE.group == "(0.5,1]" & lps.0h.genomic$t0_t1.Diff == 0]
ks.test(int2.x, int2.y, alternative="less")

# (1,1.5]
int3.x <- lps.0h.genomic$CV2[lps.0h.genomic$CGI_SIZE.group == "(1,1.5]" & lps.0h.genomic$t0_t1.Diff == 1]
int3.y <- lps.0h.genomic$CV2[lps.0h.genomic$CGI_SIZE.group == "(1,1.5]" & lps.0h.genomic$t0_t1.Diff == 0]
ks.test(int3.x, int3.y, alternative="less")

# (1.5,2]
int4.x <- lps.0h.genomic$CV2[lps.0h.genomic$CGI_SIZE.group == "(1.5,2]" & lps.0h.genomic$t0_t1.Diff == 1]
int4.y <- lps.0h.genomic$CV2[lps.0h.genomic$CGI_SIZE.group == "(1.5,2]" & lps.0h.genomic$t0_t1.Diff == 0]
ks.test(int4.x, int4.y, alternative="less")

# (2,2.5]
int5.x <- lps.0h.genomic$CV2[lps.0h.genomic$CGI_SIZE.group == "(2,4.585]" & lps.0h.genomic$t0_t1.Diff == 1]
int5.y <- lps.0h.genomic$CV2[lps.0h.genomic$CGI_SIZE.group == "(2,4.585]" & lps.0h.genomic$t0_t1.Diff == 0]
ks.test(int5.x, int5.y, alternative="less")

# do you see the same effect with the absolute residual CV^2?
## set up a kolmogorov-smirnoff test between the undiff and up-diff for each CpG island size interval
# (0,0.5]
rcv2.int1.x <- lps.0h.genomic$rCV2[lps.0h.genomic$CGI_SIZE.group == "(0,0.5]" & lps.0h.genomic$t0_t1.Diff == 1]
rcv2.int1.y <- lps.0h.genomic$rCV2[lps.0h.genomic$CGI_SIZE.group == "(0,0.5]" & lps.0h.genomic$t0_t1.Diff == 0]
ks.test(rcv2.int1.x, rcv2.int1.y, alternative="less")

# (0.5,1]
rcv2.int2.x <- lps.0h.genomic$rCV2[lps.0h.genomic$CGI_SIZE.group == "(0.5,1]" & lps.0h.genomic$t0_t1.Diff == 1]
rcv2.int2.y <- lps.0h.genomic$rCV2[lps.0h.genomic$CGI_SIZE.group == "(0.5,1]" & lps.0h.genomic$t0_t1.Diff == 0]
ks.test(rcv2.int2.x, rcv2.int2.y, alternative="less")

# (1,1.5]
rcv2.int3.x <- lps.0h.genomic$rCV2[lps.0h.genomic$CGI_SIZE.group == "(1,1.5]" & lps.0h.genomic$t0_t1.Diff == 1]
rcv2.int3.y <- lps.0h.genomic$rCV2[lps.0h.genomic$CGI_SIZE.group == "(1,1.5]" & lps.0h.genomic$t0_t1.Diff == 0]
ks.test(rcv2.int3.x, rcv2.int3.y, alternative="less")

# (1.5,2]
rcv2.int4.x <- lps.0h.genomic$rCV2[lps.0h.genomic$CGI_SIZE.group == "(1.5,2]" & lps.0h.genomic$t0_t1.Diff == 1]
rcv2.int4.y <- lps.0h.genomic$rCV2[lps.0h.genomic$CGI_SIZE.group == "(1.5,2]" & lps.0h.genomic$t0_t1.Diff == 0]
ks.test(rcv2.int4.x, rcv2.int4.y, alternative="less")

# (2,4.585]
rcv2.int5.x <- lps.0h.genomic$rCV2[lps.0h.genomic$CGI_SIZE.group == "(2,4.585]" & lps.0h.genomic$t0_t1.Diff == 1]
rcv2.int5.y <- lps.0h.genomic$rCV2[lps.0h.genomic$CGI_SIZE.group == "(2,4.585]" & lps.0h.genomic$t0_t1.Diff == 0]
ks.test(rcv2.int5.x, rcv2.int5.y, alternative="less")

int_all.x <- lps.0h.genomic$CV2[lps.0h.genomic$t0_t1.Diff == 1]
int_all.y <- lps.0h.genomic$CV2[lps.0h.genomic$t0_t1.Diff == 0]
int_all.d <- lps.0h.genomic$CV2[lps.0h.genomic$t0_t1.Diff == -1]
ks.test(int_all.x, int_all.y, alternative="less")
ks.test(int_all.x, int_all.d, alternative="less")

# fit a single model, peform ANOVA


#########
## PIC ##
#########
# calculate the variability in expression at the 0h

pic.0h.cells <- na.omit(bmdc.meta$Sample[bmdc.meta$Timepoint == "0h"])
pic.0h.mean <- rowMeans(bmdc.pic.exprs[, pic.0h.cells])
pic.0h.var <- apply(bmdc.pic.exprs[, pic.0h.cells], 1, var)
pic.0h.cv2 <- pic.0h.var/(pic.0h.mean**2)

# calculate the residual CV^2
useForFit <- is.na(pic.0h.cv2) | pic.0h.mean < 0.1 | is.na(1/pic.0h.mean)

smoothScatter(x=pic.0h.mean[!useForFit], y=pic.0h.cv2[!useForFit])

# fit with a gamma-distributed GLM
fit <- glmgam.fit(cbind(a0 = 1, a1tilde=1/pic.0h.mean[!useForFit]), 
                  pic.0h.cv2[!useForFit])
pic.0h.rCV2 <- abs(pic.0h.cv2[!useForFit] - fitted.values(fit))

# plot gene expression varibility Vs logFC
pic.0h.sum <- do.call(cbind.data.frame,
                      list("Mean"=pic.0h.mean[!useForFit],
                           "CV2"=pic.0h.cv2[!useForFit],
                           "rCV2"=pic.0h.rCV2,
                           "gene_id"=rownames(bmdc.pic.exprs)[!useForFit]))

# merge with genomic features
pic.0h.demerge <- merge(pic.0h.sum, bmdc.pic.de.res.0v1, by='gene_id')
pic.0h.genomic <- merge(pic.0h.demerge, mouse.genomic.features, by.x='gene_id', by.y='GENE')

# split CGIs into size bins on quartiles
pic.0h.genomic$CGI_SIZE.group <- as.character(cut(pic.0h.genomic$CGI_SIZE.kb,
                                                  breaks=c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0)))
pic.0h.genomic$CGI_SIZE.group[(pic.0h.genomic$CGI_SIZE.kb <= 0.5)] <- "(0,0.5]"
pic.0h.genomic$CGI_SIZE.group[(pic.0h.genomic$CGI_SIZE.kb > 3.0)] <- "(3,4.585]"
pic.0h.genomic$CGI_SIZE.group[is.na(pic.0h.genomic$CGI_SIZE.group) | pic.0h.genomic$CGI_SIZE.kb == 0] <- "Absent"
pic.0h.genomic$CGI_SIZE.group <- factor(pic.0h.genomic$CGI_SIZE.group,
                                        labels=c("Absent", "(0,0.5]", "(0.5,1]", "(1,1.5]",  "(1.5,2]", "(2,2.5]",
                                                 "(2.5,3]", "(3,4.585]"),
                                        levels=c("Absent", "(0,0.5]", "(0.5,1]", "(1,1.5]",  "(1.5,2]", "(2,2.5]",
                                                 "(2.5,3]", "(3,4.585]"))

pic.cv2.by.cgi <- ggplot(pic.0h.genomic[pic.0h.genomic$N_CpG == 1,],
       aes(y=CV2, x=CGI_SIZE.group, fill=t0_t1.Diff,
           colour=t0_t1.Diff)) +
  geom_jitter(alpha=0.5, 
              position=position_jitterdodge(jitter.width=0.5)) +
  #geom_boxplot(colour="black") + 
  theme_mike() +
  #scale_fill_Publication(values=c("darkblue", "yellow", "darkred")) +
  scale_colour_manual(values=c("darkblue", "yellow", "darkred")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
  labs(x="CpG island size interval (kb)",
       y=expression(paste("Unstimulated CV"^2)))

pic.cv2.by.cgi.cdf <- ggplot(pic.0h.genomic[pic.0h.genomic$N_CpG == 1,],
       aes(x=CV2,
           colour=t0_t1.Diff)) +
  stat_ecdf() +
  theme_mike() +
  #scale_fill_Publication(values=c("darkblue", "yellow", "darkred")) +
  scale_colour_manual(values=c("darkblue", "yellow", "darkred")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
  labs(x="CpG island size interval (kb)",
       y=expression(paste("Unstimulated CV"^2))) +
  facet_wrap(~CGI_SIZE.group,
             ncol=3)

## set up a kolmogorov-smirnoff test between the undiff and up-diff for each CpG island size interval
# (0,0.5]
int1.x <- pic.0h.genomic$CV2[pic.0h.genomic$CGI_SIZE.group == "(0,0.5]" & pic.0h.genomic$t0_t1.Diff == 1]
int1.y <- pic.0h.genomic$CV2[pic.0h.genomic$CGI_SIZE.group == "(0,0.5]" & pic.0h.genomic$t0_t1.Diff == 0]
ks.test(int1.x, int1.y)

# (0.5,1]
int2.x <- pic.0h.genomic$CV2[pic.0h.genomic$CGI_SIZE.group == "(0.5,1]" & pic.0h.genomic$t0_t1.Diff == 1]
int2.y <- pic.0h.genomic$CV2[pic.0h.genomic$CGI_SIZE.group == "(0.5,1]" & pic.0h.genomic$t0_t1.Diff == 0]
ks.test(int2.x, int2.y)

# (1,1.5]
int3.x <- pic.0h.genomic$CV2[pic.0h.genomic$CGI_SIZE.group == "(1,1.5]" & pic.0h.genomic$t0_t1.Diff == 1]
int3.y <- pic.0h.genomic$CV2[pic.0h.genomic$CGI_SIZE.group == "(1,1.5]" & pic.0h.genomic$t0_t1.Diff == 0]
ks.test(int3.x, int3.y)

# (1.5,2]
int4.x <- pic.0h.genomic$CV2[pic.0h.genomic$CGI_SIZE.group == "(1.5,2]" & pic.0h.genomic$t0_t1.Diff == 1]
int4.y <- pic.0h.genomic$CV2[pic.0h.genomic$CGI_SIZE.group == "(1.5,2]" & pic.0h.genomic$t0_t1.Diff == 0]
ks.test(int4.x, int4.y)

# (2,2.5]
int5.x <- pic.0h.genomic$CV2[pic.0h.genomic$CGI_SIZE.group == "(2,2.5]" & pic.0h.genomic$t0_t1.Diff == 1]
int5.y <- pic.0h.genomic$CV2[pic.0h.genomic$CGI_SIZE.group == "(2,2.5]" & pic.0h.genomic$t0_t1.Diff == 0]
ks.test(int5.x, int5.y)


#########
## PAM ##
#########
# calculate the variability in expression at the 0h

pam.0h.cells <- na.omit(bmdc.meta$Sample[bmdc.meta$Timepoint == "0h"])
pam.0h.mean <- rowMeans(bmdc.pam.exprs[, pam.0h.cells])
pam.0h.var <- apply(bmdc.pam.exprs[, pam.0h.cells], 1, var)
pam.0h.cv2 <- pam.0h.var/(pam.0h.mean**2)

# calculate the residual CV^2
useForFit <- is.na(pam.0h.cv2) | pam.0h.mean < 0.1 | is.na(1/pam.0h.mean)

smoothScatter(x=pam.0h.mean[!useForFit], y=pam.0h.cv2[!useForFit])

# fit with a gamma-distributed GLM
fit <- glmgam.fit(cbind(a0 = 1, a1tilde=1/pam.0h.mean[!useForFit]), 
                  pam.0h.cv2[!useForFit])
pam.0h.rCV2 <- abs(pam.0h.cv2[!useForFit] - fitted.values(fit))

# plot gene expression varibility Vs logFC
pam.0h.sum <- do.call(cbind.data.frame,
                      list("Mean"=pam.0h.mean[!useForFit],
                           "CV2"=pam.0h.cv2[!useForFit],
                           "rCV2"=pam.0h.rCV2,
                           "gene_id"=rownames(bmdc.pam.exprs)[!useForFit]))

# merge with genomic features
pam.0h.demerge <- merge(pam.0h.sum, bmdc.pam.de.res.0v1, by='gene_id')
pam.0h.genomic <- merge(pam.0h.demerge, mouse.genomic.features, by.x='gene_id', by.y='GENE')

# split CGIs into size bins on quartiles
pam.0h.genomic$CGI_SIZE.group <- as.character(cut(pam.0h.genomic$CGI_SIZE.kb,
                                                  breaks=c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0)))
pam.0h.genomic$CGI_SIZE.group[(pam.0h.genomic$CGI_SIZE.kb <= 0.5)] <- "(0,0.5]"
pam.0h.genomic$CGI_SIZE.group[(pam.0h.genomic$CGI_SIZE.kb > 3.0)] <- "(3,4.585]"
pam.0h.genomic$CGI_SIZE.group[is.na(pam.0h.genomic$CGI_SIZE.group) | pam.0h.genomic$CGI_SIZE.kb == 0] <- "Absent"
pam.0h.genomic$CGI_SIZE.group <- factor(pam.0h.genomic$CGI_SIZE.group,
                                        labels=c("Absent", "(0,0.5]", "(0.5,1]", "(1,1.5]",  "(1.5,2]", "(2,2.5]",
                                                 "(2.5,3]", "(3,4.585]"),
                                        levels=c("Absent", "(0,0.5]", "(0.5,1]", "(1,1.5]",  "(1.5,2]", "(2,2.5]",
                                                 "(2.5,3]", "(3,4.585]"))

pam.cv2.by.cgi <- ggplot(pam.0h.genomic[pam.0h.genomic$N_CpG == 1,],
       aes(y=CV2, x=CGI_SIZE.group, fill=t0_t1.Diff,
           colour=t0_t1.Diff)) +
  geom_jitter(alpha=0.5, 
              position=position_jitterdodge(jitter.width=0.5)) +
  #geom_boxplot(colour="black") + 
  theme_mike() +
  #scale_fill_Publication(values=c("darkblue", "yellow", "darkred")) +
  scale_colour_manual(values=c("darkblue", "yellow", "darkred")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
  labs(x="CpG island size interval (kb)",
       y=expression(paste("Unstimulated CV"^2)))

pam.rcv2.by.cgi <- ggplot(pam.0h.genomic[pam.0h.genomic$N_CpG == 1,],
                         aes(y=rCV2, x=CGI_SIZE.group, fill=t0_t1.Diff,
                             colour=t0_t1.Diff)) +
  geom_jitter(alpha=0.5, 
              position=position_jitterdodge(jitter.width=0.5)) +
  #geom_boxplot(colour="black") + 
  theme_mike() +
  #scale_fill_Publication(values=c("darkblue", "yellow", "darkred")) +
  scale_colour_manual(values=c("darkblue", "yellow", "darkred")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
  labs(x="CpG island size interval (kb)",
       y=expression(paste("Unstimulated Residual CV"^2)))

pam.cv2.by.cgi.cdf <- ggplot(pam.0h.genomic[pam.0h.genomic$N_CpG == 1,],
       aes(x=CV2,
           colour=t0_t1.Diff)) +
  stat_ecdf() +
  theme_mike() +
  #scale_fill_Publication(values=c("darkblue", "yellow", "darkred")) +
  scale_colour_manual(values=c("darkblue", "yellow", "darkred")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
  labs(x="CpG island size interval (kb)",
       y=expression(paste("Unstimulated CV"^2))) +
  facet_wrap(~CGI_SIZE.group,
             ncol=3)

## set up a kolmogorov-smirnoff test between the undiff and up-diff for each CpG island size interval
# (0,0.5]
int1.x <- pam.0h.genomic$rCV2[pam.0h.genomic$CGI_SIZE.group == "(0,0.5]" & pam.0h.genomic$t0_t1.Diff == 1]
int1.y <- pam.0h.genomic$rCV2[pam.0h.genomic$CGI_SIZE.group == "(0,0.5]" & pam.0h.genomic$t0_t1.Diff == 0]
ks.test(int1.x, int1.y, alternative="two.sided")

# (0.5,1]
int2.x <- pam.0h.genomic$rCV2[pam.0h.genomic$CGI_SIZE.group == "(0.5,1]" & pam.0h.genomic$t0_t1.Diff == 1]
int2.y <- pam.0h.genomic$rCV2[pam.0h.genomic$CGI_SIZE.group == "(0.5,1]" & pam.0h.genomic$t0_t1.Diff == 0]
ks.test(int2.x, int2.y, alternative="two.sided")

# (1,1.5]
int3.x <- pam.0h.genomic$rCV2[pam.0h.genomic$CGI_SIZE.group == "(1,1.5]" & pam.0h.genomic$t0_t1.Diff == 1]
int3.y <- pam.0h.genomic$rCV2[pam.0h.genomic$CGI_SIZE.group == "(1,1.5]" & pam.0h.genomic$t0_t1.Diff == 0]
ks.test(int3.x, int3.y, alternative="two.sided")

# (1.5,2]
int4.x <- pam.0h.genomic$rCV2[pam.0h.genomic$CGI_SIZE.group == "(1.5,2]" & pam.0h.genomic$t0_t1.Diff == 1]
int4.y <- pam.0h.genomic$rCV2[pam.0h.genomic$CGI_SIZE.group == "(1.5,2]" & pam.0h.genomic$t0_t1.Diff == 0]
ks.test(int4.x, int4.y, alternative="two.sided")

# (2,2.5]
int5.x <- pam.0h.genomic$rCV2[pam.0h.genomic$CGI_SIZE.group == "(2,2.5]" & pam.0h.genomic$t0_t1.Diff == 1]
int5.y <- pam.0h.genomic$rCV2[pam.0h.genomic$CGI_SIZE.group == "(2,2.5]" & pam.0h.genomic$t0_t1.Diff == 0]
ks.test(int5.x, int5.y, alternative="two.sided")


int_all.x <- pam.0h.genomic$rCV2[pam.0h.genomic$t0_t1.Diff == 1]
int_all.y <- pam.0h.genomic$rCV2[pam.0h.genomic$t0_t1.Diff == 0]
ks.test(int_all.x, int_all.y, alternative="two.sided")

# plot heatmap of ranked CpG island size
# these are displayed under the density plots
# to illustrate how there is an enrichment of short CpG island genes
#########
## LPS ##
#########
lps.diff.df <- do.call(cbind.data.frame,
                       list("SizeDiff"=bmdc.lps.sign.size_rank,
                            "Size"=bmdc.lps.res01.size_rank))
lps.diff.df$Comparison <- "Test"

lps.null.df <- do.call(cbind.data.frame,
                       list("SizeDiff"=rbinom(n=length(bmdc.lps.res01.size_rank),
                                           size=1, prob=0.5),
                            "Size"=bmdc.lps.res01.size_rank))
lps.null.df$Comparison <- "Null"

lps.binom.df <- do.call(rbind.data.frame,
                        list("0v1"=lps.diff.df,
                             "null"=lps.null.df))

lps.binom.df$Comparison <- factor(lps.binom.df$Comparison,
                                  levels=c("Null", "Test"),
                                  labels=c("Null", "Test"))

lps.0v1.heat <- ggplot(lps.binom.df,
                       aes(x=Size, y=Comparison)) +
  geom_tile(aes(fill=SizeDiff)) +
  scale_fill_gradient(low="grey", high="darkred") +
  scale_x_continuous(limits=c(0, 3), oob=censor) +
  theme_mike()  +
  theme(panel.grid=element_blank()) +
  labs(x="CpG island size (kb)", y="Comparison") +
  guides(fill=FALSE)

ggsave(lps.0v1.heat,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/BMDC_LPS-binom_heat.png",
       width=8.25, height=2.25, dpi=300)

##################################
# also check 0 vs 2 against 0v1? #
##################################

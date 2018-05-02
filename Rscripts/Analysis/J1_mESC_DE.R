# find the top 1000 DE genes from differentiating J1 mESCs
library(biomaRt)
library(limma)
library(Matrix)
j1.data <- read.table("~/Dropbox/mESC/GSE3749_AFFY.txt",
                      h=TRUE, sep="\t", stringsAsFactors=FALSE)

ensembl <- useEnsembl(biomart='ensembl', dataset='mmusculus_gene_ensembl', GRCh=37)

j1.gene_symbol <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 
                                     'affy_moe430a'),
                           filters='affy_moe430a', mart=ensembl,
                           values=j1.data$X)

# pull out the expression data
j1.exprs.cols <- colnames(j1.data)[grepl(colnames(j1.data), pattern="Expression")]
j1.exprs <- j1.data[, c("X", j1.exprs.cols)]

# S242 - 0 hour, S243 - 1 hour, ...
time.point <- unlist(lapply(strsplit(j1.exprs.cols, split=".", fixed=TRUE),
                            FUN=function(X) paste0(X[3])))
rep.stem <- unlist(lapply(strsplit(j1.exprs.cols, split=".", fixed=TRUE),
                          FUN=function(X) paste0(X[4])))

j1.meta <- do.call(cbind.data.frame,
                   list("Sample"=j1.exprs.cols,
                        "TimeID"=time.point,
                        "Replicate"=rep.stem))
j1.meta$TimeID <- as.character(j1.meta$TimeID)

j1.meta$Timepoint <- "0h"
j1.meta$Timepoint[j1.meta$TimeID == "S242"] <- "0h"
j1.meta$Timepoint[j1.meta$TimeID == "S243"] <- "6h"
j1.meta$Timepoint[j1.meta$TimeID == "S244"] <- "12h"
j1.meta$Timepoint[j1.meta$TimeID == "S245"] <- "36h"
j1.meta$Timepoint[j1.meta$TimeID == "S246"] <- "4d"
j1.meta$Timepoint[j1.meta$TimeID == "S247"] <- "7d"
j1.meta$Timepoint[j1.meta$TimeID == "S248"] <- "9d"
j1.meta$Timepoint[j1.meta$TimeID == "S249"] <- "14d"
j1.meta$Timepoint[j1.meta$TimeID == "S250"] <- "48h"
j1.meta$Timepoint[j1.meta$TimeID == "S251"] <- "18h"
j1.meta$Timepoint[j1.meta$TimeID == "S252"] <- "24h"

j1.meta$Timepoint <- ordered(j1.meta$Timepoint,
                             levels=c("0h", "6h", "12h", "18h", "24h", "36h", "48h", "4d", "7d", "9d", "14d"))

rownames(j1.meta) <- j1.meta$Sample
rownames(j1.exprs) <- j1.exprs$X

# log10 transform
j1.exprs <- log10(j1.exprs[, -1])

# DE testing with limma
j1.design <- model.matrix(~ Timepoint, data=j1.meta)

# fit a linear model
j1.fit <- lmFit(j1.exprs, j1.design)
j1.fit <- eBayes(j1.fit)
sum.res.j1 <- summary(decideTests(j1.fit))

de.res.j1 <- topTable(j1.fit, coef=2, n=Inf, sort="p", p=1.0)
de.res.j1$Sig <- 0
de.res.j1$Sig[de.res.j1$adj.P.Val <= 0.01] <- 1
de.res.j1$Sig <- as.factor(de.res.j1$Sig)

de.res.j1$Diff <- 0
de.res.j1$Diff[de.res.j1$logFC < 0 & de.res.j1$Sig == 1] <- -1
de.res.j1$Diff[de.res.j1$logFC > 0 & de.res.j1$Sig == 1] <- 1

de.res.j1$affy_id <- rownames(de.res.j1)

de.res.j1.merge <- merge(de.res.j1, j1.gene_symbol, by.x='affy_id', by.y='affy_moe430a')

devel_genes <- de.res.j1.merge[de.res.j1.merge$Diff == 1, ]





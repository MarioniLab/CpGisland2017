## Investigating genomic factors that influence gene expression noise
library(ggplot2)
library(reshape2)
library(Rtsne)
library(glmnet)
library(biomaRt)
source("~/Dropbox/R_sessions/GGMike/palette_256.R")
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

th.cells <- read.table("~/Dropbox/Th_data/ThCell_norm.tsv",
                        sep="\t", h=T, stringsAsFactors=F)

rownames(th.cells) <- th.cells$gene_id
th.cells <- th.cells[grepl(rownames(th.cells), pattern="ENS"), 1:(dim(th.cells)[2]-1)]

th.meta <- read.table("~/Dropbox/Th_data/ThCell-meta.tsv",
                       sep="\t", h=T, stringsAsFactors=F)

# need to calculate noise features within each stimulation time point and condition
# ignore the cross information for now as the number of cells for some is v.low
# might just have to adjust in DE linear models
clusters <- expand.grid(unique(th.meta$Condition), unique(th.meta$Timepoint))
clusters[, 1] <- as.character(clusters[, 1])
clusters[, 2] <- as.character(clusters[, 2])
cluster.list <- list()

for(i in seq_along(clusters[, 1])){
  clust <- paste(clusters[i, ], collapse="_")
  clust.cells <- th.meta$Sample[(th.meta$Condition == clusters[i, 1]) & (th.meta$Timepoint == clusters[i, 2])]
  clust.df <- th.cells[, colnames(th.cells) %in% clust.cells]
  print(dim(clust.df))
  if(dim(clust.df)[2] > 10){
    th.means <- rowMeans(clust.df)
    th.vars <- apply(clust.df,
                      1, var)
    th.median <- apply(clust.df,
                        1, median)
    th.mad <- apply(clust.df,
                     1, FUN=function(M) mad(M, constant=1))

    # calculate correlation between SP1 expression and all other genes.
    sp1.exprs <- clust.df["ENSMUSG00000001280",]
    sp3.exprs <- clust.df["ENSMUSG00000027109",]

    sp1.cor <- apply(clust.df,
                     1, FUN=function(Q) cor(Q, t(sp1.exprs)))

    sp1.sp3.ratio <- sp1.exprs/sp3.exprs
    sp1.sp3.ratio[is.infinite(unlist(sp1.sp3.ratio))] <- 0
    sp1.sp3.ratio[is.na(unlist(sp1.sp3.ratio))] <- 0
    sp.ratio.cor <- apply(clust.df,
                          1, FUN=function(Q) cor(Q, t(sp1.sp3.ratio)))

    # create gene expression groups based on average expression over cells
    th.gene.summary <- as.data.frame(cbind(th.means, th.vars, th.median, th.mad, sp1.cor, sp.ratio.cor))
    colnames(th.gene.summary) <- c("Mean", "Var", "Median", "MAD", "SP1.COR", "SP.RATIO")
    th.gene.summary$CV2 <- th.gene.summary$Var/(th.gene.summary$Mean ** 2)
    th.gene.summary$GENE <- rownames(th.gene.summary)
    th.gene.summary <- th.gene.summary[(!th.gene.summary$Mean <= 0.05), ]
    cluster.list[[clust]] <- th.gene.summary
  }
}



## Investigating genomic factors that influence gene expression noise
library(ggplot2)
library(reshape2)
library(Rtsne)
library(glmnet)
library(biomaRt)
source("~/Dropbox/R_sessions/GGMike/palette_256.R")
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

cd8.cells <- read.table("~/Dropbox/Tcell/CD8_Tcell_norm.tsv",
                         sep="\t", h=T, stringsAsFactors=F)
rownames(cd8.cells) <- cd8.cells$gene_id
cd8.cells <- cd8.cells[grepl(rownames(cd8.cells), pattern="ENS"), 1:(dim(cd8.cells)[2]-1)]

cd8.meta <- read.table("~/Dropbox/Tcell/CD8_Tcell-meta.tsv",
                       sep="\t", h=T, stringsAsFactors=F)
cd8.meta$Sample <- gsub(cd8.meta$Sample, pattern="-", replacement=".")

# need to calculate noise features within each stimulation time point
clusters <- unique(cd8.meta$Condition)
cluster.list <- list()

for(i in seq_along(clusters)){
  clust <- clusters[i]
  clust.cells <- cd8.meta$Sample[cd8.meta$Condition == clust]
  clust.df <- cd8.cells[, colnames(cd8.cells) %in% clust.cells]
  if(dim(clust.df)[2] > 10){
    cd8.means <- rowMeans(clust.df)
    cd8.vars <- apply(clust.df,
                       1, var)
    cd8.median <- apply(clust.df,
                         1, median)
    cd8.mad <- apply(clust.df,
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
    cd8.gene.summary <- as.data.frame(cbind(cd8.means, cd8.vars, cd8.median, cd8.mad, sp1.cor, sp.ratio.cor))
    colnames(cd8.gene.summary) <- c("Mean", "Var", "Median", "MAD", "SP1.COR", "SP.RATIO")
    cd8.gene.summary$CV2 <- cd8.gene.summary$Var/(cd8.gene.summary$Mean ** 2)
    cd8.gene.summary$GENE <- rownames(cd8.gene.summary)
    cd8.gene.summary <- cd8.gene.summary[(!cd8.gene.summary$Mean <= 0.05), ]
    cluster.list[[clust]] <- cd8.gene.summary
  }
}

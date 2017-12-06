source("~/Dropbox/R_sessions/Noise/mESC_chromHMM.R")
source("~/Dropbox/R_sessions/Noise/human_genomic_noise_features.R")
source("~/Dropbox/R_sessions/Noise/genomic_noise_features.R")

library(pheatmap)
mouse.genomic.features$CGI_SIZE.kb <- mouse.genomic.features$CGI_SIZE/1000
human.genomic.features$CGI_SIZE.kb <- human.genomic.features$CGI_SIZE/1000

# plot CpG island features as a heatmap
mouse.cor <- cor(mouse.genomic.features[, c("CGI_SIZE.kb", "cpg_Overlap", "cpg_RATIO", "SP1")], method="spearman")
colnames(mouse.cor) <- c("CpG island size", "CpG island overlap", "CpG dinucleotide ratio", "Number SP1 motifs")
rownames(mouse.cor) <- c("CpG island size", "CpG island overlap", "CpG dinucleotide ratio", "Number SP1 motifs")

human.cor <- cor(human.genomic.features[, c("CGI_SIZE.kb", "cpg_Overlap", "cpg_RATIO", "SP1")], method="spearman")
colnames(human.cor) <- c("CpG island size", "CpG island overlap", "CpG dinucleotide ratio", "Number SP1 motifs")
rownames(human.cor) <- c("CpG island size", "CpG island overlap", "CpG dinucleotide ratio", "Number SP1 motifs")

pheatmap(mouse.cor,
         filename="~/Dropbox/Noise_genomics/Figures/ms_figures/Supp_mouse_CGI_features-heatmap.png",
         height=5.75, width=5.75, res=300)
pheatmap(human.cor,
         filename="~/Dropbox/Noise_genomics/Figures/ms_figures/Supp_human_CGI_features-heatmap.png",
         height=5.75, width=5.75, res=300)

source("~/Dropbox/R_sessions/Noise/mESC_chromHMM.R")
source("~/Dropbox/R_sessions/Noise/human_genomic_noise_features.R")
source("~/Dropbox/R_sessions/Noise/hESC_genomic_noise_features.R")

library(ggplot2)
library(robustbase)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

human.genomic.features$CGI_SIZE.kb <- human.genomic.features$CGI_SIZE/1000

hesc.vars <- colnames(human.genomic.features)
hesc.vars <- hesc.vars[grepl(hesc.vars, pattern="(CGI_SIZE.kb)|(TOTLEN)|(SP1)")]

hesc.match <- merge(hesc.gene.summary, human.genomic.features,
                    by='GENE')
hesc.match$Recip.Mean <- 1/hesc.match$Mean

# select only those genes with a CpG island
#hesc.match <- hesc.match[hesc.match$N_CpG == 1, ]

#########################################
## univariate robust linear regression ##
#########################################
# this is plotting for presentation/manuscript figures
hesc.univariate_list <- list()
model.control <- lmrob.control(max.it=500, k.max=500, rel.tol=1e-7)

for(x in seq_along(hesc.vars)){
  #.variable <- paste(c("Recip.Mean", hesc.vars[x]), collapse=" + ")
  .variable <- hesc.vars[x]
  .glm.form <- as.formula(paste("Residual.CV2", .variable, sep=" ~ "))
  
  m.rlm <- lmrob(.glm.form, data=hesc.match, control=model.control)
  m.robust <- summary(m.rlm)
  m.rlm.res <- as.data.frame(m.robust$coefficients)
  m.rlm.res$Sig <- as.numeric(m.rlm.res$`Pr(>|t|)` <= 0.05)
  m.res.mat <- as.matrix(m.rlm.res)
  m.res.var <- m.res.mat[2, ]
  names(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
  hesc.univariate_list[[hesc.vars[x]]] <- as.list(m.res.var)
}

# # add the mean expression on it's own
# .glm.form <- as.formula(paste("Residual.CV2", "Recip.Mean", sep=" ~ "))
# 
# m.rlm <- rlm(.glm.form, data=hesc.match)
# m.robust <- summary(m.rlm)
# m.rlm.res <- as.data.frame(m.robust$coefficients)
# m.rlm.res$Pval <- 2*pt(-abs(m.rlm.res[, 3]), df=3)
# m.rlm.res$Sig <- as.numeric(m.rlm.res$Pval <= 0.05)
# m.res.mat <- as.matrix(m.rlm.res)
# m.res.var <- m.res.mat[2, ]
# names(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
# hesc.univariate_list[["Recip.Mean"]] <- as.list(m.res.var)

hesc.rlm.df <- do.call(rbind.data.frame, hesc.univariate_list)
hesc.rlm.df$Predictor <- rownames(hesc.rlm.df)
hesc.rlm.df$Padjust <- p.adjust(hesc.rlm.df$P)
hesc.rlm.df$Sig <- as.numeric(hesc.rlm.df$Padjust <= 0.01)
hesc.rlm.df$Direction <- "NoEffect"
hesc.rlm.df$Direction[hesc.rlm.df$COEFF < 0 & hesc.rlm.df$Sig == 1] <- "Less"
hesc.rlm.df$Direction[hesc.rlm.df$COEFF > 0 & hesc.rlm.df$Sig == 1] <- "More"

# give the features more informative/better formated names
hesc.rlm.df$Predictor[hesc.rlm.df$Predictor == "Recip.Mean"] <- "Mean expression"
hesc.rlm.df$Predictor[hesc.rlm.df$Predictor == "N_CpG"] <- "CpG island"
hesc.rlm.df$Predictor[hesc.rlm.df$Predictor == "cpg_GCNUM"] <- "CpG island GC content"
hesc.rlm.df$Predictor[hesc.rlm.df$Predictor == "SP1"] <- "Number SP1 motifs"
hesc.rlm.df$Predictor[hesc.rlm.df$Predictor == "cpg_Overlap"] <- "CpG island overlap"
hesc.rlm.df$Predictor[hesc.rlm.df$Predictor == "cpg_RATIO"] <- "CpG dinucleotide ratio"
hesc.rlm.df$Predictor[hesc.rlm.df$Predictor == "CGI_SIZE.kb"] <- "CpG island size"
hesc.rlm.df$Predictor[hesc.rlm.df$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"

hesc.rlm.df$Tissue <- "ESC"
hesc.rlm.df$Species <- "Human"
hesc.rlm.df$Model <- "Univariate"

effect.cols <- c("#62148f", "#878787", "#feaf10")
names(effect.cols) <- c("Less", "NoEffect", "More")

write.table(hesc.rlm.df,
            file="~/Dropbox/Noise_genomics/Model_results/human_ESC_CGI-univariateRLM.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")

cpg.plot <- ggplot(hesc.rlm.df,
                   aes(x=reorder(Predictor, -STAT),
                       y=STAT, fill=Direction)) +
  geom_point(alpha=0.55, shape=21, size=5) + theme_mike() +
  scale_fill_manual(values=effect.cols) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Annotation", y="t-statistic") +
  guides(fill=FALSE) +
  scale_y_continuous(limits=c(-50, 50), oob=squish) +
  geom_hline(mapping=aes(yintercept=0), linetype="dashed", colour="grey") 

ggsave(cpg.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/hESC_CGI-univariateLM.png",
       height=4.75, width=6.75, dpi=300)

###########################################
## multivariate robust linear regression ##
###########################################
hesc.genomic.vars <- paste(hesc.vars,
                           collapse=" + ")

hesc.glm.form <- as.formula(paste("Residual.CV2",
                                  hesc.genomic.vars, sep=" ~ "))

hesc.rlm <- lmrob(hesc.glm.form, data=hesc.match, control=model.control)
hesc.robust <- summary(hesc.rlm)
hesc.rlm.res <- as.data.frame(hesc.robust$coefficients)
hesc.rlm.res$Padjust <- p.adjust(hesc.rlm.res$`Pr(>|t|)`)
hesc.rlm.res$Sig <- as.numeric(hesc.rlm.res$Padjust <= 0.01)
hesc.rlm.res$Predictor <- rownames(hesc.rlm.res)
colnames(hesc.rlm.res) <- c("COEFF", "SE", "STAT", "P", "Padjust", "Sig", "Predictor")
hesc.rlm.res$Direction <- "NoEffect"
hesc.rlm.res$Direction[hesc.rlm.res$COEFF < 0 & hesc.rlm.res$Sig == 1] <- "Less"
hesc.rlm.res$Direction[hesc.rlm.res$COEFF > 0 & hesc.rlm.res$Sig == 1] <- "More"

# give the features more informative/better formated names
hesc.rlm.res$Predictor[hesc.rlm.res$Predictor == "Recip.Mean"] <- "Mean expression"
hesc.rlm.res$Predictor[hesc.rlm.res$Predictor == "N_CpG"] <- "CpG island"
hesc.rlm.res$Predictor[hesc.rlm.res$Predictor == "cpg_GCNUM"] <- "CpG island GC content"
hesc.rlm.res$Predictor[hesc.rlm.res$Predictor == "SP1"] <- "Number SP1 motifs"
hesc.rlm.res$Predictor[hesc.rlm.res$Predictor == "cpg_Overlap"] <- "CpG island overlap"
hesc.rlm.res$Predictor[hesc.rlm.res$Predictor == "cpg_RATIO"] <- "CpG dinucleotide ratio"
hesc.rlm.res$Predictor[hesc.rlm.res$Predictor == "CGI_SIZE.kb"] <- "CpG island size"
hesc.rlm.res$Predictor[hesc.rlm.res$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"

# remove the intercept term
hesc.rlm.res <- hesc.rlm.res[!grepl(hesc.rlm.res$Predictor, pattern="Intercept"), ]

hesc.rlm.res$Tissue <- "ESC"
hesc.rlm.res$Species <- "Human"
hesc.rlm.res$Model <- "Multivariate"

write.table(hesc.rlm.res,
            file="~/Dropbox/Noise_genomics/Model_results/human_ESC_CGI-multivariateRLM.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")

cpg.plot <- ggplot(hesc.rlm.res,
                   aes(x=reorder(Predictor, -STAT),
                       y=STAT, fill=Direction)) +
  geom_point(alpha=0.55, shape=21, size=5) + theme_mike() +
  scale_fill_manual(values=effect.cols) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Annotation", y="t-statistic") +
  guides(fill=FALSE) +
  scale_y_continuous(limits=c(-50, 50), oob=squish) +
  geom_hline(mapping=aes(yintercept=0), linetype="dashed", colour="grey") 

ggsave(cpg.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/hESC_CGI-multivariateLM.png",
       height=4.75, width=6.75, dpi=300)

# plot uni and multivariate results together
hesc.rlm.all <- do.call(rbind.data.frame, list("Univariate"=hesc.rlm.df,
                                               "Multivariate"=hesc.rlm.res))

all.plot <- ggplot(hesc.rlm.all,
                   aes(x=reorder(Predictor, -STAT),
                       y=STAT, fill=Direction, shape=Model)) +
  geom_hline(mapping=aes(yintercept=0), linetype="dashed", colour="grey") +
  geom_jitter(size=4, alpha=0.8,
              position=position_jitterdodge(jitter.height=0,
                                            jitter.width=0.1)) +
  theme_mike() +
  scale_fill_manual(values=c("#62148f", "#feaf10", "#878787")) +
  scale_shape_manual(values=c(21, 23)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Annotation", y="t-statistic") +
  guides(fill=FALSE, shape=FALSE) +
  scale_y_continuous(limits=c(-30, 30), oob=squish)

ggsave(all.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-ESC_CGI-allLM.png",
       height=5.6, width=7.75, dpi=300)
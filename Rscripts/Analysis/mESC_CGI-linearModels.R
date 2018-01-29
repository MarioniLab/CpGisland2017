source("~/Dropbox/R_sessions/Noise/mESC_chromHMM.R")
source("~/Dropbox/R_sessions/Noise/genomic_noise_features.R")
source("~/Dropbox/R_sessions/Noise/mESC_genomic_noise_features.R")

library(ggplot2)
library(scales)
library(robustbase)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

genomic.features$CGI_SIZE.kb <- genomic.features$CGI_SIZE/1000

mesc.vars <- colnames(genomic.features)

# there is considerable correlation between CGI features,
# just pick the island length for now
mesc.vars <- mesc.vars[grepl(mesc.vars, pattern="(CGI_SIZE.kb)|(TOTLEN)|(SP1)")]

mesc.match <- merge(mesc.gene.summary, genomic.features,
                     by='GENE')
mesc.match$Recip.Mean <- 1/mesc.match$Mean

# select only those genes with a CpG island
#mesc.match <- mesc.match[mesc.match$N_CpG == 1, ]

#########################################
## univariate robust linear regression ##
#########################################
# this is plotting for presentation/manuscript figures
mesc.univariate_list <- list()
model.control <- lmrob.control(max.it=500, k.max=500, rel.tol=1e-7)

for(x in seq_along(mesc.vars)){
  #.variable <- paste(c("Mean", mesc.vars[x]), collapse=" + ")
  .variable <- mesc.vars[x]
  .glm.form <- as.formula(paste("Residual.CV2", .variable, sep=" ~ "))
  
  m.rlm <- lmrob(.glm.form, data=mesc.match, control=model.control)
  m.robust <- summary(m.rlm)
  m.rlm.res <- as.data.frame(m.robust$coefficients)
  m.rlm.res$Sig <- as.numeric(m.rlm.res$`Pr(>|t|)` <= 0.05)
  m.res.mat <- as.matrix(m.rlm.res)
  m.res.var <- m.res.mat[2, ]
  names(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
  mesc.univariate_list[[mesc.vars[x]]] <- as.list(m.res.var)
}

# # add the mean expression on it's own
# .glm.form <- as.formula(paste("Residual.CV2", "Recip.Mean", sep=" ~ "))
# 
# m.rlm <- lmrob(.glm.form, data=mesc.match, control=model.control)
# m.robust <- summary(m.rlm)
# m.rlm.res <- as.data.frame(m.robust$coefficients)
# m.rlm.res$Sig <- as.numeric(m.rlm.res$`Pr(>|t|)` <= 0.05)
# m.res.mat <- as.matrix(m.rlm.res)
# m.res.var <- m.res.mat[2, ]
# names(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
# mesc.univariate_list[["Recip.Mean"]] <- as.list(m.res.var)


mesc.rlm.df <- do.call(rbind.data.frame, mesc.univariate_list)
mesc.rlm.df$Predictor <- rownames(mesc.rlm.df)
mesc.rlm.df$Padjust <- p.adjust(mesc.rlm.df$P)
mesc.rlm.df$Sig <- as.numeric(mesc.rlm.df$Padjust <= 0.01)
mesc.rlm.df$Direction <- "NoEffect"
mesc.rlm.df$Direction[mesc.rlm.df$COEFF < 0 & mesc.rlm.df$Sig == 1] <- "Less"
mesc.rlm.df$Direction[mesc.rlm.df$COEFF > 0 & mesc.rlm.df$Sig == 1] <- "More"

# give the features more informative/better formated names
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "Recip.Mean"] <- "Mean expression"
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "N_CpG"] <- "CpG island"
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "cpg_GCNUM"] <- "CpG island GC content"
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "SP1"] <- "Number SP1 motifs"
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "cpg_Overlap"] <- "CpG island overlap"
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "cpg_RATIO"] <- "CpG dinucleotide ratio"
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "CGI_SIZE.kb"] <- "CpG island size"
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"

mesc.rlm.df$Tissue <- "ESC"
mesc.rlm.df$Species <- "Mouse"
mesc.rlm.df$Model <- "Univariate"

effect.cols <- c("#62148f", "#878787", "#feaf10")
names(effect.cols) <- c("Less", "NoEffect", "More")

write.table(mesc.rlm.df,
            file="~/Dropbox/Noise_genomics/Model_results/mouse_ESC_CGI-univariateRLM.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")

cpg.plot <- ggplot(mesc.rlm.df,
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
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mESC_CGI-univariateLM.png",
       height=4.75, width=6.75, dpi=300)

###########################################
## multivariate robust linear regression ##
###########################################
mesc.genomic.vars <- paste(c(mesc.vars),
                            collapse=" + ")

mesc.glm.form <- as.formula(paste("Residual.CV2",
                                   mesc.genomic.vars, sep=" ~ "))

mesc.rlm <- lmrob(mesc.glm.form, data=mesc.match, control=model.control)
mesc.robust <- summary(mesc.rlm)
mesc.rlm.res <- as.data.frame(mesc.robust$coefficients)
#mesc.rlm.res$Pval <- 2*pt(-abs(mesc.rlm.res[, 3]), df=dim(mesc.match)[2]-1)
mesc.rlm.res$Padjust <- p.adjust(mesc.rlm.res$`Pr(>|t|)`)
mesc.rlm.res$Sig <- as.numeric(mesc.rlm.res$Padjust <= 0.01)
mesc.rlm.res$Predictor <- rownames(mesc.rlm.res)
colnames(mesc.rlm.res) <- c("COEFF", "SE", "STAT", "P", "Padjust", "Sig", "Predictor")
mesc.rlm.res$Direction <- "NoEffect"
mesc.rlm.res$Direction[mesc.rlm.res$COEFF < 0 & mesc.rlm.res$Sig == 1] <- "Less"
mesc.rlm.res$Direction[mesc.rlm.res$COEFF > 0 & mesc.rlm.res$Sig == 1] <- "More"

# give the features more informative/better formated names
mesc.rlm.res$Predictor[mesc.rlm.res$Predictor == "Recip.Mean"] <- "Mean expression"
mesc.rlm.res$Predictor[mesc.rlm.res$Predictor == "N_CpG"] <- "CpG island"
mesc.rlm.res$Predictor[mesc.rlm.res$Predictor == "cpg_GCNUM"] <- "CpG island GC content"
mesc.rlm.res$Predictor[mesc.rlm.res$Predictor == "SP1"] <- "Number SP1 motifs"
mesc.rlm.res$Predictor[mesc.rlm.res$Predictor == "cpg_Overlap"] <- "CpG island overlap"
mesc.rlm.res$Predictor[mesc.rlm.res$Predictor == "cpg_RATIO"] <- "CpG dinucleotide ratio"
mesc.rlm.res$Predictor[mesc.rlm.res$Predictor == "CGI_SIZE.kb"] <- "CpG island size"
mesc.rlm.res$Predictor[mesc.rlm.res$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"

# remove the intercept term
mesc.rlm.res <- mesc.rlm.res[!grepl(mesc.rlm.res$Predictor, pattern="Intercept"), ]

mesc.rlm.res$Tissue <- "ESC"
mesc.rlm.res$Species <- "Mouse"
mesc.rlm.res$Model <- "Multivariate"

write.table(mesc.rlm.res,
            file="~/Dropbox/Noise_genomics/Model_results/mouse_ESC_CGI-multivariateRLM.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")

cpg.plot <- ggplot(mesc.rlm.res,
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
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mESC_CGI-univariateLM.png",
       height=4.75, width=6.75, dpi=300)

# plot uni and multivariate results together
mesc.rlm.all <- do.call(rbind.data.frame, list("Univariate"=mesc.rlm.df,
                                                "Multivariate"=mesc.rlm.res))

all.plot <- ggplot(mesc.rlm.all,
                   aes(x=reorder(Predictor, -STAT),
                       y=STAT, fill=Direction, shape=Model)) +
  geom_hline(mapping=aes(yintercept=0), linetype="dashed", colour="grey") +
  geom_jitter(size=4, alpha=0.8,
              position=position_jitterdodge(jitter.height=0,
                                            jitter.width=0.1)) +
  theme_mike() +
  scale_fill_manual(values=c("#62148f", "#feaf10", "#878787")) +
  scale_shape_manual(values=c(21, 23)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        axis.text=element_text(size=16)) +
  labs(x="Annotation", y="t-statistic") +
  guides(fill=FALSE, shape=FALSE) +
  scale_y_continuous(limits=c(-30, 30), oob=squish)

ggsave(all.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mESC_CGI-allLM.png",
       height=5.6, width=7.75, dpi=300)

# CpG island features across different cell types
#source("~/Dropbox/R_sessions/Noise/mESC_chromHMM.R")
source("~/Dropbox/R_sessions/Noise/human_genomic_noise_features.R")
source("~/Dropbox/R_sessions/Noise/human_beta_islet_genomic_noise_features.R")
library(scales)
library(ggplot2)
library(robustbase)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

human.genomic.features$CGI_SIZE.kb <- human.genomic.features$CGI_SIZE/1000

beta.vars <- colnames(human.genomic.features)
beta.vars <- beta.vars[grepl(beta.vars, pattern="(CGI_SIZE.kb)|(TOTLEN)|(SP1)")]

beta.match <- merge(beta.gene.summary, human.genomic.features,
                    by='GENE')
beta.match$Recip.Mean <- 1/beta.match$Mean

# select only those genes with a CpG island
#beta.match <- beta.match[beta.match$N_CpG == 1, ]

#########################################
## univariate robust linear regression ##
#########################################
# this is plotting for presentation/manuscript figures
beta.univariate_list <- list()
model.control <- lmrob.control(max.it=500, k.max=500, rel.tol=1e-7)

for(x in seq_along(beta.vars)){
  #.variable <- paste(c("Recip.Mean", beta.vars[x]), collapse=" + ")
  .variable <- beta.vars[x]
  .glm.form <- as.formula(paste("Residual.CV2", .variable, sep=" ~ "))
  
  m.rlm <- lmrob(.glm.form, data=beta.match, control=model.control)
  m.robust <- summary(m.rlm)
  m.rlm.res <- as.data.frame(m.robust$coefficients)
  #m.rlm.res$Pval <- 2*pt(-abs(m.rlm.res[, 3]), df=3)
  m.rlm.res$Sig <- as.numeric(m.rlm.res$`Pr(>|t|)` <= 0.05)
  m.res.mat <- as.matrix(m.rlm.res)
  m.res.var <- m.res.mat[2, ]
  names(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
  beta.univariate_list[[beta.vars[x]]] <- as.list(m.res.var)
}

# # add the mean expression on it's own
# .glm.form <- as.formula(paste("Residual.CV2", "Recip.Mean", sep=" ~ "))
# 
# m.rlm <- rlm(.glm.form, data=beta.match)
# m.robust <- summary(m.rlm)
# m.rlm.res <- as.data.frame(m.robust$coefficients)
# m.rlm.res$Pval <- 2*pt(-abs(m.rlm.res[, 3]), df=3)
# m.rlm.res$Sig <- as.numeric(m.rlm.res$Pval <= 0.05)
# m.res.mat <- as.matrix(m.rlm.res)
# m.res.var <- m.res.mat[2, ]
# names(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
# beta.univariate_list[["Recip.Mean"]] <- as.list(m.res.var)

beta.rlm.df <- do.call(rbind.data.frame, beta.univariate_list)
beta.rlm.df$Predictor <- rownames(beta.rlm.df)
beta.rlm.df$Padjust <- p.adjust(beta.rlm.df$P)
beta.rlm.df$Sig <- as.numeric(beta.rlm.df$Padjust <= 0.01)
beta.rlm.df$Direction <- "NoEffect"
beta.rlm.df$Direction[beta.rlm.df$COEFF < 0 & beta.rlm.df$Sig == 1] <- "Less"
beta.rlm.df$Direction[beta.rlm.df$COEFF > 0 & beta.rlm.df$Sig == 1] <- "More"

# give the features more informative/better formated names
beta.rlm.df$Predictor[beta.rlm.df$Predictor == "Recip.Mean"] <- "Mean expression"
beta.rlm.df$Predictor[beta.rlm.df$Predictor == "N_CpG"] <- "CpG island"
beta.rlm.df$Predictor[beta.rlm.df$Predictor == "cpg_GCNUM"] <- "CpG island GC content"
beta.rlm.df$Predictor[beta.rlm.df$Predictor == "SP1"] <- "Number SP1 motifs"
beta.rlm.df$Predictor[beta.rlm.df$Predictor == "cpg_Overlap"] <- "CpG island overlap"
beta.rlm.df$Predictor[beta.rlm.df$Predictor == "cpg_RATIO"] <- "CpG dinucleotide ratio"
beta.rlm.df$Predictor[beta.rlm.df$Predictor == "CGI_SIZE.kb"] <- "CpG island size"
beta.rlm.df$Predictor[beta.rlm.df$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"

beta.rlm.df$Tissue <- "BetaIslet"
beta.rlm.df$Species <- "Human"
beta.rlm.df$Model <- "Univariate"

effect.cols <- c("#62148f", "#878787", "#feaf10")
names(effect.cols) <- c("Less", "NoEffect", "More")

write.table(beta.rlm.df,
            file="~/Dropbox/Noise_genomics/Model_results/human_BetaIslet_CGI-univariateRLM.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

cpg.plot <- ggplot(beta.rlm.df,
                   aes(x=reorder(Predictor, -STAT),
                       y=STAT, fill=Direction)) +
  geom_point(alpha=0.55, shape=21, size=5) + theme_mike() +
  scale_fill_manual(values=effect.cols) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Annotation", y="t-statistic") +
  guides(fill=FALSE) +
  scale_y_continuous(limits=c(-10, 10), oob=squish) +
  geom_hline(mapping=aes(yintercept=0), linetype="dashed", colour="grey") 

ggsave(cpg.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-BetaIslet_CGI-univariateLM.png",
       height=5.25, width=11.75, dpi=300)

###########################################
## multivariate robust linear regression ##
###########################################
beta.genomic.vars <- paste(beta.vars,
                           collapse=" + ")

beta.glm.form <- as.formula(paste("Residual.CV2",
                                  beta.genomic.vars, sep=" ~ "))

beta.rlm <- lmrob(beta.glm.form, data=beta.match, control=model.control)
beta.robust <- summary(beta.rlm)
beta.rlm.res <- as.data.frame(beta.robust$coefficients)
beta.rlm.res$Padjust <- p.adjust(beta.rlm.res$`Pr(>|t|)`)
beta.rlm.res$Sig <- as.numeric(beta.rlm.res$Padjust <= 0.01)
beta.rlm.res$Predictor <- rownames(beta.rlm.res)
colnames(beta.rlm.res) <- c("COEFF", "SE", "STAT", "P", "Padjust", "Sig", "Predictor")
beta.rlm.res$Direction <- "NoEffect"
beta.rlm.res$Direction[beta.rlm.res$COEFF < 0 & beta.rlm.res$Sig == 1] <- "Less"
beta.rlm.res$Direction[beta.rlm.res$COEFF > 0 & beta.rlm.res$Sig == 1] <- "More"

# give the features more informative/better formated names
beta.rlm.res$Predictor[beta.rlm.res$Predictor == "Recip.Mean"] <- "Mean expression"
beta.rlm.res$Predictor[beta.rlm.res$Predictor == "N_CpG"] <- "CpG island"
beta.rlm.res$Predictor[beta.rlm.res$Predictor == "cpg_GCNUM"] <- "CpG island GC content"
beta.rlm.res$Predictor[beta.rlm.res$Predictor == "SP1"] <- "Number SP1 motifs"
beta.rlm.res$Predictor[beta.rlm.res$Predictor == "cpg_Overlap"] <- "CpG island overlap"
beta.rlm.res$Predictor[beta.rlm.res$Predictor == "cpg_RATIO"] <- "CpG dinucleotide ratio"
beta.rlm.res$Predictor[beta.rlm.res$Predictor == "CGI_SIZE.kb"] <- "CpG island size"
beta.rlm.res$Predictor[beta.rlm.res$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"

# remove the intercept term
beta.rlm.res <- beta.rlm.res[!grepl(beta.rlm.res$Predictor, pattern="Intercept"), ]

beta.rlm.res$Tissue <- "BetaIslet"
beta.rlm.res$Species <- "Human"
beta.rlm.res$Model <- "Multivariate"

write.table(beta.rlm.res,
            file="~/Dropbox/Noise_genomics/Model_results/human_BetaIslet_CGI-multivariateRLM.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

cpg.plot <- ggplot(beta.rlm.res,
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
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-BetaIslet_CGI-multiivariateLM.png",
       height=5.25, width=11.75, dpi=300)

# plot uni and multivariate results together
beta.rlm.all <- do.call(rbind.data.frame, list("Univariate"=beta.rlm.df,
                                               "Multivariate"=beta.rlm.res))

all.plot <- ggplot(beta.rlm.all,
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
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-BetaIslet_CGI-allLM.png",
       height=5.6, width=7.75, dpi=300)
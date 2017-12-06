#source("~/Dropbox/R_sessions/Noise/mESC_chromHMM.R")
source("~/Dropbox/R_sessions/Noise/human_genomic_noise_features.R")
source("~/Dropbox/R_sessions/Noise/human_pancreas_genomic_noise_features.R")

library(scales)
library(ggplot2)
library(MASS)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

human.genomic.features$CGI_SIZE.kb <- human.genomic.features$CGI_SIZE/1000

panc.vars <- colnames(human.genomic.features)
panc.vars <- panc.vars[grepl(panc.vars, pattern="(CGI_SIZE.kb)|(cpg_[GC|O|R])|(TOTLEN)
                             |(SP1)")]

panc.match <- merge(panc.gene.summary, human.genomic.features,
                    by='GENE')
panc.match$Recip.Mean <- 1/panc.match$Mean

# select only those genes with a CpG island
#panc.match <- panc.match[panc.match$N_CpG == 1, ]

#########################################
## univariate robust linear regression ##
#########################################
# this is plotting for presentation/manuscript figures
panc.univariate_list <- list()

for(x in seq_along(panc.vars)){
  .variable <- paste(c("Recip.Mean", panc.vars[x]), collapse=" + ")
  .glm.form <- as.formula(paste("Residual.CV2", .variable, sep=" ~ "))
  
  m.rlm <- rlm(.glm.form, data=panc.match)
  m.robust <- summary(m.rlm)
  m.rlm.res <- as.data.frame(m.robust$coefficients)
  m.rlm.res$Pval <- 2*pt(-abs(m.rlm.res[, 3]), df=3)
  m.rlm.res$Sig <- as.numeric(m.rlm.res$Pval <= 0.05)
  m.res.mat <- as.matrix(m.rlm.res)
  m.res.var <- m.res.mat[3, ]
  names(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
  panc.univariate_list[[panc.vars[x]]] <- as.list(m.res.var)
}

# add the mean expression on it's own
.glm.form <- as.formula(paste("Residual.CV2", "Recip.Mean", sep=" ~ "))

m.rlm <- rlm(.glm.form, data=panc.match)
m.robust <- summary(m.rlm)
m.rlm.res <- as.data.frame(m.robust$coefficients)
m.rlm.res$Pval <- 2*pt(-abs(m.rlm.res[, 3]), df=3)
m.rlm.res$Sig <- as.numeric(m.rlm.res$Pval <= 0.05)
m.res.mat <- as.matrix(m.rlm.res)
m.res.var <- m.res.mat[2, ]
names(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
panc.univariate_list[["Recip.Mean"]] <- as.list(m.res.var)


panc.rlm.df <- do.call(rbind.data.frame, panc.univariate_list)
panc.rlm.df$Predictor <- rownames(panc.rlm.df)
panc.rlm.df$Direction <- "NoEffect"
panc.rlm.df$Direction[panc.rlm.df$COEFF < 0 & panc.rlm.df$Sig == 1] <- "Less"
panc.rlm.df$Direction[panc.rlm.df$COEFF > 0 & panc.rlm.df$Sig == 1] <- "More"

# give the features more informative/better formated names
panc.rlm.df$Predictor[panc.rlm.df$Predictor == "Recip.Mean"] <- "Mean expression"
panc.rlm.df$Predictor[panc.rlm.df$Predictor == "N_CpG"] <- "CpG island"
panc.rlm.df$Predictor[panc.rlm.df$Predictor == "cpg_GCNUM"] <- "CpG island GC content"
panc.rlm.df$Predictor[panc.rlm.df$Predictor == "SP1"] <- "Number SP1 motifs"
panc.rlm.df$Predictor[panc.rlm.df$Predictor == "cpg_Overlap"] <- "CpG island overlap"
panc.rlm.df$Predictor[panc.rlm.df$Predictor == "cpg_RATIO"] <- "CpG dinucleotide ratio"
panc.rlm.df$Predictor[panc.rlm.df$Predictor == "CGI_SIZE.kb"] <- "CpG island size"
panc.rlm.df$Predictor[panc.rlm.df$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"

panc.rlm.df$Tissue <- "AlphaIslet"
panc.rlm.df$Species <- "Human"
panc.rlm.df$Model <- "Univariate"

effect.cols <- c("#62148f", "#878787", "#feaf10")
names(effect.cols) <- c("Less", "NoEffect", "More")

write.table(panc.rlm.df,
            file="~/Dropbox/Noise_genomics/Model_results/human_AlphaIslet_CGI-univariateRLM.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

cpg.plot <- ggplot(panc.rlm.df,
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
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-AlphaIslet_CGI-univariateLM.png",
       height=5.25, width=11.75, dpi=300)

###########################################
## multivariate robust linear regression ##
###########################################
panc.genomic.vars <- paste(c("Recip.Mean",
                             panc.vars),
                           collapse=" + ")

panc.glm.form <- as.formula(paste("Residual.CV2",
                                  panc.genomic.vars, sep=" ~ "))

panc.rlm <- rlm(panc.glm.form, data=panc.match)
panc.robust <- summary(panc.rlm)
panc.rlm.res <- as.data.frame(panc.robust$coefficients)
panc.rlm.res$Pval <- 2*pt(-abs(panc.rlm.res[, 3]), df=dim(panc.match)[2]-1)
panc.rlm.res$Sig <- as.numeric(panc.rlm.res$Pval <= 0.05)
panc.rlm.res$Predictor <- rownames(panc.rlm.res)
colnames(panc.rlm.res) <- c("COEFF", "SE", "STAT", "P", "Sig", "Predictor")
panc.rlm.res$Direction <- "NoEffect"
panc.rlm.res$Direction[panc.rlm.res$COEFF < 0 & panc.rlm.res$Sig == 1] <- "Less"
panc.rlm.res$Direction[panc.rlm.res$COEFF > 0 & panc.rlm.res$Sig == 1] <- "More"

# give the features more informative/better formated names
panc.rlm.res$Predictor[panc.rlm.res$Predictor == "Recip.Mean"] <- "Mean expression"
panc.rlm.res$Predictor[panc.rlm.res$Predictor == "N_CpG"] <- "CpG island"
panc.rlm.res$Predictor[panc.rlm.res$Predictor == "cpg_GCNUM"] <- "CpG island GC content"
panc.rlm.res$Predictor[panc.rlm.res$Predictor == "SP1"] <- "Number SP1 motifs"
panc.rlm.res$Predictor[panc.rlm.res$Predictor == "cpg_Overlap"] <- "CpG island overlap"
panc.rlm.res$Predictor[panc.rlm.res$Predictor == "cpg_RATIO"] <- "CpG dinucleotide ratio"
panc.rlm.res$Predictor[panc.rlm.res$Predictor == "CGI_SIZE.kb"] <- "CpG island size"
panc.rlm.res$Predictor[panc.rlm.res$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"

# remove the intercept term
panc.rlm.res <- panc.rlm.res[!grepl(panc.rlm.res$Predictor, pattern="Intercept"), ]

panc.rlm.res$Tissue <- "AlphaIslet"
panc.rlm.res$Species <- "Human"
panc.rlm.res$Model <- "Multivariate"

write.table(panc.rlm.res,
            file="~/Dropbox/Noise_genomics/Model_results/human_AlphaIslet_CGI-multivariateRLM.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

cpg.plot <- ggplot(panc.rlm.res,
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
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-AlphaIslet_CGI-multiivariateLM.png",
       height=5.25, width=11.75, dpi=300)

# plot uni and multivariate results together
panc.rlm.all <- do.call(rbind.data.frame, list("Univariate"=panc.rlm.df,
                                               "Multivariate"=panc.rlm.res))

all.plot <- ggplot(panc.rlm.all,
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
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-AlphaIslet_CGI-allLM.png",
       height=5.6, width=7.75, dpi=300)
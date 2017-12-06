source("~/Dropbox/R_sessions/Noise/mESC_chromHMM.R")
source("~/Dropbox/R_sessions/Noise/genomic_noise_features.R")
source("~/Dropbox/R_sessions/Noise/mouse_liver_noise_features.R")

library(ggplot2)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")
genomic.features$CGI_SIZE.kb <- genomic.features$CGI_SIZE/1000

liver.vars <- colnames(genomic.features)
liver.vars <- liver.vars[grepl(liver.vars, pattern="(CGI_SIZE.kb)|(cpg_[GC|O|R])|(TOTLEN)
                             |(SP1)")]

liver.match <- merge(liver.gene.summary, genomic.features,
                    by='GENE')
liver.match$Recip.Mean <- 1/liver.match$Mean

# select only those genes with a CpG island
#liver.match <- liver.match[liver.match$N_CpG == 1, ]

#########################################
## univariate robust linear regression ##
#########################################
# this is plotting for presentation/manuscript figures
liver.univariate_list <- list()

for(x in seq_along(liver.vars)){
  .variable <- paste(c("Recip.Mean", liver.vars[x]), collapse=" + ")
  .glm.form <- as.formula(paste("Residual.CV2", .variable, sep=" ~ "))
  
  m.rlm <- rlm(.glm.form, data=liver.match)
  m.robust <- summary(m.rlm)
  m.rlm.res <- as.data.frame(m.robust$coefficients)
  m.rlm.res$Pval <- 2*pt(-abs(m.rlm.res[, 3]), df=3)
  m.rlm.res$Sig <- as.numeric(m.rlm.res$Pval <= 0.05)
  m.res.mat <- as.matrix(m.rlm.res)
  m.res.var <- m.res.mat[3, ]
  names(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
  liver.univariate_list[[liver.vars[x]]] <- as.list(m.res.var)
}

# add the mean expression on it's own
.glm.form <- as.formula(paste("Residual.CV2", "Recip.Mean", sep=" ~ "))

m.rlm <- rlm(.glm.form, data=liver.match)
m.robust <- summary(m.rlm)
m.rlm.res <- as.data.frame(m.robust$coefficients)
m.rlm.res$Pval <- 2*pt(-abs(m.rlm.res[, 3]), df=3)
m.rlm.res$Sig <- as.numeric(m.rlm.res$Pval <= 0.05)
m.res.mat <- as.matrix(m.rlm.res)
m.res.var <- m.res.mat[2, ]
names(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
liver.univariate_list[["Recip.Mean"]] <- as.list(m.res.var)


liver.rlm.df <- do.call(rbind.data.frame, liver.univariate_list)
liver.rlm.df$Predictor <- rownames(liver.rlm.df)
liver.rlm.df$Direction <- "NoEffect"
liver.rlm.df$Direction[liver.rlm.df$COEFF < 0 & liver.rlm.df$Sig == 1] <- "Less"
liver.rlm.df$Direction[liver.rlm.df$COEFF > 0 & liver.rlm.df$Sig == 1] <- "More"

# give the features more informative/better formated names
liver.rlm.df$Predictor[liver.rlm.df$Predictor == "Recip.Mean"] <- "Mean expression"
liver.rlm.df$Predictor[liver.rlm.df$Predictor == "N_CpG"] <- "CpG island"
liver.rlm.df$Predictor[liver.rlm.df$Predictor == "cpg_GCNUM"] <- "CpG island GC content"
liver.rlm.df$Predictor[liver.rlm.df$Predictor == "SP1"] <- "Number SP1 motifs"
liver.rlm.df$Predictor[liver.rlm.df$Predictor == "cpg_Overlap"] <- "CpG island overlap"
liver.rlm.df$Predictor[liver.rlm.df$Predictor == "cpg_RATIO"] <- "CpG dinucleotide ratio"
liver.rlm.df$Predictor[liver.rlm.df$Predictor == "CGI_SIZE.kb"] <- "CpG island size"
liver.rlm.df$Predictor[liver.rlm.df$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"

liver.rlm.df$Tissue <- "liver"
liver.rlm.df$Species <- "Mouse"
liver.rlm.df$Model <- "Univariate"

effect.cols <- c("#62148f", "#878787", "#feaf10")
names(effect.cols) <- c("Less", "NoEffect", "More")

write.table(liver.rlm.df,
            file="~/Dropbox/Noise_genomics/Model_results/mouse_liver_CGI-univariateRLM.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")

cpg.plot <- ggplot(liver.rlm.df,
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
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mouse-liver_CGI-univariateLM.png",
       height=4.75, width=6.75, dpi=300)

###########################################
## multivariate robust linear regression ##
###########################################
liver.genomic.vars <- paste(c("Recip.Mean",
                             liver.vars),
                           collapse=" + ")

liver.glm.form <- as.formula(paste("Residual.CV2",
                                  liver.genomic.vars, sep=" ~ "))

liver.rlm <- rlm(liver.glm.form, data=liver.match)
liver.robust <- summary(liver.rlm)
liver.rlm.res <- as.data.frame(liver.robust$coefficients)
liver.rlm.res$Pval <- 2*pt(-abs(liver.rlm.res[, 3]), df=dim(liver.match)[2]-1)
liver.rlm.res$Sig <- as.numeric(liver.rlm.res$Pval <= 0.05)
liver.rlm.res$Predictor <- rownames(liver.rlm.res)
colnames(liver.rlm.res) <- c("COEFF", "SE", "STAT", "P", "Sig", "Predictor")
liver.rlm.res$Direction <- "NoEffect"
liver.rlm.res$Direction[liver.rlm.res$COEFF < 0 & liver.rlm.res$Sig == 1] <- "Less"
liver.rlm.res$Direction[liver.rlm.res$COEFF > 0 & liver.rlm.res$Sig == 1] <- "More"

# give the features more informative/better formated names
liver.rlm.res$Predictor[liver.rlm.res$Predictor == "Recip.Mean"] <- "Mean expression"
liver.rlm.res$Predictor[liver.rlm.res$Predictor == "N_CpG"] <- "CpG island"
liver.rlm.res$Predictor[liver.rlm.res$Predictor == "cpg_GCNUM"] <- "CpG island GC content"
liver.rlm.res$Predictor[liver.rlm.res$Predictor == "SP1"] <- "Number SP1 motifs"
liver.rlm.res$Predictor[liver.rlm.res$Predictor == "cpg_Overlap"] <- "CpG island overlap"
liver.rlm.res$Predictor[liver.rlm.res$Predictor == "cpg_RATIO"] <- "CpG dinucleotide ratio"
liver.rlm.res$Predictor[liver.rlm.res$Predictor == "CGI_SIZE.kb"] <- "CpG island size"
liver.rlm.res$Predictor[liver.rlm.res$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"

# remove the intercept term
liver.rlm.res <- liver.rlm.res[!grepl(liver.rlm.res$Predictor, pattern="Intercept"), ]

liver.rlm.res$Tissue <- "liver"
liver.rlm.res$Species <- "Mouse"
liver.rlm.res$Model <- "Multivariate"

write.table(liver.rlm.res,
            file="~/Dropbox/Noise_genomics/Model_results/mouse_liver_CGI-multivariateRLM.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")

cpg.plot <- ggplot(liver.rlm.res,
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
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mouse-liver_CGI-multivariateLM.png",
       height=4.75, width=6.75, dpi=300)

# plot uni and multivariate results together
liver.rlm.all <- do.call(rbind.data.frame, list("Univariate"=liver.rlm.df,
                                               "Multivariate"=liver.rlm.res))

all.plot <- ggplot(liver.rlm.all,
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
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mouse-liver_CGI-allLM.png",
       height=5.6, width=7.75, dpi=300)

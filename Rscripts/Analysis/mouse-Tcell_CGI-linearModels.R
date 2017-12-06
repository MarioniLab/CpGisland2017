source("~/Dropbox/R_sessions/Noise/mESC_chromHMM.R")
source("~/Dropbox/R_sessions/Noise/genomic_noise_features.R")
source("~/Dropbox/R_sessions/Noise/tcell_genomic_noise_features.R")
library(scales)
library(ggplot2)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

genomic.features$CGI_SIZE.kb <- genomic.features$CGI_SIZE/1000

tcell.vars <- colnames(genomic.features)
tcell.vars <- tcell.vars[grepl(tcell.vars, pattern="(CGI_SIZE.kb)|(cpg_[GC|O|R])|(TOTLEN)
                             |(SP1)")]

tcell.match <- merge(tcell.gene.summary, genomic.features,
                     by='GENE')
tcell.match$Recip.Mean <- 1/tcell.match$Mean

# select only those genes with a CpG island
#tcell.match <- tcell.match[tcell.match$N_CpG == 1, ]

#########################################
## univariate robust linear regression ##
#########################################
# this is plotting for presentation/manuscript figures
tcell.univariate_list <- list()

for(x in seq_along(tcell.vars)){
  .variable <- paste(c("Recip.Mean", tcell.vars[x]), collapse=" + ")
  .glm.form <- as.formula(paste("Residual.CV2", .variable, sep=" ~ "))
  
  m.rlm <- rlm(.glm.form, data=tcell.match)
  m.robust <- summary(m.rlm)
  m.rlm.res <- as.data.frame(m.robust$coefficients)
  m.rlm.res$Pval <- 2*pt(-abs(m.rlm.res[, 3]), df=3)
  m.rlm.res$Sig <- as.numeric(m.rlm.res$Pval <= 0.05)
  m.res.mat <- as.matrix(m.rlm.res)
  m.res.var <- m.res.mat[3, ]
  names(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
  tcell.univariate_list[[tcell.vars[x]]] <- as.list(m.res.var)
}

# add the mean expression on it's own
.glm.form <- as.formula(paste("Residual.CV2", "Recip.Mean", sep=" ~ "))

m.rlm <- rlm(.glm.form, data=tcell.match)
m.robust <- summary(m.rlm)
m.rlm.res <- as.data.frame(m.robust$coefficients)
m.rlm.res$Pval <- 2*pt(-abs(m.rlm.res[, 3]), df=3)
m.rlm.res$Sig <- as.numeric(m.rlm.res$Pval <= 0.05)
m.res.mat <- as.matrix(m.rlm.res)
m.res.var <- m.res.mat[2, ]
names(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
tcell.univariate_list[["Recip.Mean"]] <- as.list(m.res.var)


tcell.rlm.df <- do.call(rbind.data.frame, tcell.univariate_list)
tcell.rlm.df$Predictor <- rownames(tcell.rlm.df)
tcell.rlm.df$Direction <- "NoEffect"
tcell.rlm.df$Direction[tcell.rlm.df$COEFF < 0 & tcell.rlm.df$Sig == 1] <- "Less"
tcell.rlm.df$Direction[tcell.rlm.df$COEFF > 0 & tcell.rlm.df$Sig == 1] <- "More"

# give the features more informative/better formated names
tcell.rlm.df$Predictor[tcell.rlm.df$Predictor == "Recip.Mean"] <- "Mean expression"
tcell.rlm.df$Predictor[tcell.rlm.df$Predictor == "N_CpG"] <- "CpG island"
tcell.rlm.df$Predictor[tcell.rlm.df$Predictor == "cpg_GCNUM"] <- "CpG island GC content"
tcell.rlm.df$Predictor[tcell.rlm.df$Predictor == "SP1"] <- "Number SP1 motifs"
tcell.rlm.df$Predictor[tcell.rlm.df$Predictor == "cpg_Overlap"] <- "CpG island overlap"
tcell.rlm.df$Predictor[tcell.rlm.df$Predictor == "cpg_RATIO"] <- "CpG dinucleotide ratio"
tcell.rlm.df$Predictor[tcell.rlm.df$Predictor == "CGI_SIZE.kb"] <- "CpG island size"
tcell.rlm.df$Predictor[tcell.rlm.df$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"

tcell.rlm.df$Tissue <- "Tcell"
tcell.rlm.df$Species <- "Mouse"
tcell.rlm.df$Model <- "Univariate"

effect.cols <- c("#62148f", "#878787", "#feaf10")
names(effect.cols) <- c("Less", "NoEffect", "More")

write.table(tcell.rlm.df,
            file="~/Dropbox/Noise_genomics/Model_results/mouse_Tcell_CGI-univariateRLM.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")

cpg.plot <- ggplot(tcell.rlm.df,
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
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/Tcell_CGI-univariateLM.png",
       height=4.75, width=6.75, dpi=300)

###########################################
## multivariate robust linear regression ##
###########################################
tcell.genomic.vars <- paste(c("Recip.Mean",
                              tcell.vars),
                            collapse=" + ")

tcell.glm.form <- as.formula(paste("Residual.CV2",
                                   tcell.genomic.vars, sep=" ~ "))

tcell.rlm <- rlm(tcell.glm.form, data=tcell.match)
tcell.robust <- summary(tcell.rlm)
tcell.rlm.res <- as.data.frame(tcell.robust$coefficients)
tcell.rlm.res$Pval <- 2*pt(-abs(tcell.rlm.res[, 3]), df=dim(tcell.match)[2]-1)
tcell.rlm.res$Sig <- as.numeric(tcell.rlm.res$Pval <= 0.05)
tcell.rlm.res$Predictor <- rownames(tcell.rlm.res)
colnames(tcell.rlm.res) <- c("COEFF", "SE", "STAT", "P", "Sig", "Predictor")
tcell.rlm.res$Direction <- "NoEffect"
tcell.rlm.res$Direction[tcell.rlm.res$COEFF < 0 & tcell.rlm.res$Sig == 1] <- "Less"
tcell.rlm.res$Direction[tcell.rlm.res$COEFF > 0 & tcell.rlm.res$Sig == 1] <- "More"

# give the features more informative/better formated names
tcell.rlm.res$Predictor[tcell.rlm.res$Predictor == "Recip.Mean"] <- "Mean expression"
tcell.rlm.res$Predictor[tcell.rlm.res$Predictor == "N_CpG"] <- "CpG island"
tcell.rlm.res$Predictor[tcell.rlm.res$Predictor == "cpg_GCNUM"] <- "CpG island GC content"
tcell.rlm.res$Predictor[tcell.rlm.res$Predictor == "SP1"] <- "Number SP1 motifs"
tcell.rlm.res$Predictor[tcell.rlm.res$Predictor == "cpg_Overlap"] <- "CpG island overlap"
tcell.rlm.res$Predictor[tcell.rlm.res$Predictor == "cpg_RATIO"] <- "CpG dinucleotide ratio"
tcell.rlm.res$Predictor[tcell.rlm.res$Predictor == "CGI_SIZE.kb"] <- "CpG island size"
tcell.rlm.res$Predictor[tcell.rlm.res$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"

# remove the intercept term
tcell.rlm.res <- tcell.rlm.res[!grepl(tcell.rlm.res$Predictor, pattern="Intercept"), ]

tcell.rlm.res$Tissue <- "Tcell"
tcell.rlm.res$Species <- "Mouse"
tcell.rlm.res$Model <- "Multivariate"

write.table(tcell.rlm.res,
            file="~/Dropbox/Noise_genomics/Model_results/mouse_Tcell_CGI-multivariateRLM.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")

cpg.plot <- ggplot(tcell.rlm.res,
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
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/Tcell_CGI-multiivariateLM.png",
       height=4.75, width=6.75, dpi=300)

# plot uni and multivariate results together
tcell.rlm.all <- do.call(rbind.data.frame, list("Univariate"=tcell.rlm.df,
                                                "Multivariate"=tcell.rlm.res))

all.plot <- ggplot(tcell.rlm.all,
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
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/Tcell_CGI-allLM.png",
       height=5.6, width=7.75, dpi=300)


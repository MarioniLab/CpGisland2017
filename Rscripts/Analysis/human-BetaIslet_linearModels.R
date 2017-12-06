source("~/Dropbox/R_sessions/Noise/mESC_chromHMM.R")
source("~/Dropbox/R_sessions/Noise/human_genomic_noise_features.R")
source("~/Dropbox/R_sessions/Noise/human_beta_islet_genomic_noise_features.R")

library(ggplot2)
library(MASS)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

beta.vars <- colnames(genomic.features)
beta.vars <- beta.vars[!grepl(beta.vars, pattern="(NMI)|(CGI_SIZE)|(cpg_)")]

beta.genomic.vars <- paste(c("Mean",
                             beta.vars),
                           collapse=" + ")

beta.match <- merge(beta.gene.summary, genomic.features,
                    by='GENE')

#########################################
## univariate robust linear regression ##
#########################################
# this is plotting for presentation/manuscript figures
beta.univariate_list <- list()

beta.var.names <- unlist(strsplit(beta.genomic.vars, split=" + ", fixed=T))
# remove redundant variables, i.e CpG island AND NMI
# remove CpG island characteristics
beta.var.names <- beta.var.names[!grepl(beta.var.names, pattern="(NMI)|(CGI_SIZE)|(cpg_)|(CGI)|(GENE)")]

for(x in seq_along(beta.var.names)){
  .variable <- paste(beta.var.names[x], sep=" + ")
  .glm.form <- as.formula(paste("Residual.CV2", .variable, sep=" ~ "))
  
  m.rlm <- rlm(.glm.form, data=beta.match)
  m.robust <- summary(m.rlm)
  m.rlm.res <- as.data.frame(m.robust$coefficients)
  m.rlm.res$Pval <- 2*pt(-abs(m.rlm.res[, 3]), df=3)
  m.rlm.res$Sig <- as.numeric(m.rlm.res$Pval <= 0.05)
  m.res.mat <- as.matrix(m.rlm.res)
  m.res.var <- m.res.mat[2, ]
  names(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
  beta.univariate_list[[beta.var.names[x]]] <- as.list(m.res.var)
}

beta.rlm.df <- do.call(rbind.data.frame, beta.univariate_list)
beta.rlm.df$Predictor <- rownames(beta.rlm.df)
beta.rlm.df$Direction <- "NoEffect"
beta.rlm.df$Direction[beta.rlm.df$COEFF < 0 & beta.rlm.df$Sig == 1] <- "Less"
beta.rlm.df$Direction[beta.rlm.df$COEFF > 0 & beta.rlm.df$Sig == 1] <- "More"

# give the features more informative/better formated names
beta.rlm.df$Predictor[beta.rlm.df$Predictor == "Mean"] <- "Mean expression"
beta.rlm.df$Predictor[beta.rlm.df$Predictor == "N_CpG"] <- "CpG island"
beta.rlm.df$Predictor[beta.rlm.df$Predictor == "GC"] <- "Promoter GC %"
beta.rlm.df$Predictor[beta.rlm.df$Predictor == "SP1"] <- "Number SP1 motifs"
beta.rlm.df$Predictor[beta.rlm.df$Predictor == "TBP"] <- "Number TBP motifs"
beta.rlm.df$Predictor[beta.rlm.df$Predictor == "MED_PHAST"] <- "Median PhastCons"
beta.rlm.df$Predictor[beta.rlm.df$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
beta.rlm.df$Predictor[beta.rlm.df$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
beta.rlm.df$Predictor[beta.rlm.df$Predictor == "NALIGN_PHAST"] <- "Numer aligned bases"
beta.rlm.df$Predictor[beta.rlm.df$Predictor == "EXON_AVLENGTH"] <- "Mean exon length"
beta.rlm.df$Predictor[beta.rlm.df$Predictor == "EXON_COUNT"] <- "Number of exons"
beta.rlm.df$Predictor[beta.rlm.df$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"
beta.rlm.df$Predictor[beta.rlm.df$Predictor == "EXON_VARLENGTH"] <- "Exon length variance"

beta.rlm.df$Tissue <- "BetaIslet"
beta.rlm.df$Species <- "Human"

univar.plot <- ggplot(beta.rlm.df,
                      aes(x=reorder(Predictor, -STAT),
                          y=STAT, fill=Direction)) +
  geom_point(alpha=0.55, shape=21, size=5) + theme_mike() +
  scale_fill_manual(values=c("#62148f", "#feaf10", "#878787")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Annotation", y="t-statistic") +
  guides(fill=FALSE) +
  geom_hline(mapping=aes(yintercept=0), linetype="dashed", colour="grey") 

ggsave(univar.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-BetaIslet_univariateLM.png",
       height=4.75, width=6.75, dpi=300)


###########################################
## multivariate robust linear regression ##
###########################################
beta.genomic.vars <- paste(c("Mean",
                             beta.var.names),
                           collapse=" + ")

beta.glm.form <- as.formula(paste("Residual.CV2",
                                  beta.genomic.vars, sep=" ~ "))

beta.rlm <- rlm(beta.glm.form, data=beta.match)
beta.robust <- summary(beta.rlm)
beta.rlm.res <- as.data.frame(beta.robust$coefficients)
beta.rlm.res$Pval <- 2*pt(-abs(beta.rlm.res[, 3]), df=dim(beta.match)[2]-1)
beta.rlm.res$Sig <- as.numeric(beta.rlm.res$Pval <= 0.05)
beta.rlm.res$Predictor <- rownames(beta.rlm.res)
colnames(beta.rlm.res) <- c("COEFF", "SE", "STAT", "P", "Sig", "Predictor")
beta.rlm.res$Direction <- "NoEffect"
beta.rlm.res$Direction[beta.rlm.res$COEFF < 0 & beta.rlm.res$Sig == 1] <- "Less"
beta.rlm.res$Direction[beta.rlm.res$COEFF > 0 & beta.rlm.res$Sig == 1] <- "More"


# give the features more informative/better formated names
beta.rlm.res$Predictor[beta.rlm.res$Predictor == "N_CpG"] <- "CpG island"
beta.rlm.res$Predictor[beta.rlm.res$Predictor == "GC"] <- "Promoter GC %"
beta.rlm.res$Predictor[beta.rlm.res$Predictor == "SP1"] <- "Number SP1 motifs"
beta.rlm.res$Predictor[beta.rlm.res$Predictor == "TBP"] <- "Number TBP motifs"
beta.rlm.res$Predictor[beta.rlm.res$Predictor == "MED_PHAST"] <- "Median PhastCons"
beta.rlm.res$Predictor[beta.rlm.res$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
beta.rlm.res$Predictor[beta.rlm.res$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
beta.rlm.res$Predictor[beta.rlm.res$Predictor == "NALIGN_PHAST"] <- "Numer aligned bases"
beta.rlm.res$Predictor[beta.rlm.res$Predictor == "EXON_AVLENGTH"] <- "Mean exon length"
beta.rlm.res$Predictor[beta.rlm.res$Predictor == "EXON_COUNT"] <- "Number of exons"
beta.rlm.res$Predictor[beta.rlm.res$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"
beta.rlm.res$Predictor[beta.rlm.res$Predictor == "EXON_VARLENGTH"] <- "Exon length variance"
beta.rlm.res$Predictor[beta.rlm.res$Predictor == "Mean"] <- "Mean expression"

beta.rlm.res <- beta.rlm.res[!grepl(beta.rlm.res$Predictor, pattern="Intercept"), ]
beta.rlm.res$Tissue <- "BetaIslet"
beta.rlm.res$Species <- "Human"

write.table(beta.rlm.res,
            file="~/Dropbox/Noise_genomics/Model_results/human_BetaIslet_multivariateRLM.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")

multivar.plot <- ggplot(beta.rlm.res,
                        aes(x=reorder(Predictor, -STAT),
                            y=STAT, fill=Direction)) +
  geom_point(alpha=0.55, shape=21, size=5) + theme_mike() +
  scale_fill_manual(values=c("#62148f", "#feaf10", "#878787")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Annotation", y="t-statistic") +
  guides(fill=FALSE) +
  geom_hline(mapping=aes(yintercept=0), linetype="dashed", colour="grey")

ggsave(multivar.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-BetaIslet_multiivariateLM.png",
       height=4.75, width=6.75, dpi=300)

# plot the uni- and multivariate on the same graph?
beta.rlm.df$Model <- "Univariate"
beta.rlm.res$Model <- "Multivariate"

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
  guides(fill=FALSE, shape=FALSE) 

ggsave(all.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-BetaIslet_allLM.png",
       height=5.6, width=7.75, dpi=300)

###############################
## Plot CpG islands vs alpha ##
###############################
beta.match$CpGisland <- factor(beta.match$N_CpG,
                               levels=c(0, 1),
                               labels=c("Non-CpG island", "CpG island"))
cpg.cols <- c("#027E00", "#00DDEC")
names(cpg.cols) <- levels(beta.match$CpGisland)

tc_cpg <- ggplot(beta.match, aes(x=CpGisland, y=Residual.CV2, colour=CpGisland)) + 
  theme_mike() + 
  geom_jitter(position=position_jitterdodge(jitter.height=0,
                                            jitter.width=1.5),
              alpha=0.7) +
  geom_boxplot(width=0.5, fill='white', colour='black') +
  scale_colour_manual(values=cpg.cols) +
  labs(x="Overlapping CpG island", y=expression(paste(alpha["r"], " Overdispersion"))) +
  guides(colour=FALSE)

ggsave(tc_cpg,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-BetaIslet_boxplot_CpGislands-overdispersion.png",
       height=3.75, width=4.75, dpi=300)

tc_cpgMe <- ggplot(beta.match, aes(x=CpGisland, y=Mean, colour=CpGisland)) + 
  theme_mike() + 
  geom_jitter(position=position_jitterdodge(jitter.height=0,
                                            jitter.width=1.5),
              alpha=0.7) +
  geom_boxplot(width=0.5, fill='white', colour='black') +
  scale_colour_manual(values=cpg.cols) +
  labs(x="Overlapping CpG island", y=expression(paste("Mean log"[2], " Expression"))) +
  guides(colour=FALSE)

ggsave(tc_cpgMe,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-BetaIslet_boxplot_CpGislands-mean.png",
       height=3.75, width=4.75, dpi=300)
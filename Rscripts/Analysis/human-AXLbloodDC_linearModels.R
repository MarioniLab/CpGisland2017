source("~/Dropbox/R_sessions/Noise/human_genomic_noise_features.R")
source("~/Dropbox/R_sessions/Noise/human_bloodDC_AXL_genomic_noise_features.R")

library(ggplot2)
library(MASS)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

axl.vars <- colnames(genomic.features)
axl.vars <- axl.vars[!grepl(axl.vars, pattern="(NMI)|(CGI_SIZE)|(cpg_)")]

axl.genomic.vars <- paste(c("Mean",
                             axl.vars),
                           collapse=" + ")

axl.match <- merge(axl.gene.summary, genomic.features,
                    by='GENE')

#########################################
## univariate robust linear regression ##
#########################################
# this is plotting for presentation/manuscript figures
axl.univariate_list <- list()

axl.var.names <- unlist(strsplit(axl.genomic.vars, split=" + ", fixed=T))
# remove redundant variables, i.e CpG island AND NMI
# remove CpG island characteristics
axl.var.names <- axl.var.names[!grepl(axl.var.names, pattern="(NMI)|(CGI_SIZE)|(cpg_)|(CGI)|(GENE)")]

for(x in seq_along(axl.var.names)){
  .variable <- paste(axl.var.names[x], sep=" + ")
  .glm.form <- as.formula(paste("Alpha_r", .variable, sep=" ~ "))
  
  m.rlm <- rlm(.glm.form, data=axl.match)
  m.robust <- summary(m.rlm)
  m.rlm.res <- as.data.frame(m.robust$coefficients)
  m.rlm.res$Pval <- 2*pt(-abs(m.rlm.res[, 3]), df=3)
  m.rlm.res$Sig <- as.numeric(m.rlm.res$Pval <= 0.05)
  m.res.mat <- as.matrix(m.rlm.res)
  m.res.var <- m.res.mat[2, ]
  names(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
  axl.univariate_list[[axl.var.names[x]]] <- as.list(m.res.var)
}

axl.rlm.df <- do.call(rbind.data.frame, axl.univariate_list)
axl.rlm.df$Predictor <- rownames(axl.rlm.df)
axl.rlm.df$Direction <- "NoEffect"
axl.rlm.df$Direction[axl.rlm.df$COEFF < 0 & axl.rlm.df$Sig == 1] <- "Less"
axl.rlm.df$Direction[axl.rlm.df$COEFF > 0 & axl.rlm.df$Sig == 1] <- "More"

# give the features more informative/better formated names
axl.rlm.df$Predictor[axl.rlm.df$Predictor == "Mean"] <- "Mean expression"
axl.rlm.df$Predictor[axl.rlm.df$Predictor == "N_CpG"] <- "CpG island"
axl.rlm.df$Predictor[axl.rlm.df$Predictor == "GC"] <- "Promoter GC %"
axl.rlm.df$Predictor[axl.rlm.df$Predictor == "SP1"] <- "Number SP1 motifs"
axl.rlm.df$Predictor[axl.rlm.df$Predictor == "TBP"] <- "Number TBP motifs"
axl.rlm.df$Predictor[axl.rlm.df$Predictor == "MED_PHAST"] <- "Median PhastCons"
axl.rlm.df$Predictor[axl.rlm.df$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
axl.rlm.df$Predictor[axl.rlm.df$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
axl.rlm.df$Predictor[axl.rlm.df$Predictor == "NALIGN_PHAST"] <- "Numer aligned bases"
axl.rlm.df$Predictor[axl.rlm.df$Predictor == "EXON_AVLENGTH"] <- "Mean exon length"
axl.rlm.df$Predictor[axl.rlm.df$Predictor == "EXON_COUNT"] <- "Number of exons"
axl.rlm.df$Predictor[axl.rlm.df$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"
axl.rlm.df$Predictor[axl.rlm.df$Predictor == "EXON_VARLENGTH"] <- "Exon length variance"

univar.plot <- ggplot(axl.rlm.df,
                      aes(x=reorder(Predictor, -STAT),
                          y=STAT, fill=Direction)) +
  geom_point(alpha=0.55, shape=21, size=5) + theme_mike() +
  scale_fill_manual(values=c("#62148f", "#feaf10", "#878787")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Annotation", y="t-statistic") +
  guides(fill=FALSE) +
  geom_hline(mapping=aes(yintercept=0), linetype="dashed", colour="grey") +
  scale_y_continuous(limits=c(-50, 50))

ggsave(univar.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-AXLbloodDC_univariateLM.png",
       height=4.75, width=6.75, dpi=300)


###########################################
## multivariate robust linear regression ##
###########################################
axl.genomic.vars <- paste(c("Mean",
                             axl.var.names),
                           collapse=" + ")

axl.glm.form <- as.formula(paste("Alpha_r",
                                  axl.genomic.vars, sep=" ~ "))

axl.rlm <- rlm(axl.glm.form, data=axl.match)
axl.robust <- summary(axl.rlm)
axl.rlm.res <- as.data.frame(axl.robust$coefficients)
axl.rlm.res$Pval <- 2*pt(-abs(axl.rlm.res[, 3]), df=dim(axl.match)[2]-1)
axl.rlm.res$Sig <- as.numeric(axl.rlm.res$Pval <= 0.05)
axl.rlm.res$Predictor <- rownames(axl.rlm.res)
colnames(axl.rlm.res) <- c("COEFF", "SE", "STAT", "P", "Sig", "Predictor")
axl.rlm.res$Direction <- "NoEffect"
axl.rlm.res$Direction[axl.rlm.res$COEFF < 0 & axl.rlm.res$Sig == 1] <- "Less"
axl.rlm.res$Direction[axl.rlm.res$COEFF > 0 & axl.rlm.res$Sig == 1] <- "More"


# give the features more informative/better formated names
axl.rlm.res$Predictor[axl.rlm.res$Predictor == "N_CpG"] <- "CpG island"
axl.rlm.res$Predictor[axl.rlm.res$Predictor == "GC"] <- "Promoter GC %"
axl.rlm.res$Predictor[axl.rlm.res$Predictor == "SP1"] <- "Number SP1 motifs"
axl.rlm.res$Predictor[axl.rlm.res$Predictor == "TBP"] <- "Number TBP motifs"
axl.rlm.res$Predictor[axl.rlm.res$Predictor == "MED_PHAST"] <- "Median PhastCons"
axl.rlm.res$Predictor[axl.rlm.res$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
axl.rlm.res$Predictor[axl.rlm.res$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
axl.rlm.res$Predictor[axl.rlm.res$Predictor == "NALIGN_PHAST"] <- "Numer aligned bases"
axl.rlm.res$Predictor[axl.rlm.res$Predictor == "EXON_AVLENGTH"] <- "Mean exon length"
axl.rlm.res$Predictor[axl.rlm.res$Predictor == "EXON_COUNT"] <- "Number of exons"
axl.rlm.res$Predictor[axl.rlm.res$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"
axl.rlm.res$Predictor[axl.rlm.res$Predictor == "EXON_VARLENGTH"] <- "Exon length variance"
axl.rlm.res$Predictor[axl.rlm.res$Predictor == "Mean"] <- "Mean expression"

axl.rlm.res <- axl.rlm.res[!grepl(axl.rlm.res$Predictor, pattern="Intercept"), ]


multivar.plot <- ggplot(axl.rlm.res,
                        aes(x=reorder(Predictor, -STAT),
                            y=STAT, fill=Direction)) +
  geom_point(alpha=0.55, shape=21, size=5) + theme_mike() +
  scale_fill_manual(values=c("#62148f", "#feaf10", "#878787")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Annotation", y="t-statistic") +
  guides(fill=FALSE) +
  geom_hline(mapping=aes(yintercept=0), linetype="dashed", colour="grey") +
  scale_y_continuous(limits=c(-50, 50))

ggsave(multivar.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-AXLbloodDC_multiivariateLM.png",
       height=4.75, width=6.75, dpi=300)

# plot the uni- and multivariate on the same graph?
axl.rlm.df$Model <- "Univariate"
axl.rlm.res$Model <- "Multivariate"

axl.rlm.all <- do.call(rbind.data.frame, list("Univariate"=axl.rlm.df,
                                               "Multivariate"=axl.rlm.res))

all.plot <- ggplot(axl.rlm.all,
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
  scale_y_continuous(limits=c(-50, 50))

ggsave(all.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-AXLbloodDC_allLM.png",
       height=5.6, width=7.75, dpi=300)

###############################
## Plot CpG islands vs alpha ##
###############################
axl.match$CpGisland <- factor(axl.match$N_CpG,
                               levels=c(0, 1),
                               labels=c("Non-CpG island", "CpG island"))
cpg.cols <- c("#027E00", "#00DDEC")
names(cpg.cols) <- levels(axl.match$CpGisland)

tc_cpg <- ggplot(axl.match, aes(x=CpGisland, y=Alpha_r, colour=CpGisland)) + 
  theme_mike() + 
  geom_jitter(position=position_jitterdodge(jitter.height=0,
                                            jitter.width=1.5),
              alpha=0.7) +
  geom_boxplot(width=0.5, fill='white', colour='black') +
  scale_colour_manual(values=cpg.cols) +
  labs(x="Overlapping CpG island", y=expression(paste(alpha["r"], " Overdispersion"))) +
  guides(colour=FALSE)

ggsave(tc_cpg,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-AXLbloodDC_boxplot_CpGislands-overdispersion.png",
       height=3.75, width=4.75, dpi=300)

tc_cpgMe <- ggplot(axl.match, aes(x=CpGisland, y=Mean, colour=CpGisland)) + 
  theme_mike() + 
  geom_jitter(position=position_jitterdodge(jitter.height=0,
                                            jitter.width=1.5),
              alpha=0.7) +
  geom_boxplot(width=0.5, fill='white', colour='black') +
  scale_colour_manual(values=cpg.cols) +
  labs(x="Overlapping CpG island", y=expression(paste("Mean log"[2], " Expression"))) +
  guides(colour=FALSE)

ggsave(tc_cpgMe,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-AXLbloodDC_boxplot_CpGislands-mean.png",
       height=3.75, width=4.75, dpi=300)
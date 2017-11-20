source("~/Dropbox/R_sessions/Noise/mESC_chromHMM.R")
source("~/Dropbox/R_sessions/Noise/human_genomic_noise_features.R")
source("~/Dropbox/R_sessions/Noise/human_monocyte_genomic_noise_features.R")

library(ggplot2)
library(MASS)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

monocyte.vars <- colnames(genomic.features)
monocyte.vars <- monocyte.vars[!grepl(monocyte.vars, pattern="(NMI)|(CGI_SIZE)|(cpg_)")]

monocyte.genomic.vars <- paste(c("Mean",
                              monocyte.vars),
                            collapse=" + ")

monocyte.match <- merge(monocyte.gene.summary, genomic.features,
                     by='GENE')

#########################################
## univariate robust linear regression ##
#########################################
# this is plotting for presentation/manuscript figures
monocyte.univariate_list <- list()

monocyte.var.names <- unlist(strsplit(monocyte.genomic.vars, split=" + ", fixed=T))
# remove redundant variables, i.e CpG island AND NMI
# remove CpG island characteristics
monocyte.var.names <- monocyte.var.names[!grepl(monocyte.var.names, pattern="(NMI)|(CGI_SIZE)|(cpg_)|(CGI)|(GENE)")]

for(x in seq_along(monocyte.var.names)){
  .variable <- paste(monocyte.var.names[x], sep=" + ")
  .glm.form <- as.formula(paste("Alpha_r", .variable, sep=" ~ "))
  
  m.rlm <- rlm(.glm.form, data=monocyte.match)
  m.robust <- summary(m.rlm)
  m.rlm.res <- as.data.frame(m.robust$coefficients)
  m.rlm.res$Pval <- 2*pt(-abs(m.rlm.res[, 3]), df=3)
  m.rlm.res$Sig <- as.numeric(m.rlm.res$Pval <= 0.05)
  m.res.mat <- as.matrix(m.rlm.res)
  m.res.var <- m.res.mat[2, ]
  names(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
  monocyte.univariate_list[[monocyte.var.names[x]]] <- as.list(m.res.var)
}

monocyte.rlm.df <- do.call(rbind.data.frame, monocyte.univariate_list)
monocyte.rlm.df$Predictor <- rownames(monocyte.rlm.df)
monocyte.rlm.df$Direction <- "NoEffect"
monocyte.rlm.df$Direction[monocyte.rlm.df$COEFF < 0 & monocyte.rlm.df$Sig == 1] <- "Less"
monocyte.rlm.df$Direction[monocyte.rlm.df$COEFF > 0 & monocyte.rlm.df$Sig == 1] <- "More"

# give the features more informative/better formated names
monocyte.rlm.df$Predictor[monocyte.rlm.df$Predictor == "Mean"] <- "Mean expression"
monocyte.rlm.df$Predictor[monocyte.rlm.df$Predictor == "N_CpG"] <- "CpG island"
monocyte.rlm.df$Predictor[monocyte.rlm.df$Predictor == "GC"] <- "Promoter GC %"
monocyte.rlm.df$Predictor[monocyte.rlm.df$Predictor == "SP1"] <- "Number SP1 motifs"
monocyte.rlm.df$Predictor[monocyte.rlm.df$Predictor == "TBP"] <- "Number TBP motifs"
monocyte.rlm.df$Predictor[monocyte.rlm.df$Predictor == "MED_PHAST"] <- "Median PhastCons"
monocyte.rlm.df$Predictor[monocyte.rlm.df$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
monocyte.rlm.df$Predictor[monocyte.rlm.df$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
monocyte.rlm.df$Predictor[monocyte.rlm.df$Predictor == "NALIGN_PHAST"] <- "Numer aligned bases"
monocyte.rlm.df$Predictor[monocyte.rlm.df$Predictor == "EXON_AVLENGTH"] <- "Mean exon length"
monocyte.rlm.df$Predictor[monocyte.rlm.df$Predictor == "EXON_COUNT"] <- "Number of exons"
monocyte.rlm.df$Predictor[monocyte.rlm.df$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"
monocyte.rlm.df$Predictor[monocyte.rlm.df$Predictor == "EXON_VARLENGTH"] <- "Exon length variance"

univar.plot <- ggplot(monocyte.rlm.df,
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
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-monocyte_univariateLM.png",
       height=4.75, width=6.75, dpi=300)


###########################################
## multivariate robust linear regression ##
###########################################
monocyte.genomic.vars <- paste(c("Mean",
                              monocyte.var.names),
                            collapse=" + ")

monocyte.glm.form <- as.formula(paste("Alpha_r",
                                   monocyte.genomic.vars, sep=" ~ "))

monocyte.rlm <- rlm(monocyte.glm.form, data=monocyte.match)
monocyte.robust <- summary(monocyte.rlm)
monocyte.rlm.res <- as.data.frame(monocyte.robust$coefficients)
monocyte.rlm.res$Pval <- 2*pt(-abs(monocyte.rlm.res[, 3]), df=dim(monocyte.match)[2]-1)
monocyte.rlm.res$Sig <- as.numeric(monocyte.rlm.res$Pval <= 0.05)
monocyte.rlm.res$Predictor <- rownames(monocyte.rlm.res)
colnames(monocyte.rlm.res) <- c("COEFF", "SE", "STAT", "P", "Sig", "Predictor")
monocyte.rlm.res$Direction <- "NoEffect"
monocyte.rlm.res$Direction[monocyte.rlm.res$COEFF < 0 & monocyte.rlm.res$Sig == 1] <- "Less"
monocyte.rlm.res$Direction[monocyte.rlm.res$COEFF > 0 & monocyte.rlm.res$Sig == 1] <- "More"


# give the features more informative/better formated names
monocyte.rlm.res$Predictor[monocyte.rlm.res$Predictor == "N_CpG"] <- "CpG island"
monocyte.rlm.res$Predictor[monocyte.rlm.res$Predictor == "GC"] <- "Promoter GC %"
monocyte.rlm.res$Predictor[monocyte.rlm.res$Predictor == "SP1"] <- "Number SP1 motifs"
monocyte.rlm.res$Predictor[monocyte.rlm.res$Predictor == "TBP"] <- "Number TBP motifs"
monocyte.rlm.res$Predictor[monocyte.rlm.res$Predictor == "MED_PHAST"] <- "Median PhastCons"
monocyte.rlm.res$Predictor[monocyte.rlm.res$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
monocyte.rlm.res$Predictor[monocyte.rlm.res$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
monocyte.rlm.res$Predictor[monocyte.rlm.res$Predictor == "NALIGN_PHAST"] <- "Numer aligned bases"
monocyte.rlm.res$Predictor[monocyte.rlm.res$Predictor == "EXON_AVLENGTH"] <- "Mean exon length"
monocyte.rlm.res$Predictor[monocyte.rlm.res$Predictor == "EXON_COUNT"] <- "Number of exons"
monocyte.rlm.res$Predictor[monocyte.rlm.res$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"
monocyte.rlm.res$Predictor[monocyte.rlm.res$Predictor == "EXON_VARLENGTH"] <- "Exon length variance"
monocyte.rlm.res$Predictor[monocyte.rlm.res$Predictor == "Mean"] <- "Mean expression"

monocyte.rlm.res <- monocyte.rlm.res[!grepl(monocyte.rlm.res$Predictor, pattern="Intercept"), ]


multivar.plot <- ggplot(monocyte.rlm.res,
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
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-monocyte_multiivariateLM.png",
       height=4.75, width=6.75, dpi=300)

# plot the uni- and multivariate on the same graph?
monocyte.rlm.df$Model <- "Univariate"
monocyte.rlm.res$Model <- "Multivariate"

monocyte.rlm.all <- do.call(rbind.data.frame, list("Univariate"=monocyte.rlm.df,
                                                "Multivariate"=monocyte.rlm.res))

all.plot <- ggplot(monocyte.rlm.all,
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
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-monocyte_allLM.png",
       height=5.6, width=7.75, dpi=300)

###############################
## Plot CpG islands vs alpha ##
###############################
monocyte.match$CpGisland <- factor(monocyte.match$N_CpG,
                                levels=c(0, 1),
                                labels=c("Non-CpG island", "CpG island"))
cpg.cols <- c("#027E00", "#00DDEC")
names(cpg.cols) <- levels(monocyte.match$CpGisland)

tc_cpg <- ggplot(monocyte.match, aes(x=CpGisland, y=Alpha_r, colour=CpGisland)) + 
  theme_mike() + 
  geom_jitter(position=position_jitterdodge(jitter.height=0,
                                            jitter.width=1.5),
              alpha=0.7) +
  geom_boxplot(width=0.5, fill='white', colour='black') +
  scale_colour_manual(values=cpg.cols) +
  labs(x="Overlapping CpG island", y=expression(paste(alpha["r"], " Overdispersion"))) +
  guides(colour=FALSE)

ggsave(tc_cpg,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-monocyte_boxplot_CpGislands-overdispersion.png",
       height=3.75, width=4.75, dpi=300)

tc_cpgMe <- ggplot(monocyte.match, aes(x=CpGisland, y=Mean, colour=CpGisland)) + 
  theme_mike() + 
  geom_jitter(position=position_jitterdodge(jitter.height=0,
                                            jitter.width=1.5),
              alpha=0.7) +
  geom_boxplot(width=0.5, fill='white', colour='black') +
  scale_colour_manual(values=cpg.cols) +
  labs(x="Overlapping CpG island", y=expression(paste("Mean log"[2], " Expression"))) +
  guides(colour=FALSE)

ggsave(tc_cpgMe,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-monocyte_boxplot_CpGislands-mean.png",
       height=3.75, width=4.75, dpi=300)
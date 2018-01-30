source("~/Dropbox/R_sessions/Noise/human_genomic_noise_features.R")
source("~/Dropbox/R_sessions/Noise/human_pancreas_genomic_noise_features-alt_data.R")

library(ggplot2)
library(robustbase)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

panc.vars <- colnames(genomic.features)
panc.vars <- panc.vars[!grepl(panc.vars, pattern="(NMI)|(CGI_SIZE)|(cpg_)|(Mean)")]

panc.genomic.vars <- paste(panc.vars,
                           collapse=" + ")

panc.match <- merge(panc.gene.summary, genomic.features,
                    by='GENE')

#########################################
## univariate robust linear regression ##
#########################################
# this is plotting for presentation/manuscript figures
panc.univariate_list <- list()
model.control <- lmrob.control(max.it=500, k.max=500, rel.tol=1e-7)

panc.var.names <- unlist(strsplit(panc.genomic.vars, split=" + ", fixed=T))
# remove redundant variables, i.e CpG island AND NMI
# remove CpG island characteristics
panc.var.names <- panc.var.names[!grepl(panc.var.names, 
                                        pattern="(NMI)|(CGI_SIZE)|(cpg_)|(CGI)|(GENE)|(Mean)")]

for(x in seq_along(panc.var.names)){
  .variable <- paste(panc.var.names[x], sep=" + ")
  .glm.form <- as.formula(paste("Residual.CV2", .variable, sep=" ~ "))
  
  m.rlm <- lmrob(.glm.form, data=panc.match, control=model.control)
  m.robust <- summary(m.rlm)
  m.rlm.res <- as.data.frame(m.robust$coefficients)
  m.res.mat <- as.matrix(m.rlm.res)
  m.res.var <- m.res.mat[2, ]
  names(m.res.var) <- c("COEFF", "SE", "STAT", "P")
  panc.univariate_list[[panc.var.names[x]]] <- as.list(m.res.var)
}

panc.rlm.df <- do.call(rbind.data.frame, panc.univariate_list)
panc.rlm.df$Predictor <- rownames(panc.rlm.df)
panc.rlm.df$Padjust <- p.adjust(panc.rlm.df$P)
panc.rlm.df$Sig <- as.numeric(panc.rlm.df$Padjust <= 0.01)
panc.rlm.df$Direction <- "NoEffect"
panc.rlm.df$Direction[panc.rlm.df$COEFF < 0 & panc.rlm.df$Sig == 1] <- "Less"
panc.rlm.df$Direction[panc.rlm.df$COEFF > 0 & panc.rlm.df$Sig == 1] <- "More"

# give the features more informative/better formated names
panc.rlm.df$Predictor[panc.rlm.df$Predictor == "Mean"] <- "Mean expression"
panc.rlm.df$Predictor[panc.rlm.df$Predictor == "N_CpG"] <- "CpG island"
panc.rlm.df$Predictor[panc.rlm.df$Predictor == "GC"] <- "Promoter GC %"
panc.rlm.df$Predictor[panc.rlm.df$Predictor == "SP1"] <- "Number SP1 motifs"
panc.rlm.df$Predictor[panc.rlm.df$Predictor == "TBP"] <- "Number TBP motifs"
panc.rlm.df$Predictor[panc.rlm.df$Predictor == "MED_PHAST"] <- "Median PhastCons"
panc.rlm.df$Predictor[panc.rlm.df$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
panc.rlm.df$Predictor[panc.rlm.df$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
panc.rlm.df$Predictor[panc.rlm.df$Predictor == "NALIGN_PHAST"] <- "Numer aligned bases"
panc.rlm.df$Predictor[panc.rlm.df$Predictor == "EXON_AVLENGTH"] <- "Mean exon length"
panc.rlm.df$Predictor[panc.rlm.df$Predictor == "EXON_COUNT"] <- "Number of exons"
panc.rlm.df$Predictor[panc.rlm.df$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"
panc.rlm.df$Predictor[panc.rlm.df$Predictor == "EXON_VARLENGTH"] <- "Exon length variance"

panc.rlm.df$Tissue <- "AlphaIslet"
panc.rlm.df$Species <- "Human"

panc.rlm.df$Direction <- factor(panc.rlm.df$Direction,
                                 levels=c("Less", "NoEffect", "More"),
                                 labels=c("Less", "NoEffect", "More"))

effect.cols <- c("#62148f" , "#878787", "#feaf10")
names(effect.cols) <- levels(panc.rlm.df$Direction)

univar.plot <- ggplot(panc.rlm.df,
                      aes(x=reorder(Predictor, -STAT),
                          y=STAT, fill=Direction)) +
  geom_point(alpha=0.55, shape=21, size=5) + theme_mike() +
  scale_fill_manual(values=effect.cols) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Annotation", y="t-statistic") +
  guides(fill=FALSE) +
  geom_hline(mapping=aes(yintercept=0), linetype="dashed", colour="grey") 

ggsave(univar.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-AlphaIslet-alt_Data_univariateLM.png",
       height=4.75, width=6.75, dpi=300)


###########################################
## multivariate robust linear regression ##
###########################################
panc.genomic.vars <- paste(panc.var.names,
                           collapse=" + ")

panc.glm.form <- as.formula(paste("Residual.CV2",
                                  panc.genomic.vars, sep=" ~ "))

panc.rlm <- lmrob(panc.glm.form, data=panc.match, control=model.control)
panc.robust <- summary(panc.rlm)
panc.rlm.res <- as.data.frame(panc.robust$coefficients)
panc.rlm.res$Padjust <- p.adjust(panc.rlm.res$`Pr(>|t|)`)
panc.rlm.res$Sig <- as.numeric(panc.rlm.res$Padjust <= 0.01)
panc.rlm.res$Predictor <- rownames(panc.rlm.res)
colnames(panc.rlm.res) <- c("COEFF", "SE", "STAT", "P", "Padjust", "Sig", "Predictor")
panc.rlm.res$Direction <- "NoEffect"
panc.rlm.res$Direction[panc.rlm.res$COEFF < 0 & panc.rlm.res$Sig == 1] <- "Less"
panc.rlm.res$Direction[panc.rlm.res$COEFF > 0 & panc.rlm.res$Sig == 1] <- "More"

# give the features more informative/better formated names
panc.rlm.res$Predictor[panc.rlm.res$Predictor == "N_CpG"] <- "CpG island"
panc.rlm.res$Predictor[panc.rlm.res$Predictor == "GC"] <- "Promoter GC %"
panc.rlm.res$Predictor[panc.rlm.res$Predictor == "SP1"] <- "Number SP1 motifs"
panc.rlm.res$Predictor[panc.rlm.res$Predictor == "TBP"] <- "Number TBP motifs"
panc.rlm.res$Predictor[panc.rlm.res$Predictor == "MED_PHAST"] <- "Median PhastCons"
panc.rlm.res$Predictor[panc.rlm.res$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
panc.rlm.res$Predictor[panc.rlm.res$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
panc.rlm.res$Predictor[panc.rlm.res$Predictor == "NALIGN_PHAST"] <- "Numer aligned bases"
panc.rlm.res$Predictor[panc.rlm.res$Predictor == "EXON_AVLENGTH"] <- "Mean exon length"
panc.rlm.res$Predictor[panc.rlm.res$Predictor == "EXON_COUNT"] <- "Number of exons"
panc.rlm.res$Predictor[panc.rlm.res$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"
panc.rlm.res$Predictor[panc.rlm.res$Predictor == "EXON_VARLENGTH"] <- "Exon length variance"
panc.rlm.res$Predictor[panc.rlm.res$Predictor == "Mean"] <- "Mean expression"

panc.rlm.res <- panc.rlm.res[!grepl(panc.rlm.res$Predictor, pattern="Intercept"), ]

panc.rlm.res$Tissue <- "AlphaIslet"
panc.rlm.res$Species <- "Human"
panc.rlm.res$Direction <- factor(panc.rlm.res$Direction,
                                 levels=c("Less", "NoEffect", "More"),
                                 labels=c("Less", "NoEffect", "More"))
write.table(panc.rlm.res,
            file="~/Dropbox/Noise_genomics/Model_results/human_AlphaIslet-alt_data_multivariateRLM.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")

multivar.plot <- ggplot(panc.rlm.res,
                        aes(x=reorder(Predictor, -STAT),
                            y=STAT, fill=Direction)) +
  geom_point(alpha=0.55, shape=21, size=5) + theme_mike() +
  scale_fill_manual(values=effect.cols) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Annotation", y="t-statistic") +
  guides(fill=FALSE) +
  geom_hline(mapping=aes(yintercept=0), linetype="dashed", colour="grey")

ggsave(multivar.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-AlphaIslet-alt_data_multiivariateLM.png",
       height=4.75, width=6.75, dpi=300)

# plot the uni- and multivariate on the same graph?
panc.rlm.df$Model <- "Univariate"
panc.rlm.res$Model <- "Multivariate"

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
  scale_fill_manual(values=effect.cols) +
  scale_shape_manual(values=c(21, 23)) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Annotation", y="t-statistic") +
  guides(fill=FALSE, shape=FALSE)

ggsave(all.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-AlphaIslet-alt_Data_allLM.png",
       height=5.6, width=7.75, dpi=300)

###############################
## Plot CpG islands vs alpha ##
###############################
panc.match$CpGisland <- factor(panc.match$N_CpG,
                               levels=c(0, 1),
                               labels=c("Non-CpG island", "CpG island"))
cpg.cols <- c("#027E00", "#00DDEC")
names(cpg.cols) <- levels(panc.match$CpGisland)

tc_cpg <- ggplot(panc.match, aes(x=CpGisland, y=Residual.CV2, colour=CpGisland)) + 
  theme_mike() + 
  geom_jitter(position=position_jitterdodge(jitter.height=0,
                                            jitter.width=1.5),
              alpha=0.7) +
  geom_boxplot(width=0.5, fill='white', colour='black') +
  scale_colour_manual(values=cpg.cols) +
  labs(x="Overlapping CpG island", y=expression(paste(alpha["r"], " Overdispersion"))) +
  guides(colour=FALSE)

ggsave(tc_cpg,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-AlphaIslet-alt_data-boxplot_CpGislands-overdispersion.png",
       height=3.75, width=4.75, dpi=300)

tc_cpgMe <- ggplot(panc.match, aes(x=CpGisland, y=Mean, colour=CpGisland)) + 
  theme_mike() + 
  geom_jitter(position=position_jitterdodge(jitter.height=0,
                                            jitter.width=1.5),
              alpha=0.7) +
  geom_boxplot(width=0.5, fill='white', colour='black') +
  scale_colour_manual(values=cpg.cols) +
  labs(x="Overlapping CpG island", y=expression(paste("Mean log"[2], " Expression"))) +
  guides(colour=FALSE)

ggsave(tc_cpgMe,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-AlphaIslet-alt_data-boxplot_CpGislands-mean.png",
       height=3.75, width=4.75, dpi=300)
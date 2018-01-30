source("~/Dropbox/R_sessions/Noise/mESC_chromHMM.R")
source("~/Dropbox/R_sessions/Noise/genomic_noise_features.R")
source("~/Dropbox/R_sessions/Noise/mESC_genomic_noise_features.R")

library(ggplot2)
library(robustbase)
library(gtools)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

mesc.vars <- colnames(genomic.features)[2:16]
mesc.vars <- mesc.vars[!grepl(mesc.vars, pattern="(NMI)|(CGI_SIZE)|(cpg_)|(PHAST)")]

mesc.genomic.vars <- paste(mesc.vars,
                           collapse=" + ")

mesc.match <- merge(mesc.gene.summary, genomic.features,
                    by='GENE')

# what about using the ranks?
mesc.match$Noise.Rank <- order(mesc.match$CV2, decreasing=TRUE)
#########################################
## univariate robust linear regression ##
#########################################
# this is plotting for presentation/manuscript figures
mesc.univariate_list <- list()
model.control <- lmrob.control(max.it=500, k.max=500, rel.tol=1e-7)

mesc.var.names <- unlist(strsplit(mesc.genomic.vars, split=" + ", fixed=T))
# remove redundant variables, i.e CpG island AND NMI
# remove CpG island characteristics
mesc.var.names <- mesc.var.names[!grepl(mesc.var.names, pattern="(NMI)|(CGI_SIZE)|(cpg_)|(PHAST)")]

for(x in seq_along(mesc.var.names)){
  .variable <- paste(mesc.var.names[x], collapse=" + ")
  .glm.form <- as.formula(paste("Residual.CV2", .variable, sep=" ~ "))

  m.rlm <- glm(.glm.form, data=mesc.match)
  m.robust <- summary(m.rlm)
  m.rlm.res <- as.data.frame(m.robust$coefficients)
  m.rlm.res$Sig <- as.numeric(m.rlm.res$`Pr(>|t|)` <= 0.05)
  m.res.mat <- as.matrix(m.rlm.res)
  m.res.var <- m.res.mat[2, ]
  names(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
  mesc.univariate_list[[mesc.var.names[x]]] <- as.list(m.res.var)
}

mesc.rlm.df <- do.call(rbind.data.frame, mesc.univariate_list)
mesc.rlm.df$Predictor <- rownames(mesc.rlm.df)
mesc.rlm.df$Padjust <- p.adjust(mesc.rlm.df$P)
mesc.rlm.df$Sig <- as.numeric(mesc.rlm.df$Padjust <= 0.01)
mesc.rlm.df$Direction <- "NoEffect"
mesc.rlm.df$Direction[mesc.rlm.df$COEFF < 0 & mesc.rlm.df$Sig == 1] <- "Less"
mesc.rlm.df$Direction[mesc.rlm.df$COEFF > 0 & mesc.rlm.df$Sig == 1] <- "More"

# give the features more informative/better formated names
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "Mean"] <- "Mean expression"
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "N_CpG"] <- "CpG island"
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "GC"] <- "Promoter GC %"
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "SP1"] <- "Number SP1 motifs"
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "TBP"] <- "Number TBP motifs"
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "MED_PHAST"] <- "Median PhastCons"
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "NALIGN_PHAST"] <- "Numer aligned bases"
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "EXON_AVLENGTH"] <- "Mean exon length"
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "EXON_COUNT"] <- "Number of exons"
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "EXON_VARLENGTH"] <- "Exon length variance"

mesc.rlm.df$Tissue <- "ESC"
mesc.rlm.df$Species <- "Mouse"

univar.plot <- ggplot(mesc.rlm.df,
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
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mESC_univariateLM.png",
       height=4.75, width=6.75, dpi=300)


###########################################
## multivariate robust linear regression ##
###########################################

mesc.glm.form <- as.formula(paste("Residual.CV2",
                                  mesc.genomic.vars, sep=" ~ "))

mesc.rlm <- lmrob(mesc.glm.form, data=mesc.match, control=model.control)
mesc.robust <- summary(mesc.rlm)
mesc.rlm.res <- as.data.frame(mesc.robust$coefficients)
mesc.rlm.res$Padjust <- p.adjust(mesc.rlm.res$`Pr(>|t|)`)
mesc.rlm.res$Sig <- as.numeric(mesc.rlm.res$Padjust <= 0.01)
mesc.rlm.res$Predictor <- rownames(mesc.rlm.res)
colnames(mesc.rlm.res) <- c("COEFF", "SE", "STAT", "P", "Padjust", "Sig", "Predictor")
mesc.rlm.res$Direction <- "NoEffect"
mesc.rlm.res$Direction[mesc.rlm.res$COEFF < 0 & mesc.rlm.res$Sig == 1] <- "Less"
mesc.rlm.res$Direction[mesc.rlm.res$COEFF > 0 & mesc.rlm.res$Sig == 1] <- "More"


# give the features more informative/better formated names
mesc.rlm.res$Predictor[mesc.rlm.res$Predictor == "N_CpG"] <- "CpG island"
mesc.rlm.res$Predictor[mesc.rlm.res$Predictor == "GC"] <- "Promoter GC %"
mesc.rlm.res$Predictor[mesc.rlm.res$Predictor == "SP1"] <- "Number SP1 motifs"
mesc.rlm.res$Predictor[mesc.rlm.res$Predictor == "TBP"] <- "Number TBP motifs"
mesc.rlm.res$Predictor[mesc.rlm.res$Predictor == "MED_PHAST"] <- "Median PhastCons"
mesc.rlm.res$Predictor[mesc.rlm.res$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
mesc.rlm.res$Predictor[mesc.rlm.res$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
mesc.rlm.res$Predictor[mesc.rlm.res$Predictor == "NALIGN_PHAST"] <- "Numer aligned bases"
mesc.rlm.res$Predictor[mesc.rlm.res$Predictor == "EXON_AVLENGTH"] <- "Mean exon length"
mesc.rlm.res$Predictor[mesc.rlm.res$Predictor == "EXON_COUNT"] <- "Number of exons"
mesc.rlm.res$Predictor[mesc.rlm.res$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"
mesc.rlm.res$Predictor[mesc.rlm.res$Predictor == "EXON_VARLENGTH"] <- "Exon length variance"
mesc.rlm.res$Predictor[mesc.rlm.res$Predictor == "Mean"] <- "Mean expression"

mesc.rlm.res <- mesc.rlm.res[!grepl(mesc.rlm.res$Predictor, pattern="Intercept"), ]

mesc.rlm.res$Tissue <- "ESC"
mesc.rlm.res$Species <- "Mouse"

write.table(mesc.rlm.res,
            file="~/Dropbox/Noise_genomics/Model_results/mouse_ESC_multivariateRLM.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")


multivar.plot <- ggplot(mesc.rlm.res,
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
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mESC_multiivariateLM.png",
       height=4.75, width=6.75, dpi=300)

# plot the uni- and multivariate on the same graph?
mesc.rlm.df$Model <- "Univariate"
mesc.rlm.res$Model <- "Multivariate"

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
  scale_y_continuous(limits=c(-25, 25), oob=squish)

ggsave(all.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mESC_allLM.png",
       height=5.6, width=7.75, dpi=300)

##############################
## Plot CpG islands vs rCV2 ##
##############################
mesc.match$CpGisland <- factor(mesc.match$N_CpG,
                                levels=c(0, 1),
                                labels=c("Non-CpG island", "CpG island"))
cpg.cols <- c("#027E00", "#00DDEC")
names(cpg.cols) <- levels(mesc.match$CpGisland)

tc_cpg <- ggplot(mesc.match, aes(x=CpGisland, y=Residual.CV2, colour=CpGisland)) + 
  theme_mike() + 
  geom_jitter(position=position_jitterdodge(jitter.height=0,
                                            jitter.width=1.5),
              alpha=0.7) +
  geom_boxplot(width=0.5, fill='white', colour='black') +
  scale_colour_manual(values=cpg.cols) +
  labs(x="Overlapping CpG island", y=expression(paste("Residual CV"^2))) +
  guides(colour=FALSE) +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16))

ggsave(tc_cpg,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mESC_boxplot_CpGislands-overdispersion.png",
       height=3.75, width=4.75, dpi=300)

tc_cpgMe <- ggplot(mesc.match, aes(x=CpGisland, y=Mean, colour=CpGisland)) + 
  theme_mike() + 
  geom_jitter(position=position_jitterdodge(jitter.height=0,
                                            jitter.width=1.5),
              alpha=0.7) +
  geom_boxplot(width=0.5, fill='white', colour='black') +
  scale_colour_manual(values=cpg.cols) +
  labs(x="Overlapping CpG island", y=expression(paste("Mean log"[2], " Expression"))) +
  guides(colour=FALSE)

ggsave(tc_cpgMe,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mESC_boxplot_CpGislands-mean.png",
       height=3.75, width=4.75, dpi=300)

source("~/Dropbox/R_sessions/Noise/mESC_chromHMM.R")
source("~/Dropbox/R_sessions/Noise/human_genomic_noise_features.R")
source("~/Dropbox/R_sessions/Noise/hESC_genomic_noise_features.R")

library(ggplot2)
library(MASS)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

hesc.vars <- colnames(genomic.features)
hesc.vars <- hesc.vars[!grepl(hesc.vars, pattern="(NMI)|(CGI_SIZE)|(cpg_)")]

hesc.genomic.vars <- paste(c("Mean",
                             hesc.vars),
                           collapse=" + ")

hesc.match <- merge(hesc.gene.summary, genomic.features,
                    by='GENE')

#########################################
## univariate robust linear regression ##
#########################################
# this is plotting for presentation/manuscript figures
hesc.univariate_list <- list()

hesc.var.names <- unlist(strsplit(hesc.genomic.vars, split=" + ", fixed=T))
# remove redundant variables, i.e CpG island AND NMI
# remove CpG island characteristics
hesc.var.names <- hesc.var.names[!grepl(hesc.var.names, pattern="(NMI)|(CGI_SIZE)|(cpg_)|(CGI)|(GENE)")]

for(x in seq_along(hesc.var.names)){
  .variable <- paste(hesc.var.names[x], sep=" + ")
  .glm.form <- as.formula(paste("Residual.CV2", .variable, sep=" ~ "))
  
  m.rlm <- rlm(.glm.form, data=hesc.match)
  m.robust <- summary(m.rlm)
  m.rlm.res <- as.data.frame(m.robust$coefficients)
  m.rlm.res$Pval <- 2*pt(-abs(m.rlm.res[, 3]), df=3)
  m.rlm.res$Sig <- as.numeric(m.rlm.res$Pval <= 0.05)
  m.res.mat <- as.matrix(m.rlm.res)
  m.res.var <- m.res.mat[2, ]
  names(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
  hesc.univariate_list[[hesc.var.names[x]]] <- as.list(m.res.var)
}

hesc.rlm.df <- do.call(rbind.data.frame, hesc.univariate_list)
hesc.rlm.df$Predictor <- rownames(hesc.rlm.df)
hesc.rlm.df$Direction <- "NoEffect"
hesc.rlm.df$Direction[hesc.rlm.df$COEFF < 0 & hesc.rlm.df$Sig == 1] <- "Less"
hesc.rlm.df$Direction[hesc.rlm.df$COEFF > 0 & hesc.rlm.df$Sig == 1] <- "More"

# give the features more informative/better formated names
hesc.rlm.df$Predictor[hesc.rlm.df$Predictor == "Mean"] <- "Mean expression"
hesc.rlm.df$Predictor[hesc.rlm.df$Predictor == "N_CpG"] <- "CpG island"
hesc.rlm.df$Predictor[hesc.rlm.df$Predictor == "GC"] <- "Promoter GC %"
hesc.rlm.df$Predictor[hesc.rlm.df$Predictor == "SP1"] <- "Number SP1 motifs"
hesc.rlm.df$Predictor[hesc.rlm.df$Predictor == "TBP"] <- "Number TBP motifs"
hesc.rlm.df$Predictor[hesc.rlm.df$Predictor == "MED_PHAST"] <- "Median PhastCons"
hesc.rlm.df$Predictor[hesc.rlm.df$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
hesc.rlm.df$Predictor[hesc.rlm.df$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
hesc.rlm.df$Predictor[hesc.rlm.df$Predictor == "NALIGN_PHAST"] <- "Numer aligned bases"
hesc.rlm.df$Predictor[hesc.rlm.df$Predictor == "EXON_AVLENGTH"] <- "Mean exon length"
hesc.rlm.df$Predictor[hesc.rlm.df$Predictor == "EXON_COUNT"] <- "Number of exons"
hesc.rlm.df$Predictor[hesc.rlm.df$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"
hesc.rlm.df$Predictor[hesc.rlm.df$Predictor == "EXON_VARLENGTH"] <- "Exon length variance"

hesc.rlm.df$Tissue <- "ESC"
hesc.rlm.df$Species <- "Human"

univar.plot <- ggplot(hesc.rlm.df,
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
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/hESC_univariateLM.png",
       height=4.75, width=6.75, dpi=300)


###########################################
## multivariate robust linear regression ##
###########################################
hesc.genomic.vars <- paste(c("Mean",
                             hesc.var.names),
                           collapse=" + ")

hesc.glm.form <- as.formula(paste("Residual.CV2",
                                  hesc.genomic.vars, sep=" ~ "))

hesc.rlm <- rlm(hesc.glm.form, data=hesc.match)
hesc.robust <- summary(hesc.rlm)
hesc.rlm.res <- as.data.frame(hesc.robust$coefficients)
hesc.rlm.res$Pval <- 2*pt(-abs(hesc.rlm.res[, 3]), df=dim(hesc.match)[2]-1)
hesc.rlm.res$Sig <- as.numeric(hesc.rlm.res$Pval <= 0.05)
hesc.rlm.res$Predictor <- rownames(hesc.rlm.res)
colnames(hesc.rlm.res) <- c("COEFF", "SE", "STAT", "P", "Sig", "Predictor")
hesc.rlm.res$Direction <- "NoEffect"
hesc.rlm.res$Direction[hesc.rlm.res$COEFF < 0 & hesc.rlm.res$Sig == 1] <- "Less"
hesc.rlm.res$Direction[hesc.rlm.res$COEFF > 0 & hesc.rlm.res$Sig == 1] <- "More"


# give the features more informative/better formated names
hesc.rlm.res$Predictor[hesc.rlm.res$Predictor == "N_CpG"] <- "CpG island"
hesc.rlm.res$Predictor[hesc.rlm.res$Predictor == "GC"] <- "Promoter GC %"
hesc.rlm.res$Predictor[hesc.rlm.res$Predictor == "SP1"] <- "Number SP1 motifs"
hesc.rlm.res$Predictor[hesc.rlm.res$Predictor == "TBP"] <- "Number TBP motifs"
hesc.rlm.res$Predictor[hesc.rlm.res$Predictor == "MED_PHAST"] <- "Median PhastCons"
hesc.rlm.res$Predictor[hesc.rlm.res$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
hesc.rlm.res$Predictor[hesc.rlm.res$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
hesc.rlm.res$Predictor[hesc.rlm.res$Predictor == "NALIGN_PHAST"] <- "Numer aligned bases"
hesc.rlm.res$Predictor[hesc.rlm.res$Predictor == "EXON_AVLENGTH"] <- "Mean exon length"
hesc.rlm.res$Predictor[hesc.rlm.res$Predictor == "EXON_COUNT"] <- "Number of exons"
hesc.rlm.res$Predictor[hesc.rlm.res$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"
hesc.rlm.res$Predictor[hesc.rlm.res$Predictor == "EXON_VARLENGTH"] <- "Exon length variance"
hesc.rlm.res$Predictor[hesc.rlm.res$Predictor == "Mean"] <- "Mean expression"

hesc.rlm.res <- hesc.rlm.res[!grepl(hesc.rlm.res$Predictor, pattern="Intercept"), ]
hesc.rlm.res$Tissue <- "ESC"
hesc.rlm.res$Species <- "Human"

write.table(hesc.rlm.res,
            file="~/Dropbox/Noise_genomics/Model_results/human_ESC_multivariateRLM.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")

multivar.plot <- ggplot(hesc.rlm.res,
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
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/hESC_multiivariateLM.png",
       height=4.75, width=6.75, dpi=300)

# plot the uni- and multivariate on the same graph?
hesc.rlm.df$Model <- "Univariate"
hesc.rlm.res$Model <- "Multivariate"

hesc.rlm.all <- do.call(rbind.data.frame, list("Univariate"=hesc.rlm.df,
                                               "Multivariate"=hesc.rlm.res))

all.plot <- ggplot(hesc.rlm.all,
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
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/hESC_allLM.png",
       height=5.6, width=7.75, dpi=300)

###############################
## Plot CpG islands vs alpha ##
###############################
hesc.match$CpGisland <- factor(hesc.match$N_CpG,
                               levels=c(0, 1),
                               labels=c("Non-CpG island", "CpG island"))
cpg.cols <- c("#027E00", "#00DDEC")
names(cpg.cols) <- levels(hesc.match$CpGisland)

tc_cpg <- ggplot(hesc.match, aes(x=CpGisland, y=Residual.CV2, colour=CpGisland)) + 
  theme_mike() + 
  geom_jitter(position=position_jitterdodge(jitter.height=0,
                                            jitter.width=1.5),
              alpha=0.7) +
  geom_boxplot(width=0.5, fill='white', colour='black') +
  scale_colour_manual(values=cpg.cols) +
  labs(x="Overlapping CpG island", y=expression(paste(alpha["r"], " Overdispersion"))) +
  guides(colour=FALSE)

ggsave(tc_cpg,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/hESC_boxplot_CpGislands-overdispersion.png",
       height=3.75, width=4.75, dpi=300)

tc_cpgMe <- ggplot(hesc.match, aes(x=CpGisland, y=Mean, colour=CpGisland)) + 
  theme_mike() + 
  geom_jitter(position=position_jitterdodge(jitter.height=0,
                                            jitter.width=1.5),
              alpha=0.7) +
  geom_boxplot(width=0.5, fill='white', colour='black') +
  scale_colour_manual(values=cpg.cols) +
  labs(x="Overlapping CpG island", y=expression(paste("Mean log"[2], " Expression"))) +
  guides(colour=FALSE)

ggsave(tc_cpgMe,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/hESC_boxplot_CpGislands-mean.png",
       height=3.75, width=4.75, dpi=300)
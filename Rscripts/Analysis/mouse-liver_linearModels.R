source("~/Dropbox/R_sessions/Noise/mESC_chromHMM.R")
source("~/Dropbox/R_sessions/Noise/genomic_noise_features.R")
source("~/Dropbox/R_sessions/Noise/mouse_liver_noise_features.R")
library(MASS)
library(ggplot2)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

liver.vars <- colnames(genomic.features)[2:16]
liver.vars <- liver.vars[!grepl(liver.vars, pattern="(NMI)|(CGI_SIZE)|(cpg_)")]

liver.genomic.vars <- paste(c("Mean",
                             liver.vars),
                           collapse=" + ")

liver.match <- merge(liver.gene.summary, genomic.features,
                    by='GENE')

#########################################
## univariate robust linear regression ##
#########################################
# this is plotting for presentation/manuscript figures
liver.univariate_list <- list()

liver.var.names <- unlist(strsplit(liver.genomic.vars, split=" + ", fixed=T))
# remove redundant variables, i.e CpG island AND NMI
# remove CpG island characteristics
liver.var.names <- liver.var.names[!grepl(liver.var.names, pattern="(NMI)|(CGI_SIZE)|(cpg_)")]

for(x in seq_along(liver.var.names)){
  .variable <- paste(liver.var.names[x], sep=" + ")
  .glm.form <- as.formula(paste("Residual.CV2", .variable, sep=" ~ "))
  
  m.rlm <- rlm(.glm.form, data=liver.match)
  m.robust <- summary(m.rlm)
  m.rlm.res <- as.data.frame(m.robust$coefficients)
  m.rlm.res$Pval <- 2*pt(-abs(m.rlm.res[, 3]), df=3)
  m.rlm.res$Sig <- as.numeric(m.rlm.res$Pval <= 0.05)
  m.res.mat <- as.matrix(m.rlm.res)
  m.res.var <- m.res.mat[2, ]
  names(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
  liver.univariate_list[[liver.var.names[x]]] <- as.list(m.res.var)
}

liver.rlm.df <- do.call(rbind.data.frame, liver.univariate_list)
liver.rlm.df$Predictor <- rownames(liver.rlm.df)
liver.rlm.df$Direction <- "NoEffect"
liver.rlm.df$Direction[liver.rlm.df$COEFF < 0 & liver.rlm.df$Sig == 1] <- "Less"
liver.rlm.df$Direction[liver.rlm.df$COEFF > 0 & liver.rlm.df$Sig == 1] <- "More"

# give the features more informative/better formated names
liver.rlm.df$Predictor[liver.rlm.df$Predictor == "Mean"] <- "Mean expression"
liver.rlm.df$Predictor[liver.rlm.df$Predictor == "N_CpG"] <- "CpG island"
liver.rlm.df$Predictor[liver.rlm.df$Predictor == "GC"] <- "Promoter GC %"
liver.rlm.df$Predictor[liver.rlm.df$Predictor == "SP1"] <- "Number SP1 motifs"
liver.rlm.df$Predictor[liver.rlm.df$Predictor == "TBP"] <- "Number TBP motifs"
liver.rlm.df$Predictor[liver.rlm.df$Predictor == "MED_PHAST"] <- "Median PhastCons"
liver.rlm.df$Predictor[liver.rlm.df$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
liver.rlm.df$Predictor[liver.rlm.df$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
liver.rlm.df$Predictor[liver.rlm.df$Predictor == "NALIGN_PHAST"] <- "Numer aligned bases"
liver.rlm.df$Predictor[liver.rlm.df$Predictor == "EXON_AVLENGTH"] <- "Mean exon length"
liver.rlm.df$Predictor[liver.rlm.df$Predictor == "EXON_COUNT"] <- "Number of exons"
liver.rlm.df$Predictor[liver.rlm.df$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"
liver.rlm.df$Predictor[liver.rlm.df$Predictor == "EXON_VARLENGTH"] <- "Exon length variance"

liver.rlm.df$Tissue <- "Liver"
liver.rlm.df$Species <- "Mouse"

univar.plot <- ggplot(liver.rlm.df,
                      aes(x=reorder(Predictor, -STAT),
                          y=STAT, fill=Direction)) +
  geom_point(alpha=0.55, shape=21, size=5) + theme_mike() +
  scale_fill_manual(values=c("#62148f", "#feaf10", "#878787")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Annotation", y="t-statistic") +
  guides(fill=FALSE) +
  geom_hline(mapping=aes(yintercept=0), linetype="dashed", colour="grey") +
  scale_y_continuous(limits=c(-60, 60))

ggsave(univar.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mouse-liver_univariateLM.png",
       height=4.75, width=6.75, dpi=300)

###########################################
## multivariate robust linear regression ##
###########################################

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
liver.rlm.res$Predictor[liver.rlm.res$Predictor == "N_CpG"] <- "CpG island"
liver.rlm.res$Predictor[liver.rlm.res$Predictor == "GC"] <- "Promoter GC %"
liver.rlm.res$Predictor[liver.rlm.res$Predictor == "SP1"] <- "Number SP1 motifs"
liver.rlm.res$Predictor[liver.rlm.res$Predictor == "TBP"] <- "Number TBP motifs"
liver.rlm.res$Predictor[liver.rlm.res$Predictor == "MED_PHAST"] <- "Median PhastCons"
liver.rlm.res$Predictor[liver.rlm.res$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
liver.rlm.res$Predictor[liver.rlm.res$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
liver.rlm.res$Predictor[liver.rlm.res$Predictor == "NALIGN_PHAST"] <- "Numer aligned bases"
liver.rlm.res$Predictor[liver.rlm.res$Predictor == "EXON_AVLENGTH"] <- "Mean exon length"
liver.rlm.res$Predictor[liver.rlm.res$Predictor == "EXON_COUNT"] <- "Number of exons"
liver.rlm.res$Predictor[liver.rlm.res$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"
liver.rlm.res$Predictor[liver.rlm.res$Predictor == "EXON_VARLENGTH"] <- "Exon length variance"
liver.rlm.res$Predictor[liver.rlm.res$Predictor == "Mean"] <- "Mean expression"

liver.rlm.res <- liver.rlm.res[!grepl(liver.rlm.res$Predictor, pattern="Intercept"), ]
liver.rlm.res$Tissue <- "Liver"
liver.rlm.res$Species <- "Mouse"

write.table(liver.rlm.res,
            file="~/Dropbox/Noise_genomics/Model_results/mouse_liver_multivariateRLM.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")

multivar.plot <- ggplot(liver.rlm.res,
                        aes(x=reorder(Predictor, -STAT),
                            y=STAT, fill=Direction)) +
  geom_point(alpha=0.55, shape=21, size=5) + theme_mike() +
  scale_fill_manual(values=c("#62148f", "#feaf10", "#878787")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Annotation", y="t-statistic") +
  guides(fill=FALSE) +
  geom_hline(mapping=aes(yintercept=0), linetype="dashed", colour="grey") +
  scale_y_continuous(limits=c(-60, 60))

ggsave(multivar.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mouse-liver_multiivariateLM.png",
       height=4.75, width=6.75, dpi=300)

# plot the uni- and multivariate on the same graph?
liver.rlm.df$Model <- "Univariate"
liver.rlm.res$Model <- "Multivariate"

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
  scale_y_continuous(limits=c(-60, 60))

ggsave(all.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mouse-liver_allLM.png",
       height=5.6, width=7.75, dpi=300)

###############################
## Plot CpG islands vs alpha ##
###############################
liver.match$CpGisland <- factor(liver.match$N_CpG,
                               levels=c(0, 1),
                               labels=c("Non-CpG island", "CpG island"))
cpg.cols <- c("#027E00", "#00DDEC")
names(cpg.cols) <- levels(liver.match$CpGisland)

tc_cpg <- ggplot(liver.match, aes(x=CpGisland, y=Residual.CV2, colour=CpGisland)) + 
  theme_mike() + 
  geom_jitter(position=position_jitterdodge(jitter.height=0,
                                            jitter.width=1.5),
              alpha=0.7) +
  geom_boxplot(width=0.5, fill='white', colour='black') +
  scale_colour_manual(values=cpg.cols) +
  labs(x="Overlapping CpG island", y=expression(paste("Residual CV"^2))) +
  guides(colour=FALSE)

ggsave(tc_cpg,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mouse-liver_boxplot_CpGislands-overdispersion.png",
       height=3.75, width=4.75, dpi=300)

tc_cpgMe <- ggplot(liver.match, aes(x=CpGisland, y=Mean, colour=CpGisland)) + 
  theme_mike() + 
  geom_jitter(position=position_jitterdodge(jitter.height=0,
                                            jitter.width=1.5),
              alpha=0.7) +
  geom_boxplot(width=0.5, fill='white', colour='black') +
  scale_colour_manual(values=cpg.cols) +
  labs(x="Overlapping CpG island", y=expression(paste("Mean log"[2], " Expression"))) +
  guides(colour=FALSE)

ggsave(tc_cpgMe,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mouse-liver_boxplot_CpGislands-mean.png",
       height=3.75, width=4.75, dpi=300)
source("~/Dropbox/R_sessions/Noise/mESC_chromHMM.R")
source("~/Dropbox/R_sessions/Noise/genomic_noise_features.R")
source("~/Dropbox/R_sessions/Noise/tcell_genomic_noise_features.R")
library(MASS)
library(ggplot2)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

tcell.vars <- colnames(genomic.features)[2:16]
tcell.vars <- tcell.vars[!grepl(tcell.vars, pattern="(NMI)|(CGI_SIZE)|(cpg_)|(PHAST)")]

tcell.genomic.vars <- paste(c("Mean",
                             tcell.vars),
                           collapse=" + ")

tcell.match <- merge(tcell.gene.summary, genomic.features,
                    by='GENE')
#########################################
## univariate robust linear regression ##
#########################################
# this is plotting for presentation/manuscript figures
tcell.univariate_list <- list()

tcell.var.names <- unlist(strsplit(tcell.genomic.vars, split=" + ", fixed=T))
# remove redundant variables, i.e CpG island AND NMI
# remove CpG island characteristics
tcell.var.names <- tcell.var.names[!grepl(tcell.var.names, pattern="(NMI)|(CGI_SIZE)|(cpg_)|(PHAST)")]

for(x in seq_along(tcell.var.names)){
  .variable <- paste(tcell.var.names[x], sep=" + ")
  .glm.form <- as.formula(paste("Residual.CV2", .variable, sep=" ~ "))
  
  m.rlm <- rlm(.glm.form, data=tcell.match)
  m.robust <- summary(m.rlm)
  m.rlm.res <- as.data.frame(m.robust$coefficients)
  m.rlm.res$Pval <- 2*pt(-abs(m.rlm.res[, 3]), df=3)
  m.rlm.res$Sig <- as.numeric(m.rlm.res$Pval <= 0.05)
  m.res.mat <- as.matrix(m.rlm.res)
  m.res.var <- m.res.mat[2, ]
  names(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
  tcell.univariate_list[[tcell.var.names[x]]] <- as.list(m.res.var)
}

tcell.rlm.df <- do.call(rbind.data.frame, tcell.univariate_list)
tcell.rlm.df$Predictor <- rownames(tcell.rlm.df)
tcell.rlm.df$Direction <- "NoEffect"
tcell.rlm.df$Direction[tcell.rlm.df$COEFF < 0 & tcell.rlm.df$Sig == 1] <- "Less"
tcell.rlm.df$Direction[tcell.rlm.df$COEFF > 0 & tcell.rlm.df$Sig == 1] <- "More"

# give the features more informative/better formated names
tcell.rlm.df$Predictor[tcell.rlm.df$Predictor == "Mean"] <- "Mean expression"
tcell.rlm.df$Predictor[tcell.rlm.df$Predictor == "N_CpG"] <- "CpG island"
tcell.rlm.df$Predictor[tcell.rlm.df$Predictor == "GC"] <- "Promoter GC %"
tcell.rlm.df$Predictor[tcell.rlm.df$Predictor == "SP1"] <- "Number SP1 motifs"
tcell.rlm.df$Predictor[tcell.rlm.df$Predictor == "TBP"] <- "Number TBP motifs"
tcell.rlm.df$Predictor[tcell.rlm.df$Predictor == "MED_PHAST"] <- "Median PhastCons"
tcell.rlm.df$Predictor[tcell.rlm.df$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
tcell.rlm.df$Predictor[tcell.rlm.df$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
tcell.rlm.df$Predictor[tcell.rlm.df$Predictor == "NALIGN_PHAST"] <- "Numer aligned bases"
tcell.rlm.df$Predictor[tcell.rlm.df$Predictor == "EXON_AVLENGTH"] <- "Mean exon length"
tcell.rlm.df$Predictor[tcell.rlm.df$Predictor == "EXON_COUNT"] <- "Number of exons"
tcell.rlm.df$Predictor[tcell.rlm.df$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"
tcell.rlm.df$Predictor[tcell.rlm.df$Predictor == "EXON_VARLENGTH"] <- "Exon length variance"

tcell.rlm.df$Tissue <- "Tcell"
tcell.rlm.df$Species <- "Mouse"

univar.plot <- ggplot(tcell.rlm.df,
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
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/Tcell_univariateLM.png",
       height=4.75, width=6.75, dpi=300)


###########################################
## multivariate robust linear regression ##
###########################################

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
tcell.rlm.res$Predictor[tcell.rlm.res$Predictor == "N_CpG"] <- "CpG island"
tcell.rlm.res$Predictor[tcell.rlm.res$Predictor == "GC"] <- "Promoter GC %"
tcell.rlm.res$Predictor[tcell.rlm.res$Predictor == "SP1"] <- "Number SP1 motifs"
tcell.rlm.res$Predictor[tcell.rlm.res$Predictor == "TBP"] <- "Number TBP motifs"
tcell.rlm.res$Predictor[tcell.rlm.res$Predictor == "MED_PHAST"] <- "Median PhastCons"
tcell.rlm.res$Predictor[tcell.rlm.res$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
tcell.rlm.res$Predictor[tcell.rlm.res$Predictor == "SUM_PHAST"] <- "Total promoter PhastCons"
tcell.rlm.res$Predictor[tcell.rlm.res$Predictor == "NALIGN_PHAST"] <- "Numer aligned bases"
tcell.rlm.res$Predictor[tcell.rlm.res$Predictor == "EXON_AVLENGTH"] <- "Mean exon length"
tcell.rlm.res$Predictor[tcell.rlm.res$Predictor == "EXON_COUNT"] <- "Number of exons"
tcell.rlm.res$Predictor[tcell.rlm.res$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"
tcell.rlm.res$Predictor[tcell.rlm.res$Predictor == "EXON_VARLENGTH"] <- "Exon length variance"
tcell.rlm.res$Predictor[tcell.rlm.res$Predictor == "Mean"] <- "Mean expression"

# remove the intercept term
tcell.rlm.res <- tcell.rlm.res[!grepl(tcell.rlm.res$Predictor, pattern="Intercept"), ]

tcell.rlm.res$Tissue <- "Tcell"
tcell.rlm.res$Species <- "Mouse"

write.table(tcell.rlm.res,
            file="~/Dropbox/Noise_genomics/Model_results/mouse_Tcell_multivariateRLM.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")

multivar.plot <- ggplot(tcell.rlm.res,
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
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/Tcell_multiivariateLM.png",
       height=4.75, width=6.75, dpi=300)

# plot the uni- and multivariate on the same graph?
tcell.rlm.df$Model <- "Univariate"
tcell.rlm.res$Model <- "Multivariate"

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
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
        axis.text=element_text(size=16)) +
  labs(x="Annotation", y="t-statistic") +
  guides(fill=FALSE, shape=FALSE) +
  scale_y_continuous(limits=c(-50, 50))

ggsave(all.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/Tcell_allLM.png",
       height=5.6, width=7.75, dpi=300)

###############################
## Plot CpG islands vs alpha ##
###############################
tcell.match$CpGisland <- factor(tcell.match$N_CpG,
                                levels=c(0, 1),
                                labels=c("Non-CpG island", "CpG island"))
cpg.cols <- c("#027E00", "#00DDEC")
names(cpg.cols) <- levels(tcell.match$CpGisland)

tc_cpg <- ggplot(tcell.match, aes(x=CpGisland, y=Residual.CV2, colour=CpGisland)) + 
  theme_mike() + 
  geom_jitter(position=position_jitterdodge(jitter.height=0,
                                            jitter.width=1.5),
              alpha=0.7) +
  geom_boxplot(width=0.5, fill='white', colour='black') +
  scale_colour_manual(values=cpg.cols) +
  labs(x="Overlapping CpG island", y=expression(paste("Residual CV"^2))) +
  guides(colour=FALSE)  +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16))

ggsave(tc_cpg,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/Tcell_boxplot_CpGislands-overdispersion.png",
       height=3.75, width=4.75, dpi=300)

tc_cpgMe <- ggplot(tcell.match, aes(x=CpGisland, y=Mean, colour=CpGisland)) + 
  theme_mike() + 
  geom_jitter(position=position_jitterdodge(jitter.height=0,
                                            jitter.width=1.5),
              alpha=0.7) +
  geom_boxplot(width=0.5, fill='white', colour='black') +
  scale_colour_manual(values=cpg.cols) +
  labs(x="Overlapping CpG island", y=expression(paste("Mean log"[2], " Expression"))) +
  guides(colour=FALSE)

ggsave(tc_cpgMe,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/Tcell_boxplot_CpGislands-mean.png",
       height=3.75, width=4.75, dpi=300)

## Investigating genomic factors that influence gene expression noise
# test whether there is any relationship between mean and other variables, with respect to
# explanatory capacity on CV^2
# for each feature fit the model with and without the feature, conditional on the mean and do an LR test
# to get a model likelihood I'll need to switch to a linear model as robust LM doesn't
# evaluate a specific likelihood function.
# Might be achievable using the robustbase package

source("~/Dropbox/R_sessions/Noise/mESC_chromHMM.R")
source("~/Dropbox/R_sessions/Noise/genomic_noise_features.R")
source("~/Dropbox/R_sessions/Noise/mESC_genomic_noise_features-spline_fit.R")

library(ggplot2)
#library(MASS)
library(robustbase)
library(scales)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

# calculate 1/mean
mesc.gene.summary$Reciprocal.Mean <- 1/mesc.gene.summary$Mean

mesc.vars <- colnames(genomic.features)[2:16]
mesc.vars <- mesc.vars[!grepl(mesc.vars, pattern="(NMI)|(CGI_SIZE)|(cpg_)|(PHAST)|(Mean)")]

mesc.genomic.vars <- paste(c("Reciprocal.Mean",
                             mesc.vars),
                           collapse=" + ")

mesc.match <- merge(mesc.gene.summary, genomic.features,
                    by='GENE')
#########################################
## univariate robust linear regression ##
#########################################
# this is plotting for presentation/manuscript figures
mesc.univariate_list <- list()

mesc.var.names <- unlist(strsplit(mesc.genomic.vars, split=" + ", fixed=T))
# remove redundant variables, i.e CpG island AND NMI
# remove CpG island characteristics
mesc.var.names <- mesc.var.names[!grepl(mesc.var.names, pattern="(NMI)|(CGI_SIZE)|(cpg_)|(PHAST)|(Mean)")]

null.glm.form <- as.formula(paste("CV2", "Spline.Fit", sep=" ~ "))
# set a higher number of iterations, as sometimes doesn't converge
model.control <- lmrob.control(max.it=500, k.max=500, rel.tol=1e-6)
null.model <- lmrob(null.glm.form, data=mesc.match, control=model.control)
  
for(x in seq_along(mesc.var.names)){
  .variable <- paste("Spline.Fit", mesc.var.names[x], sep=" + ")
  .glm.form <- as.formula(paste("CV2", .variable, sep=" ~ "))
  
  m.rlm <- lmrob(.glm.form, data=mesc.match, control=model.control)
  # get the table of deviances (log-likelihoods), test by chi-squared
  m.lrtest <- anova(m.rlm, null.model, test="Deviance")
  m.rlm.res <- list("Pval"=m.lrtest$`Pr(>chisq)`[2],
                    "PseudoDF"=m.lrtest$pseudoDf[2],
                    "Stat"=m.lrtest$Test.Stat[2],
                    "DF"=m.lrtest$Df[2],
                    "COEFF"=summary(m.rlm)$coefficients[3, 1])
  
  m.rlm.res$Sig <- as.numeric(m.rlm.res$Pval <= 0.05)
  m.res.mat <- as.matrix(m.rlm.res)
  mesc.univariate_list[[mesc.var.names[x]]] <- as.list(m.rlm.res)
}

mesc.rlm.df <- do.call(rbind.data.frame, mesc.univariate_list)
mesc.rlm.df$Predictor <- rownames(mesc.rlm.df)
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
                      aes(x=reorder(Predictor, Stat),
                          y=Stat, fill=Direction)) +
  geom_point(alpha=0.55, shape=21, size=5) + theme_mike() +
  scale_fill_manual(values=c("#62148f", "#feaf10", "#878787")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Annotation", y="t-statistic") +
  guides(fill=FALSE) +
  geom_hline(mapping=aes(yintercept=0), linetype="dashed", colour="grey") +
  scale_y_continuous(limits=c(0, 25), oob=squish)

univar.plot

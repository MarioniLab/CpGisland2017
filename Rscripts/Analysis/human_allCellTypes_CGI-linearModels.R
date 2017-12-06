# try running the multivariate model, combining all mouse data together, blocking on cell type
# test for interactions between cell type and genomic features
# this should improve power, and show if any of the CpG island features are cell-type specific

source("~/Dropbox/R_sessions/Noise/human_genomic_noise_features.R")
source("~/Dropbox/R_sessions/Noise/human_beta_islet_genomic_noise_features.R")
source("~/Dropbox/R_sessions/Noise/human_pancreas_genomic_noise_features.R")
source("~/Dropbox/R_sessions/Noise/hESC_genomic_noise_features.R")

library(MASS)
library(ggplot2)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

genomic.features$CGI_SIZE.kb <- genomic.features$CGI_SIZE/1000

## setup panc features

panc.vars <- colnames(human.genomic.features)
panc.vars <- panc.vars[grepl(panc.vars, pattern="(CGI_SIZE.kb)|(cpg_[O|R])|(TOTLEN)|
                               (SP1)")]

panc.genomic.vars <- paste(c("Mean",
                              panc.vars),
                            collapse=" + ")

panc.match <- merge(panc.gene.summary, human.genomic.features,
                     by='GENE')
panc.match <- panc.match[panc.match$N_CpG == 1, ]


# setup T cell features
beta.vars <- colnames(human.genomic.features)
beta.vars <- beta.vars[grepl(beta.vars, pattern="(CGI_SIZE.kb)|(cpg_[O|R])|(TOTLEN)|
                               (SP1)")]

beta.genomic.vars <- paste(c("Mean",
                              beta.vars),
                            collapse=" + ")

beta.match <- merge(beta.gene.summary, human.genomic.features,
                     by='GENE')

# select only those genes with a CpG island
beta.match <- beta.match[beta.match$N_CpG == 1, ]

# setup hesc features

hesc.vars <- colnames(human.genomic.features)
hesc.vars <- hesc.vars[grepl(hesc.vars, pattern="(CGI_SIZE.kb)|(cpg_[O|R])|(TOTLEN)|
                             (SP1)")]

hesc.genomic.vars <- paste(c("Mean",
                             hesc.vars),
                           collapse=" + ")

hesc.match <- merge(hesc.gene.summary, human.genomic.features,
                    by='GENE')

hesc.match <- hesc.match[hesc.match$N_CpG == 1, ]

# make sure all columns are the same, if not select all the same columns
panc.match$CellType <- "AlphaIslet"
beta.match$CellType <- "BetaIslet"
hesc.match$CellType <- "hESC"
comm.cols <- intersect(colnames(hesc.match), intersect(colnames(beta.match), colnames(panc.match)))

all.match <- do.call(rbind.data.frame, list("hesc"=hesc.match[, comm.cols],
                                            "beta"=beta.match[, comm.cols],
                                            "panc"=panc.match[, comm.cols]))

###########################################
## multivariate robust linear regression ##
###########################################
all.vars <- unique(c(hesc.vars, beta.vars, panc.vars))
all.vars <- c(all.vars, "CellType")

all.genomic.vars <- paste(c("Mean",
                            all.vars),
                          collapse=" + ")

# add an interaction term between cell type and CpG island features


all.glm.form <- as.formula(paste("Residual.CV2",
                                 all.genomic.vars, sep=" ~ "))

all.rlm <- rlm(all.glm.form, data=all.match)
all.robust <- summary(all.rlm)
all.rlm.res <- as.data.frame(all.robust$coefficients)
all.rlm.res$Pval <- 2*pt(-abs(all.rlm.res[, 3]), df=dim(all.match)[2]-1)
all.rlm.res$Sig <- as.numeric(all.rlm.res$Pval <= 0.05)
all.rlm.res$Predictor <- rownames(all.rlm.res)
colnames(all.rlm.res) <- c("COEFF", "SE", "STAT", "P", "Sig", "Predictor")
all.rlm.res$Direction <- "NoEffect"
all.rlm.res$Direction[all.rlm.res$COEFF < 0 & all.rlm.res$Sig == 1] <- "Less"
all.rlm.res$Direction[all.rlm.res$COEFF > 0 & all.rlm.res$Sig == 1] <- "More"

# give the features more informative/better formated names
all.rlm.res$Predictor[all.rlm.res$Predictor == "Mean"] <- "Mean expression"
all.rlm.res$Predictor[all.rlm.res$Predictor == "N_CpG"] <- "CpG island"
all.rlm.res$Predictor[all.rlm.res$Predictor == "cpg_GCNUM"] <- "CpG island GC content"
all.rlm.res$Predictor[all.rlm.res$Predictor == "SP1"] <- "Number SP1 motifs"
all.rlm.res$Predictor[all.rlm.res$Predictor == "cpg_Overlap"] <- "CpG island overlap"
all.rlm.res$Predictor[all.rlm.res$Predictor == "cpg_RATIO"] <- "CpG dinucleotide ratio"
all.rlm.res$Predictor[all.rlm.res$Predictor == "CGI_SIZE.kb"] <- "CpG island size"
all.rlm.res$Predictor[all.rlm.res$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"

all.rlm.res <- all.rlm.res[!grepl(all.rlm.res$Predictor, pattern="Intercept"), ]
all.rlm.res$Tissue <- "all"
all.rlm.res$Species <- "Mouse"

write.table(all.rlm.res,
            file="~/Dropbox/Noise_genomics/Model_results/human_All_CGI-multivariateRLM.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")
effect.cols <- c("#62148f", "#878787", "#feaf10")
names(effect.cols) <- c("Less", "NoEffect", "More")

cpg.plot <- ggplot(all.rlm.res,
                   aes(x=reorder(Predictor, -STAT),
                       y=STAT, fill=Direction)) +
  geom_point(alpha=0.55, shape=21, size=5) + theme_mike() +
  scale_fill_manual(values=effect.cols) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Annotation", y="t-statistic") +
  guides(fill=FALSE) +
  geom_hline(mapping=aes(yintercept=0), linetype="dashed", colour="grey")

ggsave(cpg.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-All_CGI-multivariateLM.png",
       height=4.75, width=6.75, dpi=300)


##########################################################
## test for interactions between cell type and features ##
##########################################################
int.vars <- unique(c(hesc.vars, beta.vars, panc.vars))
int.vars <- paste(int.vars, "CellType", sep="*")

int.genomic.vars <- paste(c("Mean",
                            int.vars),
                          collapse=" + ")

# add an interaction term between cell type and CpG island features


int.glm.form <- as.formula(paste("Residual.CV2",
                                 int.genomic.vars, sep=" ~ "))

int.rlm <- rlm(int.glm.form, data=all.match)
int.robust <- summary(int.rlm)
int.rlm.res <- as.data.frame(int.robust$coefficients)
int.rlm.res$Pval <- 2*pt(-abs(int.rlm.res[, 3]), df=dim(all.match)[2]-1)
int.rlm.res$Sig <- as.numeric(int.rlm.res$Pval <= 0.05)
int.rlm.res$Predictor <- rownames(int.rlm.res)
colnames(int.rlm.res) <- c("COEFF", "SE", "STAT", "P", "Sig", "Predictor")
int.rlm.res$Direction <- "NoEffect"
int.rlm.res$Direction[int.rlm.res$COEFF < 0 & int.rlm.res$Sig == 1] <- "Less"
int.rlm.res$Direction[int.rlm.res$COEFF > 0 & int.rlm.res$Sig == 1] <- "More"

# give the features more informative/better formated names
int.rlm.res$Predictor[int.rlm.res$Predictor == "Mean"] <- "Mean expression"
int.rlm.res$Predictor[int.rlm.res$Predictor == "N_CpG"] <- "CpG island"
int.rlm.res$Predictor[int.rlm.res$Predictor == "cpg_GCNUM"] <- "CpG island GC content"
int.rlm.res$Predictor[int.rlm.res$Predictor == "SP1"] <- "Number SP1 motifs"
int.rlm.res$Predictor[int.rlm.res$Predictor == "cpg_Overlap"] <- "CpG island overlap"
int.rlm.res$Predictor[int.rlm.res$Predictor == "cpg_RATIO"] <- "CpG dinucleotide ratio"
int.rlm.res$Predictor[int.rlm.res$Predictor == "CGI_SIZE.kb"] <- "CpG island size"
int.rlm.res$Predictor[int.rlm.res$Predictor == "EXON_TOTLENGTH"] <- "Transcript length"

int.rlm.res <- int.rlm.res[!grepl(int.rlm.res$Predictor, pattern="Intercept"), ]

write.table(int.rlm.res,
            file="~/Dropbox/Noise_genomics/Model_results/human_Interaction_CGI-multivariateRLM.tsv",
            quote=FALSE, row.names=FALSE, sep="\t")
effect.cols <- c("#62148f", "#878787", "#feaf10")
names(effect.cols) <- c("Less", "NoEffect", "More")

cpg.plot <- ggplot(int.rlm.res,
                   aes(x=reorder(Predictor, -STAT),
                       y=STAT, fill=Direction)) +
  geom_point(alpha=0.55, shape=21, size=5) + theme_mike() +
  scale_fill_manual(values=effect.cols) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Annotation", y="t-statistic") +
  guides(fill=FALSE) +
  geom_hline(mapping=aes(yintercept=0), linetype="dashed", colour="grey")

ggsave(cpg.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/human-Interaction_CGI-multivariateLM.png",
       height=6.75, width=6.75, dpi=300)
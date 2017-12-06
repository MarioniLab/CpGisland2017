library(ggplot2)
library(scales)
library(data.table)
library(cowplot)
source("~/Dropbox/R_sessions/Noise/mESC_chromHMM.R")
source("~/Dropbox/R_sessions/Noise/genomic_noise_features.R")
source("~/Dropbox/R_sessions/Noise/mESC_genomic_noise_features.R")
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

extra.chip <- read.table("~/Dropbox/Chromatin/ESC/GSE36114-Window_ChIP-mean.tsv.gz",
                         h=T, stringsAsFactors=F, sep="\t")
mesc.match <- merge(mesc.gene.summary, genomic.features,
                    by='GENE')

# mesc.consist <- merge(mesc.consist, genomic.features, by='GENE')
mesc.merge <- merge(mesc.match, extra.chip, by="GENE")
mesc.merge$CGI_SIZE.kb <- mesc.merge$CGI_SIZE/1000
mesc.merge$EXON_TOTLENGTH <- log10(mesc.merge$EXON_TOTLENGTH)
mesc.merge$Recip.Mean <- 1/mesc.merge$Mean

mesc.merge$GeneGroup <- 2
mesc.merge$GeneGroup[mesc.merge$GENE %in% mouse.consistent.genes] <- 3
mesc.merge$GeneGroup[mesc.merge$GENE %in% mouse.tra.genes] <- 1
mesc.merge$GeneGroup <- factor(mesc.merge$GeneGroup,
                               levels=c(1, 2, 3),
                               labels=c("TissueSpec", "Misc", "Consistent"))

# split NMIs into size bins on quartiles
mesc.merge$CGI_SIZE.group <- as.character(cut(mesc.merge$CGI_SIZE.kb,
                                              breaks=quantile(mesc.merge$CGI_SIZE.kb[mesc.merge$CGI_SIZE.kb != 0],
                                                              probs=c(0, 0.25, 0.5, 0.75, 1.0))))
mesc.merge$CGI_SIZE.group[(mesc.merge$CGI_SIZE.kb <= 0.201)] <- "(0,0.201]"
mesc.merge$CGI_SIZE.group[is.na(mesc.merge$CGI_SIZE.group) |mesc.merge$CGI_SIZE.kb == 0] <- "Absent"

# merge the V.Short and Short CGIs
mesc.merge$CGI_SIZE.group[mesc.merge$CGI_SIZE.group == "(0,0.201]"] <- "(0.201,0.422]"
mesc.merge$CGI_SIZE.group <- factor(mesc.merge$CGI_SIZE.group,
                                    labels=c("Absent", "Short", "Short.Mid", "Long.Mid", "Long"),
                                    levels=c("Absent", "(0.201,0.422]", "(0.422,0.627]", "(0.627,0.904]", "(0.904,5.13]"))

########################################
## Define bivalent promoters as those ##
## with both H3K27me3 and H3K4me3     ##
########################################

# define active, repressed and poised enhancers based on H3K27me3 and H3K4me3
minimums <- function(x) which(x - shift(x, 1) < 5e-5  & x - shift(x, 1, type='lead') < 5e-5)

# find the local minimum point of the kernel density, with default bandwidth
# promoters below are inactive and above are active
mesc.merge$ChIP_Ratio <- mesc.merge$H3K27me3/mesc.merge$H3K4me3

ratio.min.value <- density(mesc.merge$ChIP_Ratio)$x[(minimums(density(mesc.merge$ChIP_Ratio)$y)[1])]
k4me3.min.value <- density(mesc.merge$H3K4me3)$x[(minimums(density(mesc.merge$H3K4me3)$y)[1])]
k27me3.min.value <- mean(mesc.merge$H3K27me3)

mesc.merge$ChIP_Ratio.bin <- as.numeric(mesc.merge$ChIP_Ratio >= ratio.min.value)
mesc.merge$H3K4me3.bin <- as.numeric(mesc.merge$H3K4me3 >= k4me3.min.value)
mesc.merge$H3K27me3.bin <- as.numeric(mesc.merge$H3K27me3 >= k27me3.min.value)

k27.cat <- ggplot(mesc.merge,
                 aes(x=H3K27me3, colour=as.factor(H3K27me3.bin))) +
  geom_vline(mapping=aes(xintercept=k27me3.min.value),
             linetype="dashed", colour="grey", size=2) +
  geom_density() + theme_mike() +
  scale_colour_manual(values=c("darkblue", "red"), labels=c("NotRepressed", "Repressed")) +
  guides(colour=guide_legend(title="Category")) +
  labs(x="H3K27me3 ChIP (Signa/Input)", y="Density")

ggsave(k27.cat,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mESC_H3K27me3_categorization.png",
       width=4.75, height=3.75, dpi=300)

k4.cat <- ggplot(mesc.merge,
                 aes(x=H3K4me3, colour=as.factor(H3K4me3.bin))) +
  geom_vline(mapping=aes(xintercept=k4me3.min.value),
             linetype="dashed", colour="grey", size=2) +
  geom_density() + theme_mike() +
  scale_colour_manual(values=c("darkblue", "red"), labels=c("NotActive", "Active")) +
  guides(colour=guide_legend(title="Category")) +
  labs(x="H3K4me3 ChIP (Signa/Input)", y="Density")

ggsave(k4.cat,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mESC_H3K4me3_categorization.png",
       width=4.75, height=3.75, dpi=300)


mesc.merge$PromoterActivity <- "Unk"
mesc.merge$PromoterActivity[(mesc.merge$H3K27me3.bin == 0) &
                              (mesc.merge$H3K4me3.bin == 1)] <- "H3K4me3"

mesc.merge$PromoterActivity[(mesc.merge$H3K27me3.bin == 1) &
                              (mesc.merge$H3K4me3.bin == 1)] <- "Bivalent"

mesc.merge$PromoterActivity[(mesc.merge$H3K27me3.bin == 1) &
                              (mesc.merge$H3K4me3.bin == 0)] <- "H3K27me3"

mesc.merge$PromoterActivity[(mesc.merge$H3K27me3.bin == 0) &
                              (mesc.merge$H3K4me3.bin == 0)] <- "Bivalent"

mesc.merge$PromoterActivity <- factor(mesc.merge$PromoterActivity,
                                      labels=c("H3K27me3", "Bivalent", "H3K4me3"),
                                      levels=c("H3K27me3", "Bivalent", "H3K4me3"))

chip.boundary <- ggplot(mesc.merge, aes(x=H3K27me3, y=H3K4me3, colour=PromoterActivity)) +
  geom_point(alpha=0.4) + theme_mike() +
  scale_colour_manual(values=c("#a614ef", "#11d101", "#ff7306",  "#f0db87")) +
  guides(colour=guide_legend(title='Promoter State'))

ggsave(chip.boundary,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mESC_ChIP_boundary-scatter.png",
       height=3.75, width=3.75*1.618, dpi=300)

#########################################
## univariate robust linear regression ##
#########################################
# this is plotting for presentation/manuscript figures
mesc.vars <- c("H3K27me3", "H3K4me3", "cpg_RATIO", "CGI_SIZE.kb")
mesc.univariate_list <- list()

for(x in seq_along(mesc.vars)){
  .variable <- paste(c("Recip.Mean", mesc.vars[x]), collapse=" + ")
  .glm.form <- as.formula(paste("Residual.CV2", .variable, sep=" ~ "))
  
  m.rlm <- rlm(.glm.form, data=mesc.merge)
  m.robust <- summary(m.rlm)
  m.rlm.res <- as.data.frame(m.robust$coefficients)
  m.rlm.res$Pval <- 2*pt(-abs(m.rlm.res[, 3]), df=3)
  m.rlm.res$Sig <- as.numeric(m.rlm.res$Pval <= 0.05)
  m.res.mat <- as.matrix(m.rlm.res)
  m.res.var <- m.res.mat[3, ]
  names(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
  mesc.univariate_list[[mesc.vars[x]]] <- as.list(m.res.var)
}

# add the mean expression on it's own
.glm.form <- as.formula(paste("Residual.CV2", "Recip.Mean", sep=" ~ "))

m.rlm <- rlm(.glm.form, data=mesc.merge)
m.robust <- summary(m.rlm)
m.rlm.res <- as.data.frame(m.robust$coefficients)
m.rlm.res$Pval <- 2*pt(-abs(m.rlm.res[, 3]), df=3)
m.rlm.res$Sig <- as.numeric(m.rlm.res$Pval <= 0.05)
m.res.mat <- as.matrix(m.rlm.res)
m.res.var <- m.res.mat[2, ]
names(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
mesc.univariate_list[["Recip.Mean"]] <- as.list(m.res.var)

mesc.rlm.df <- do.call(rbind.data.frame, mesc.univariate_list)
mesc.rlm.df$Predictor <- rownames(mesc.rlm.df)
mesc.rlm.df$Direction <- "NoEffect"
mesc.rlm.df$Direction[mesc.rlm.df$COEFF < 0 & mesc.rlm.df$Sig == 1] <- "Less"
mesc.rlm.df$Direction[mesc.rlm.df$COEFF > 0 & mesc.rlm.df$Sig == 1] <- "More"

# give the features more informative/better formated names
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "Recip.Mean"] <- "Mean expression"
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "cpg_GCNUM"] <- "CpG island GC content"
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "cpg_RATIO"] <- "CpG dinucleotide ratio"
mesc.rlm.df$Predictor[mesc.rlm.df$Predictor == "CGI_SIZE.kb"] <- "CpG island size"
mesc.rlm.df$Model <- "Univariate"

write.table(mesc.rlm.df,
            file="~/Dropbox/Noise_genomics/Model_results/mouse_ESC_chromatin-UnivariateLM.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

### run the RLM with continuous variables
cv2.fit <- rlm(Residual.CV2 ~  Recip.Mean + H3K27me3 + H3K4me3 + cpg_RATIO + CGI_SIZE.kb,
               data=mesc.merge)
cv2.robust <- summary(cv2.fit)
cv2.rlm.res <- as.data.frame(cv2.robust$coefficients)
n.params <- dim(cv2.rlm.res)[1]
cv2.rlm.res$Pval <- 2*pt(-abs(cv2.rlm.res[, 3]), df=(dim(mesc.merge)[1]-1)-n.params)
cv2.rlm.res$Sig <- as.numeric(cv2.rlm.res$Pval <= 0.01)
cv2.rlm.res$Predictor <- rownames(cv2.rlm.res)
cv2.rlm.res <- cv2.rlm.res[-1, ]
colnames(cv2.rlm.res) <- c("COEFF", "SE", "STAT", "P", "Sig", "Predictor")

cv2.rlm.res$Direction <- "NoEffect"
cv2.rlm.res$Direction[cv2.rlm.res$COEFF < 0 & cv2.rlm.res$Sig == 1] <- "Less"
cv2.rlm.res$Direction[cv2.rlm.res$COEFF > 0 & cv2.rlm.res$Sig == 1] <- "More"

# give the features more informative/better formated names
cv2.rlm.res$Predictor[cv2.rlm.res$Predictor == "Recip.Mean"] <- "Mean expression"
cv2.rlm.res$Predictor[cv2.rlm.res$Predictor == "cpg_GCNUM"] <- "CpG island GC content"
cv2.rlm.res$Predictor[cv2.rlm.res$Predictor == "cpg_RATIO"] <- "CpG dinucleotide ratio"
cv2.rlm.res$Predictor[cv2.rlm.res$Predictor == "CGI_SIZE.kb"] <- "CpG island size"
cv2.rlm.res$Model <- "Multivariate"

write.table(cv2.rlm.res,
            file="~/Dropbox/Noise_genomics/Model_results/mouse_ESC_chromatin-MultivariateLM.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

effect.cols <- c("#62148f", "#878787", "#feaf10")
names(effect.cols) <- c("Less", "NoEffect", "More")

cpg.plot <- ggplot(cv2.rlm.res,
                   aes(x=reorder(Predictor, -STAT),
                       y=STAT, fill=Direction)) +
  geom_point(alpha=0.55, shape=21, size=5) + theme_mike() +
  scale_fill_manual(values=effect.cols) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Annotation", y="t-statistic") +
  guides(fill=FALSE) +
  scale_y_continuous(limits=c(-30, 30), oob=squish) +
  geom_hline(mapping=aes(yintercept=0), linetype="dashed", colour="grey") 

ggsave(cpg.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mouse_ESC_chromatin-multivariateLM.png",
       height=4.75, width=5.75, dpi=300)


# plot uni and multivariate results together
chip.rlm.all <- do.call(rbind.data.frame, list("Univariate"=mesc.rlm.df,
                                              "Multivariate"=cv2.rlm.res))

chip.all.plot <- ggplot(chip.rlm.all,
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
  scale_y_continuous(limits=c(-30, 30), oob=squish)

ggsave(chip.all.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mouse_ESC_chromatin-allLM.png",
       height=4.75, width=5.75, dpi=300)


###############################################
# plot uni and multivariate results together ##
###############################################
mesc.rlm.all <- do.call(rbind.data.frame, list("Univariate"=mesc.rlm.df,
                                               "Multivariate"=cv2.rlm.res))

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
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Annotation", y="t-statistic") +
  guides(fill=FALSE, shape=FALSE) +
  scale_y_continuous(limits=c(-30, 30), oob=squish)

ggsave(all.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mouse_ESC_chromatin-allLM.png",
       height=4.75, width=5.75, dpi=300)

### plot relationships

rcv.k27 <- ggplot(mesc.merge,
                  aes(x=H3K27me3, y=Residual.CV2)) +
  geom_point(alpha=0.5) + theme_mike() +
  labs(x="H3K27me3 ChIP signal", y=expression(paste("Residual CV"^2)))

rcv.k4 <- ggplot(mesc.merge,
                  aes(x=H3K4me3, y=Residual.CV2)) +
  geom_point(alpha=0.5) + theme_mike() +
  labs(x="H3K4me3 ChIP signal", y=expression(paste("Residual CV"^2)))

chip.rcv <- plot_grid(rcv.k27, rcv.k4, nrow=2)

ggsave(chip.rcv,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mESC_ChIP_rCV.png",
       height=4.75, width=4.75, dpi=300)

########################################
### run RLM using 'promoter activity' ##
########################################

#########################################
## univariate robust linear regression ##
#########################################
# this is plotting for presentation/manuscript figures
biv.vars <- c("cpg_RATIO", "CGI_SIZE.kb")
biv.univariate_list <- list()

for(x in seq_along(biv.vars)){
  .variable <- paste(c("Recip.Mean", biv.vars[x]), collapse=" + ")
  .glm.form <- as.formula(paste("Residual.CV2", .variable, sep=" ~ "))
  
  m.rlm <- rlm(.glm.form, data=mesc.merge)
  m.robust <- summary(m.rlm)
  m.rlm.res <- as.data.frame(m.robust$coefficients)
  m.rlm.res$Pval <- 2*pt(-abs(m.rlm.res[, 3]), df=3)
  m.rlm.res$Sig <- as.numeric(m.rlm.res$Pval <= 0.05)
  m.res.mat <- as.matrix(m.rlm.res)
  m.res.var <- m.res.mat[3:dim(m.res.mat)[1], ]
  names(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
  biv.univariate_list[[biv.vars[x]]] <- as.list(m.res.var)
}

# add the mean expression on it's own
.glm.form <- as.formula(paste("Residual.CV2", "Recip.Mean", sep=" ~ "))

m.rlm <- rlm(.glm.form, data=mesc.merge)
m.robust <- summary(m.rlm)
m.rlm.res <- as.data.frame(m.robust$coefficients)
m.rlm.res$Pval <- 2*pt(-abs(m.rlm.res[, 3]), df=3)
m.rlm.res$Sig <- as.numeric(m.rlm.res$Pval <= 0.05)
m.res.mat <- as.matrix(m.rlm.res)
m.res.var <- m.res.mat[2, ]
names(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
biv.univariate_list[["Recip.Mean"]] <- as.list(m.res.var)

# add the promoter activity
.variable <- paste(c("Recip.Mean", "PromoterActivity"), collapse=" + ")
.glm.form <- as.formula(paste("Residual.CV2", .variable, sep=" ~ "))

m.rlm <- rlm(.glm.form, data=mesc.merge)
m.robust <- summary(m.rlm)
m.rlm.res <- as.data.frame(m.robust$coefficients)
m.rlm.res$Pval <- 2*pt(-abs(m.rlm.res[, 3]), df=3)
m.rlm.res$Sig <- as.numeric(m.rlm.res$Pval <= 0.05)
m.res.mat <- as.matrix(m.rlm.res)
m.res.var <- m.res.mat[c(3:4), ]
colnames(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
biv.univariate_list[["PromoterActivityBivalent"]] <- as.list(m.res.var[1, ])
biv.univariate_list[["PromoterActivityH3K4me3"]] <- as.list(m.res.var[1, ])

biv.rlm.df <- do.call(rbind.data.frame, biv.univariate_list)
biv.rlm.df$Predictor <- rownames(biv.rlm.df)
biv.rlm.df$Direction <- "NoEffect"
biv.rlm.df$Direction[biv.rlm.df$COEFF < 0 & biv.rlm.df$Sig == 1] <- "Less"
biv.rlm.df$Direction[biv.rlm.df$COEFF > 0 & biv.rlm.df$Sig == 1] <- "More"

# give the features more informative/better formated names
biv.rlm.df$Predictor[biv.rlm.df$Predictor == "Recip.Mean"] <- "Mean expression"
biv.rlm.df$Predictor[biv.rlm.df$Predictor == "cpg_RATIO"] <- "CpG dinucleotide ratio"
biv.rlm.df$Predictor[biv.rlm.df$Predictor == "CGI_SIZE.kb"] <- "CpG island size"
biv.rlm.df$Predictor[biv.rlm.df$Predictor == "PromoterActivityBivalent"] <- "Bivalent Promoter"
biv.rlm.df$Predictor[biv.rlm.df$Predictor == "PromoterActivityH3K4me3"] <- "Active Promoter"
biv.rlm.df$Model <- "Univariate"

write.table(biv.rlm.df,
            file="~/Dropbox/Noise_genomics/Model_results/mouse_ESC_bivalent-UnivariateLM.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)

######################
## Multivariate fit ##
######################

bivalent.fit <- rlm(Residual.CV2 ~  Recip.Mean + PromoterActivity + cpg_RATIO + CGI_SIZE.kb,
               data=mesc.merge)
bivalent.robust <- summary(bivalent.fit)
bivalent.rlm.res <- as.data.frame(bivalent.robust$coefficients)
n.params <- dim(bivalent.rlm.res)[1]
bivalent.rlm.res$Pval <- 2*pt(-abs(bivalent.rlm.res[, 3]), df=(dim(mesc.merge)[1]-1)-n.params)
bivalent.rlm.res$Sig <- as.numeric(bivalent.rlm.res$Pval <= 0.01)
bivalent.rlm.res$Predictor <- rownames(bivalent.rlm.res)
bivalent.rlm.res <- bivalent.rlm.res[-1, ]
colnames(bivalent.rlm.res) <- c("COEFF", "SE", "STAT", "P", "Sig", "Predictor")

bivalent.rlm.res$Direction <- "NoEffect"
bivalent.rlm.res$Direction[bivalent.rlm.res$COEFF < 0 & bivalent.rlm.res$Sig == 1] <- "Less"
bivalent.rlm.res$Direction[bivalent.rlm.res$COEFF > 0 & bivalent.rlm.res$Sig == 1] <- "More"

# give the features more informative/better formated names
bivalent.rlm.res$Predictor[bivalent.rlm.res$Predictor == "Recip.Mean"] <- "Mean expression"
bivalent.rlm.res$Predictor[bivalent.rlm.res$Predictor == "PromoterActivityBivalent"] <- "Bivalent Promoter"
bivalent.rlm.res$Predictor[bivalent.rlm.res$Predictor == "PromoterActivityH3K4me3"] <- "Active Promoter"
bivalent.rlm.res$Predictor[bivalent.rlm.res$Predictor == "cpg_RATIO"] <- "CpG dinucleotide ratio"
bivalent.rlm.res$Predictor[bivalent.rlm.res$Predictor == "CGI_SIZE.kb"] <- "CpG island size"

bivalent.rlm.res$Model <- "Multivariate"

biv.plot <- ggplot(bivalent.rlm.res,
                   aes(x=reorder(Predictor, -STAT),
                       y=STAT, fill=Direction)) +
  geom_point(alpha=0.55, shape=21, size=5) + theme_mike() +
  scale_fill_manual(values=effect.cols) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
  labs(x="Annotation", y="t-statistic") +
  guides(fill=FALSE) +
  scale_y_continuous(limits=c(-30, 30), oob=squish) +
  geom_hline(mapping=aes(yintercept=0), linetype="dashed", colour="grey")

ggsave(biv.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mESC_bivalent-multivariateLM.png",
       height=4.75, width=5.75, dpi=300)

rcv.biv <- ggplot(mesc.merge,
                  aes(x=PromoterActivity, y=Residual.CV2, colour=PromoterActivity)) +
  geom_jitter(alpha=0.5) + theme_mike() +
  geom_boxplot(outlier.colour='black', outlier.size=1) +
  scale_colour_manual(values=c("#a614ef", "#11d101", "#ff7306",  "#f0db87")) +
  labs(x="Promoter State", y=expression(paste("Residual CV"^2))) +
  guides(fill=FALSE, colour=FALSE) +
  scale_x_discrete(labels=c("Repressed Promoter", "Bivalent Promoter", "Active Promoter")) +
  theme(axis.title.x=element_blank())

ggsave(rcv.biv,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mESC_rCV_bivalency-boxplot.png",
       width=6.75, height=4.75, dpi=300)

mean.biv <- ggplot(mesc.merge,
                   aes(x=PromoterActivity, y=Mean, colour=PromoterActivity)) +
  geom_jitter(alpha=0.5) + theme_mike() +
  geom_boxplot(outlier.colour='black', outlier.size=1) +
  scale_colour_manual(values=c("#a614ef", "#11d101", "#ff7306",  "#f0db87")) +
  labs(x="Promoter State", y=expression(paste("mean log"[2], " Expression"))) +
  guides(fill=FALSE, colour=FALSE) +
  scale_x_discrete(labels=c("Repressed Promoter", "Bivalent Promoter", "Active Promoter")) +
  theme(axis.title.x=element_blank())

ggsave(mean.biv,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mESC_mean_bivalency-boxplot.png",
       width=6.75, height=4.75, dpi=300)

# plot uni and multivariate results together
biv.rlm.all <- do.call(rbind.data.frame, list("Univariate"=biv.rlm.df,
                                               "Multivariate"=bivalent.rlm.res))

biv.all.plot <- ggplot(biv.rlm.all,
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
  scale_y_continuous(limits=c(-30, 30), oob=squish)

ggsave(biv.all.plot,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/mouse_ESC_bivalency-allLM.png",
       height=4.75, width=5.75, dpi=300)
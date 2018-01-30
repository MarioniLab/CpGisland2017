## Investigating genomic factors that influence gene expression noise
# use the mean expression as the model weights, check the residuals for heteroscedasticity
# still use a robust linear model to deal with the heavy tails of the CV^2 though
# I want to bin the mean, into say 100 genes per bin and calculate the variance
# within each bin

source("~/Dropbox/R_sessions/Noise/mESC_chromHMM.R")
source("~/Dropbox/R_sessions/Noise/genomic_noise_features.R")
source("~/Dropbox/R_sessions/Noise/mESC_genomic_noise_features-spline_fit.R")

library(ggplot2)
library(robustbase)
library(scales)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

# calculate 1/mean
mesc.gene.summary$Reciprocal.Mean <- 1/mesc.gene.summary$Mean

mesc.vars <- colnames(genomic.features)[2:16]
mesc.vars <- mesc.vars[!grepl(mesc.vars, pattern="(NMI)|(CGI_SIZE)|(cpg_)|(PHAST)|(Mean)")]

mesc.genomic.vars <- paste(c(mesc.vars),
                           collapse=" + ")

mesc.match <- merge(mesc.gene.summary, genomic.features,
                    by='GENE')
mesc.match$CGI_SIZE.kb <- mesc.match$CGI_SIZE/1000

#########################################
## univariate robust linear regression ##
#########################################
# this is plotting for presentation/manuscript figures
mesc.univariate_list <- list()

mesc.var.names <- unlist(strsplit(mesc.genomic.vars, split=" + ", fixed=T))
# remove redundant variables, i.e CpG island AND NMI
# remove CpG island characteristics
mesc.var.names <- mesc.var.names[!grepl(mesc.var.names, pattern="(NMI)|(CGI_SIZE)|(cpg_)|(PHAST)|(Mean)")]

# set a higher number of iterations, as sometimes doesn't converge
model.control <- lmrob.control(max.it=500, k.max=500, rel.tol=1e-6)

# check this model for the relationship between mean expression and residuals
null.glm.form <- as.formula(paste("CV2", "Mean", sep=" ~ "))
null.model <- lmrob(null.glm.form, data=mesc.match, control=model.control)

# plot the residuals vs fitted values to get a handle on the heteroscedasticity
plot(null.model$fitted.values, null.model$residuals)

# now try fitting with the reciprocal of the mean as the weights
null.model <- lmrob(null.glm.form, data=mesc.match, control=model.control,
                    weights=mesc.match$Reciprocal.Mean)

# plot the residuals vs fitted values to get a handle on the heteroscedasticity
plot(null.model$fitted.values, null.model$residuals)

## bin the mean into 100 equal sized bins, calculate the variance for each bin
mesc.match$Mean.bin <- cut(mesc.match$Mean, breaks=quantile(x = mesc.match$Mean, c(1:50)/50))
bin.var <- by(mesc.match$Mean, mesc.match$Mean.bin, FUN=function(Q) var(Q, na.rm=TRUE))

weight.df <- do.call(cbind.data.frame,
                     list("Mean.bin"=levels(mesc.match$Mean.bin),
                          "Bin.var"=as.numeric(bin.var)))

mesc.merge <- merge(mesc.match, weight.df, by='Mean.bin')

# use the inverse of the bin variance as the model weight
# now try fitting with the reciprocal of the mean as the weights
# set a higher number of iterations, as sometimes doesn't converge
model.control <- lmrob.control(max.it=500, k.max=500, rel.tol=1e-6)

null.model <- lmrob(null.glm.form, data=mesc.merge, control=model.control,
                    weights=mesc.merge$Bin.var)

# plot the residuals vs fitted values to get a handle on the heteroscedasticity
plot(null.model$fitted.values, null.model$residuals)

######################################
# Fit a model to each expression bin #
######################################
# check whether the effects of CpG island and size hold across the full expression spectrum
bins <- levels(mesc.merge$Mean.bin)
cpg.univariate_list <- list()

## CpG islands
for(q in seq_along(bins)){
  m.bin <- bins[q]
  bin.data <- mesc.merge[mesc.merge$Mean.bin == m.bin, ]
  #.variable <- paste("CGI_SIZE.kb", sep=" + ")
  .variable <- "N_CpG"
  .glm.form <- as.formula(paste("CV2", .variable, sep=" ~ "))
  
  m.rlm <- lmrob(.glm.form, data=bin.data, control=model.control) 
  
  m.robust <- summary(m.rlm)
  m.rlm.res <- as.data.frame(m.robust$coefficients)
  m.rlm.res$Sig <- as.numeric(m.rlm.res$`Pr(>|t|)` <= 0.05)
  m.res.mat <- as.matrix(m.rlm.res)
  m.res.var <- m.res.mat[2, ]
  names(m.res.var) <- c("COEFF", "SE", "STAT", "P", "Sig")
  cpg.univariate_list[[m.bin]] <- as.list(m.res.var)
}

cpg.rlm.df <- do.call(rbind.data.frame, cpg.univariate_list)
cpg.rlm.df$Predictor <- rownames(cpg.rlm.df)
cpg.rlm.df$Direction <- "NoEffect"
cpg.rlm.df$Direction[cpg.rlm.df$COEFF < 0 & cpg.rlm.df$Sig == 1] <- "Less"
cpg.rlm.df$Direction[cpg.rlm.df$COEFF > 0 & cpg.rlm.df$Sig == 1] <- "More"

ggplot(cpg.rlm.df,
       aes(x=Predictor,
           y=STAT, fill=Direction)) +
  geom_point(alpha=0.55, shape=21, size=5) + theme_mike() +
  scale_fill_manual(values=c("#62148f", "#feaf10", "#878787")) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))

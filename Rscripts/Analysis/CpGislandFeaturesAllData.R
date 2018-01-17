library(ggplot2)
library(scales)
source("~/Dropbox/R_sessions/GGMike/theme_mike.R")

# read in mouse data
mouse.esc <- read.table("~/Dropbox/Noise_genomics/Model_results/mouse_ESC_CGI-multivariateRLM.tsv",
                        h=TRUE, sep="\t", stringsAsFactors=FALSE)

mouse.tcell <- read.table("~/Dropbox/Noise_genomics/Model_results/mouse_Tcell_CGI-multivariateRLM.tsv",
                          h=TRUE, sep="\t", stringsAsFactors=FALSE)

# mouse.liver <- read.table("~/Dropbox/Noise_genomics/Model_results/mouse_liver_CGI-multivariateRLM.tsv",
#                           h=TRUE, sep="\t", stringsAsFactors=FALSE)

# read in human data
human.alpha <- read.table("~/Dropbox/Noise_genomics/Model_results/human_AlphaIslet_CGI-multivariateRLM.tsv",
                          h=TRUE, sep="\t", stringsAsFactors=FALSE)

human.beta <- read.table("~/Dropbox/Noise_genomics/Model_results/human_BetaIslet_CGI-multivariateRLM.tsv",
                         h=TRUE, sep="\t", stringsAsFactors=FALSE)

human.esc <- read.table("~/Dropbox/Noise_genomics/Model_results/human_ESC_CGI-multivariateRLM.tsv",
                        h=TRUE, sep="\t", stringsAsFactors=FALSE)

# check the same column names across tables
all(colnames(human.esc) == colnames(human.beta))
all(colnames(mouse.tcell) == colnames(mouse.esc))
all(colnames(human.esc)  == colnames(mouse.esc))

all.merge <- do.call(rbind.data.frame, list("mESC"=mouse.esc,
                                            "mTcell"=mouse.tcell,
                                            "hAlpha"=human.alpha,
                                            "hBeta"=human.beta,
                                            "hESC"=human.esc))
all.merge$Sig <- as.factor(all.merge$Sig)

dir.cols <- c("#62148f", "#feaf10", "#878787")
species.cols <- c("#0008FF", "#D90000")
names(species.cols) <- c("Human", "Mouse")

sig.alpha <- c(0.2, 1)
names(sig.alpha) <- levels(all.merge$Sig)

all.merge$Predictor <- reorder(all.merge$Predictor,
                               -all.merge$STAT)
all.merge$Id <- as.factor(with(all.merge, 
                               order(-ave(all.merge$STAT,
                                          all.merge$Predictor, FUN=max), all.merge$STAT)))

# order the predictor variables by STAT in reverse ordrer
# need to manually add the panel names as the X-axis ticks
all.lm <- ggplot(all.merge,
                 aes(x=reorder(Id, -STAT),
                     y=STAT, fill=Species, 
                     shape=Tissue,
                     alpha=Sig)) +
  geom_hline(mapping=aes(yintercept=0), linetype="dashed", colour="grey") +
  geom_point(size=6) + 
  theme_mike() +
  scale_y_continuous(limits=c(-30, 30), oob=squish) +
  scale_fill_manual(values=species.cols) +
  scale_shape_manual(values=c(21:25)) +
  scale_alpha_manual(values=sig.alpha) +
  labs(x="Annotation", y="t-statistic") +
  guides(fill=FALSE, shape=FALSE, alpha=FALSE) +
  facet_grid(~Predictor, shrink=TRUE,
             space="free_x",
             scales="free_x", switch="x") +
  theme(panel.spacing=unit(0.25, "lines"),
        strip.text=element_text(angle=90, vjust=1, hjust=0.5, size=16,
                                family='Helvetica', face='plain'),
        strip.background=element_blank()) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank())

ggsave(all.lm,
       filename="~/Dropbox/Noise_genomics/Figures/ms_figures/AllLM-CGI_figure.png",
       height=5.25, width=11.75, dpi=300)


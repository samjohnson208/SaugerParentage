## Sequoia_params_plot.R by SPJ 091825
## Purpose: I have compiled a color coded table called sequoia_params that tells
# me a bunch of info about the sequoia output for a mess of genotype matrices with
# different md, sequoia error, thinning, and number of snps. I am now interested in
# visualizing how the assignment, accuracy, and composite scores (assignment x accuracy)
# vary with respect to these dataset parameters. I made sequoia_params_toplot.csv
# that will serve as the input for this script. We're looking for a six panel plot,
# with Assign, Accuracy, and Composite all plotted against thin or n_snps, for md10
# and md5. 

## SEE NOTES AT THE END FOR STATUS

## libraries
################################################################################

# install.packages("sequoia")
# 65
library(tidyverse)

# install.packages("dplyr")
# 65
library(dplyr)

# install.packages("ggplot2")
# 65
library(ggplot2)

install.packages("ggpubr")
#65
library(ggpubr)

install.packages("gridExtra")
#65
library(gridExtra)

install.packages("grid")
#65
library(grid)

################################################################################

setwd("/Users/samjohnson/Desktop/")

################################################################################

dat <- read.csv(file = "sequoia_params_toplot.csv", header = TRUE)

################################################################################

head(dat)
tail(dat)
miss10 <- dat[1:51,]
miss5 <- dat[52:nrow(dat),]

p1 <- ggplot(miss10, aes(x = n_snps, y = assignment_rate, color = as.factor(sequoia_err))) +
          geom_point() +
          geom_line() +
          # geom_smooth(method = "lm", se = FALSE) +
          labs(x = "n SNPs", y = "Assignment Rate (%)", color = "Sequoia Error Rate") +
          theme_bw() +
          theme(axis.title = element_text(size = 14)) +
          xlim(100, 1225) + ylim(20, 95)

p2 <- ggplot(miss10, aes(x = n_snps, y = accuracy_rate, color = as.factor(sequoia_err))) +
  geom_point() +
  geom_line() +
  # geom_smooth(method = "lm", se = FALSE) +
  labs(x = "n SNPs", y = "Accuracy Rate (%)", color = "sequoia error rate") +
  theme_bw() +
  theme(axis.title = element_text(size = 14)) +
  xlim(100, 1225) + ylim(20, 95)

p3 <- ggplot(miss10, aes(x = n_snps, y = composite_score, color = as.factor(sequoia_err))) +
  geom_point() +
  geom_line() +
  # geom_smooth(method = "lm", se = FALSE) +
  labs(x = "n SNPs", y = "Composite Score", color = "sequoia error rate") +
  theme_bw() +
  theme(axis.title = element_text(size = 14)) +
  xlim(100, 1225) + ylim(20, 95)

p4 <- ggplot(miss5, aes(x = n_snps, y = assignment_rate, color = as.factor(sequoia_err))) +
  geom_point() +
  geom_line() +
  # geom_smooth(method = "lm", se = FALSE) + 
  labs(x = "n SNPs", y = "Assignment Rate (%)", color = "Sequoia Error Rate") +
  theme_bw() +
  theme(axis.title = element_text(size = 14)) +
  xlim(100, 1225) + ylim(20, 95)

p5 <- ggplot(miss5, aes(x = n_snps, y = accuracy_rate, color = as.factor(sequoia_err))) +
  geom_point() +
  geom_line() +
  # geom_smooth(method = "lm", se = FALSE) +
  labs(x = "n SNPs", y = "Accuracy Rate (%)", color = "sequoia error rate") +
  theme_bw() +
  theme(axis.title = element_text(size = 14)) +
  xlim(100, 1225) + ylim(20, 95)

p6 <- ggplot(miss5, aes(x = n_snps, y = composite_score, color = as.factor(sequoia_err))) +
  geom_point() +
  geom_line() +
  # geom_smooth(method = "lm", se = FALSE) +
  labs(x = "n SNPs", y = "Composite Score", color = "sequoia error rate") +
  theme_bw() +
  theme(axis.title = element_text(size = 14)) +
  xlim(100, 1225) + ylim(20, 95)

# ggarrange(p1, p2, p3, p4, p5, p6, 
#           nrow = 2, ncol = 3,
#           common.legend = TRUE,
#           legend = "right")

row1 <- ggarrange(p1, 
                  p2 + theme(legend.position = "none"), 
                  p3 + theme(legend.position = "none"),
                  ncol = 3, common.legend = TRUE, legend = "right")
row2 <- ggarrange(p4, 
                  p5 + theme(legend.position = "none"), 
                  p6 + theme(legend.position = "none"),
                  ncol = 3, common.legend = TRUE, legend = "right")

row1_labeled <- annotate_figure(row1, left = text_grob("10% Missing Data/Site", rot = 90, size = 16, hjust = 0.45))
row2_labeled <- annotate_figure(row2, left = text_grob("5% Missing Data/Site",  rot = 90, size = 16, hjust = 0.45))

final <- ggarrange(row1_labeled, row2_labeled, ncol = 1, heights = c(1,1), common.legend = TRUE, legend = "right")
final

col_labels <- ggarrange(
  text_grob("Assignment", size = 16, hjust = 0.75),
  text_grob("Accuracy", size = 16, hjust = 0.5),
  text_grob("Composite", size = 16, hjust = 0.5),
  ncol = 3
)

final <- ggarrange(
  col_labels,
  ggarrange(row1_labeled, row2_labeled, ncol = 1, heights = c(1,1)),
  ncol = 1,
  heights = c(0.05, 0.5))

final

# okay, got this how it is right now, put it into powerpoint, and changed the placement
# of the column labels, as well as the order of the legend. see sequoia_params_plots.pptx
# for most up to date version.




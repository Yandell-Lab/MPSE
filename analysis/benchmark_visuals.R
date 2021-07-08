#!/usr/bin/env R

library(tidyverse)

setwd("/home/bennet/Documents/PhD/Yandell/software/MPSE")
dr <- read_csv("analysis/bernoulli_nb_diagnostic_rate.csv")

line_plot <- ggplot(dr, aes(x=list_fraction, y=diag_rate)) +
	geom_line(color="red", size=1) + 
  scale_y_continuous(name="Diagnostic rate", 
                     limits=c(0.0, 1.0),
                     breaks=c(0.00, 0.18, 0.35, 0.50, 1.00),
                     labels=c("0%","18%","35%","50%","100%"),
                     minor_breaks=NULL,
                     expand=c(0,0)) + 
  scale_x_continuous(name="Top scoring fraction of probands",
                     breaks=c(0.00, 0.21, 0.50, 1.00),
                     labels=c("0%","21%","50%","100%"),
                     minor_breaks=NULL,
                     expand=c(0,0)) +
	geom_hline(yintercept=c(0.18, 0.35), linetype=3) + 
	geom_vline(xintercept=0.21, alpha=0.3, size=1) + 
  theme_bw()




rad <- read_csv("analysis/bernoulli_nb_predictions.csv")
edw <- read_csv("analysis/edw_case_predictions.csv") %>% 
  mutate(y=0)
neo <- read_csv("analysis/neo_case_predictions.csv") %>% 
  mutate(y=0)

pos <- rad %>% filter(neg_proba!=1) %>% mutate(scr=-log(neg_proba))
neg <- rad %>% filter(neg_proba==1) %>% mutate(scr=-log(1/pos_proba))
scrs <- bind_rows(pos, neg)

dens_plot <- ggplot(scrs, aes(x=scr, colour=factor(outcome))) + 
  geom_density(adjust=2) + 
    geom_jitter(data=edw, aes(-log(neg_proba), y, shape=factor(diagnostic)), 
		              height = 0.0005, size=2, inherit.aes = FALSE) + 
  scale_shape_discrete(guide = guide_legend(title=NULL), solid = FALSE,
		                              labels = c("Not Diagnostic", "Diagnostic")) +
  scale_x_continuous(name="MPSE Score", limits=c(-100,520),
		                          breaks=seq(-100,500,50)) + 
  scale_colour_discrete(guide = guide_legend(title=NULL),
			                        labels = c("Not Sequenced", "Sequenced"))


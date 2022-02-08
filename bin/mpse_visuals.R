#!/usr/bin/env R

library(readr)
library(dplyr)
library(ggplot2)


args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
	stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
	args[2] = "out.txt"
}


setwd("~/Documents/PhD/Yandell/software/MPSE")
train <- read_delim("analysis/example/tables/training_predictions.csv",
                    delim = "\t", 
                    col_types = cols(
                      sex = col_factor(),
                      race = col_factor(),
                      ethnicity = col_factor(),
                      seq_status = col_factor(),
                      diagnostic = col_factor(),
                      incidental = col_factor(),
                      classes = col_factor()
                    ))

train_sample <- read_delim("analysis/example/tables/training_predictions_sample.csv",
                           delim = "\t",
                           col_types = cols(
                             sex = col_factor(),
                             race = col_factor(),
                             ethnicity = col_factor(),
                             seq_status = col_factor(),
                             #diagnostic = col_factor(),
                             incidental = col_factor(),
                             classes = col_factor()
                           )) %>% 
  mutate(neg_log_proba = -log(neg_proba),
         pos_log_proba = -log(pos_proba)) %>% 
  arrange(desc(neg_log_proba), desc(pos_proba)) %>% 
  mutate(rank_cumsum = cumsum(as.integer(diagnostic)),
         diag_rate = rank_cumsum / 1:nrow(.),
         list_fraction = 1:nrow(.) / nrow(.))


line_plot <- ggplot(train_sample, aes(x=list_fraction, y=diag_rate)) +
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


pos <- train %>% filter(neg_proba!=1) %>% mutate(scr=-log(neg_proba))
neg <- train %>% filter(neg_proba==1) %>% mutate(scr=-log(1-pos_proba))

scrs <- bind_rows(pos, neg)


dens_plot <- ggplot(scrs, aes(x=scr, colour=seq_status)) + 
  geom_density(adjust=2) +
  scale_x_continuous(name="MPSE Score", limits=c(-100,520),
                     breaks=seq(-100,500,50)) +
  scale_colour_discrete(guide = guide_legend(title=NULL),
                        labels = c("Not Sequenced", "Sequenced"))


ggplot(scrs, aes(x=scr, y=seq_status)) + 
  geom_violin() + 
  geom_jitter(aes(colour=classes), alpha=0.3)


#!/usr/bin/env R

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)



setwd("~/Documents/PhD/Yandell/software/MPSE")
train <- read_delim("analysis/javier_power/tables/training_predictions.tsv",
                    delim = "\t", 
                    col_types = cols(
                      sex = col_factor(),
                      race = col_factor(),
                      ethnicity = col_factor(),
                      seq_status = col_factor(),
                      diagnostic = col_factor(),
                      incidental = col_factor(),
                      class = col_factor()
                    ))


train_sample <- read_delim("analysis/javier_power/tables/training_predictions_sample.tsv",
                           delim = "\t",
                           col_types = cols(
                             sex = col_factor(),
                             race = col_factor(),
                             ethnicity = col_factor(),
                             seq_status = col_factor(),
                             #diagnostic = col_factor(),
                             incidental = col_factor(),
                             class = col_factor()
                           )) %>% 
  mutate(diagnostic = replace_na(diagnostic, 0)) %>% 
  arrange(desc(pos_log_proba), neg_log_proba) %>% 
  mutate(rank_cumsum = cumsum(as.integer(diagnostic)),
         diag_rate = rank_cumsum / 1:nrow(.),
         list_fraction = 1:nrow(.) / nrow(.))


valid <- read_delim("analysis/feature_select_scores.tsv",
                    delim = "\t", 
                    col_types = cols(
                      seq_status = col_factor(),
                      #diagnostic = col_factor(),
                      class = col_factor()
                    )) %>% 
  arrange(desc(pos_log_proba), neg_log_proba) %>% 
  mutate(rank_cumsum = cumsum(as.integer(diagnostic)),
         diag_rate = rank_cumsum / 1:nrow(.),
         list_fraction = 1:nrow(.) / nrow(.))


train_valid_join <- bind_rows("train"=select(train_sample, scr, pos_log_proba, neg_log_proba),
                              "valid"=select(valid, scr, pos_log_proba, neg_log_proba, diagnostic), 
                              .id = "cohort") %>% 
  arrange(desc(pos_log_proba), neg_log_proba) %>% 
  mutate(list_fraction = 1:nrow(.) / nrow(.))


line_plot <- ggplot(train_sample, aes(x=list_fraction, y=diag_rate)) +
  geom_line(color="red", size=1) + 
  scale_y_continuous(name="Diagnostic rate", 
                     limits=c(-0.05, 1.0),
                     breaks=c(0.00, 0.18, 0.35, 0.50, 1.00),
                     labels=c("0%","18%","35%","50%","100%"),
                     minor_breaks=NULL,
                     expand=c(0,0)) + 
  scale_x_continuous(name="Top scoring fraction of probands",
                     breaks=c(0.00, 0.20, 0.50, 1.00),
                     labels=c("0%","20%","50%","100%"),
                     minor_breaks=NULL,
                     expand=c(0,0)) + 
  geom_jitter(data=filter(train_valid_join, cohort=="valid"), 
              aes(x=list_fraction, y=0, shape=factor(diag_binary)),
              color="black", alpha=0.7, height=0.02, width=0.005) + 
  scale_shape_manual(guide = guide_legend(title=NULL),
                     values=c(1,17),
                     labels = c("Not Diagnostic", "Diagnostic")) +
  geom_hline(yintercept=c(0.18, 0.35), linetype=3) + 
  geom_vline(xintercept=0.20, alpha=0.3, size=1) + 
  theme_bw()


double_line_plot <- ggplot(train_sample, aes(x=list_fraction, y=diag_rate)) +
  geom_line(color="red", size=1) + 
  scale_y_continuous(name="Diagnostic rate", 
                     limits=c(-0.05, 1.0),
                     breaks=c(0.00, 0.18, 0.35, 0.50, 1.00),
                     labels=c("0%","18%","35%","50%","100%"),
                     minor_breaks=NULL,
                     expand=c(0,0)) + 
  scale_x_continuous(name="Top scoring fraction of probands",
                     breaks=c(0.00, 0.20, 0.50, 1.00),
                     labels=c("0%","20%","50%","100%"),
                     minor_breaks=NULL,
                     expand=c(0,0)) + 
  geom_line(data=valid, aes(x=list_fraction, y=diag_rate),
            color="blue", size=1) +
  geom_hline(yintercept=c(0.18, 0.35), linetype=3) + 
  geom_vline(xintercept=0.20, alpha=0.3, size=1) + 
  theme_bw()


dens_plot <- ggplot(train, aes(x=pos_log_proba, color=seq_status)) + 
  geom_density(adjust=1) + 
  scale_y_continuous(name="Density", 
                     limits=c(-0.005, 0.06),
                     minor_breaks=NULL) + 
  scale_x_continuous(name="MPSE Score",
                     limits=c(-70, 15),
                     breaks=c(-60, -40, -20, 0),
                     minor_breaks=NULL) +
  scale_color_discrete(guide = guide_legend(title=NULL),
                       labels = c("Not Sequenced", "Sequenced")) +
  geom_jitter(data=valid, aes(x=pos_log_proba, y=0, shape=factor(diag_binary)), 
              color="black", alpha=0.7, height=0.005, width=0.6) + 
  scale_shape_manual(guide = guide_legend(title=NULL),
                     values = c(1,17),
                     labels = c("Not Diagnostic", "Diagnostic")) + 
  theme_bw()


dens_plot2 <- ggplot(train, aes(x=scr, color=rev(seq_status))) + 
  geom_density(adjust=1.5) + 
  scale_y_continuous(name="Density", 
                     limits=c(-0.004, 0.032),
                     minor_breaks=NULL) + 
  scale_x_continuous(name="MPSE Score",
                     limits=c(-80, 120),
                     breaks=seq(-60, 120, 30),
                     minor_breaks=NULL) +
  scale_color_discrete(guide = guide_legend(title=NULL),
                       labels = c("Sequenced", "Not Sequenced")) +
  geom_jitter(data=valid, aes(x=scr, y=-0.002, shape=factor(diagnostic)), 
              color="black", alpha=0.7, height=0.002, width=0.001) + 
  scale_shape_manual(guide = guide_legend(title=NULL),
                     values = c(1,17),
                     labels = c("Not Diagnostic", "Diagnostic")) + 
  theme_bw()


t.test(scr ~ diagnostic, data=valid)
t.test(scr ~ diagnostic, data=valid, subset=valid$diag_multilevel != "Pending")
t.test(scr ~ diagnostic, data=valid, subset=valid$diag_multilevel %in% c("Yes","No"))







dens_plot <- ggplot(scrs, aes(x=scr, colour=seq_status)) + 
  geom_density(adjust=2) +
  scale_x_continuous(name="MPSE Score", limits=c(-100,520),
                     breaks=seq(-100,500,50)) +
  scale_colour_discrete(guide = guide_legend(title=NULL),
                        labels = c("Not Sequenced", "Sequenced"))


ggplot(scrs, aes(x=scr, y=seq_status)) + 
  geom_violin() + 
  geom_jitter(aes(colour=class), alpha=0.3)


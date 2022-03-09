#!/usr/bin/env R

library(readr)
library(tibble)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)



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


valid <- read_delim("analysis/javier_power/tables/validation_predictions.tsv",
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


ggplot(train_sample, aes(x=list_fraction, y=diag_rate)) +
  geom_smooth(method = lm, formula = y ~ poly(x, 25), se = FALSE, color="red", size=1) + 
  scale_y_continuous(name="Diagnostic rate", 
                     limits=c(0, 1.0),
                     breaks=c(0.00, 0.18, 0.35, 0.50, 1.00),
                     labels=c("0%","18%","35%","50%","100%"),
                     minor_breaks=NULL,
                     expand=c(0,0)) + 
  scale_x_continuous(name="Top scoring fraction of probands",
                     limits=c(0, 1.0),
                     breaks=c(0.00, 0.21, 0.50, 1.00),
                     labels=c("0%","21%","50%","100%"),
                     minor_breaks=NULL,
                     expand=c(0,0)) + 
  geom_line(data=valid, aes(x=list_fraction, y=diag_rate),
            color="blue", size=1) +
  geom_hline(yintercept=0.18, linetype="dashed") + 
  geom_hline(yintercept=0.35, linetype="dotted") + 
  geom_vline(xintercept=0.21, alpha=0.3, size=0.7) + 
  annotate(geom="text", label="Top 100 scoring probands", x=0.36, y=0.9) +
  annotate(geom="text", label="41% Diagnostic Rate", x=0.36, y=0.85) +
  annotate(geom="text", label="RCIGM Diag. Rate (35%)", x=0.85, y=0.4) +
  annotate(geom="text", label="NSIGHT2 rate of Mendelian", x=0.4, y=0.2) +
  annotate(geom="text", label="disease in NICU (18%)", x=0.4, y=0.16) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())


dens_plot <- ggplot(train, aes(x=scr, fill=seq_status)) + 
  geom_density(alpha=0.6) + 
  geom_segment(aes(x=-75, xend=350, y=0, yend=0), inherit.aes=FALSE, size=0.3, alpha=0.5) +
  scale_y_continuous(name="Density", 
                     limits=c(-0.004, 0.04),
                     minor_breaks=NULL) + 
  scale_x_continuous(name="MPSE Score",
                     limits=c(-75, 350),
                     breaks=seq(-60, 350, 30),
                     minor_breaks=NULL) +
  scale_fill_manual(values = c("#00BFC4","#F8766D"),
                    guide = guide_legend(title=NULL),
                    labels = c("Not Sequenced", "Sequenced")) +
  geom_jitter(data=valid, aes(x=scr, y=-0.002, shape=factor(diagnostic)), 
              inherit.aes=FALSE, color="black", alpha=0.7, height=0.002, width=0.001) + 
  scale_shape_manual(guide = guide_legend(title=NULL),
                     values = c(1,17),
                     labels = c("Not Diagnostic", "Diagnostic")) + 
  theme_bw()




binom_alphas <- c(0.001, 0.05, 0.25, 0.5, 0.75, 0.9, 1.0)
mean_scrs <- c(6.96, 9.77, 12.27, 13.58, 14.61, 14.68, 14.75)
alpha_tbl <- tibble(alpha=binom_alphas, scr=mean_scrs)
alpha_gg <- ggplot(data=alpha_tbl, aes(x=binom_alphas, y=mean_scrs)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(name="Binomial Test alpha",
                     breaks=binom_alphas,
                     labels=c("0.001","0.05","0.25","0.5","0.75","0.9","1.0"),
                     minor_breaks=NULL) +
  scale_y_continuous(name="Mean MPSE score (NeoSeq; n=36)",
                     limits=c(0,15),
                     breaks=seq(0, 15, 2.5),
                     labels=c("0","2.5","5","7.5","10","12.5","15"),
                     minor_breaks=NULL) + 
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


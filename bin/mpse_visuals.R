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
  geom_smooth(method = lm, formula = y ~ poly(x, 25), se = FALSE, color="red", size=1) + 
  scale_y_continuous(name="Diagnostic rate", 
                     limits=c(0, 1.0),
                     breaks=c(0.00, 0.18, 0.35, 0.50, 1.00),
                     labels=c("0%","18%","35%","50%","100%"),
                     minor_breaks=NULL,
                     expand=c(0,0)) + 
  scale_x_continuous(name="Top scoring fraction of probands",
                     limits=c(0, 1.0),
                     breaks=c(0.00, 0.20, 0.50, 1.00),
                     labels=c("0%","20%","50%","100%"),
                     minor_breaks=NULL,
                     expand=c(0,0)) + 
  geom_line(data=valid, aes(x=list_fraction, y=diag_rate),
            color="blue", size=1) +
  geom_hline(yintercept=0.18, linetype="dashed") + 
  geom_hline(yintercept=0.35, linetype="dotted") + 
  geom_vline(xintercept=0.20, alpha=0.3, size=0.7) + 
  geom_segment(aes(x=0.7, xend=0.6, y=0.47, yend=0.37), 
               inherit.aes=FALSE, size=0.3, 
               arrow=arrow(length=unit(0.01, "npc"))) +
  annotate(geom="text", label="Top 20% scoring probands =", x=0.36, y=0.9) +
  annotate(geom="text", label="41% Diagnostic Rate", x=0.36, y=0.85) +
  annotate(geom="text", label="RCIGM Diag. Rate (35%)", x=0.8, y=0.5) +
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



train_alpha <- read_csv("analysis/alpha_training_preds.tsv",
                        col_types = cols(
                          alpha = col_factor(levels=c("1e-06","1e-05","1e-04","1e-03",
                                                 "0.05","0.25","0.5","0.75","0.9","1.0")),
                          seq_status = col_factor(),
                          diagnostic = col_factor(),
                          class = col_factor(),
                          scr = col_double()
                        ))
alpha_summary <- train_alpha %>% 
  group_by(alpha, seq_status) %>% 
  summarise("mean_scrs"=mean(scr),
            "median_scrs"=median(scr))
alpha_summary_long <- alpha_summary %>% 
  pivot_longer(cols = c("mean_scrs","median_scrs"))

alphas <- c("1e-06","1e-05","1e-04","1e-03","0.05","0.25","0.5","0.75","0.9","1.0")
mean_scrs <- c(3.01, 3.71, 5.13, 6.96, 9.77, 12.27, 13.58, 14.61, 14.68, 14.75)
median_scrs <- c(0.05, 1.51, -0.11, 2.11, 5.10, 7.56, 8.32, 8.82, 8.91, 8.74)
neo_alpha <- tibble(alpha=factor(alphas), mean_scrs=mean_scrs, median_scrs=median_scrs)

ggplot(data=alpha_summary, aes(x=alpha, y=median_scrs, group=seq_status)) +
  geom_line(aes(color=seq_status)) + 
  geom_line(data=neo_alpha, aes(x=alpha, y=median_scrs, group=1), inherit.aes = FALSE) +
  scale_x_discrete(name="Binomial Test alpha") +
  scale_y_continuous(name="Mean MPSE score (NeoSeq; n=36)",
                     limits=c(-45,15),
                     breaks=seq(-45, 15, 5),
                     minor_breaks=NULL) +
  theme_bw()

ggplot(data=alpha_summary, aes(x=alpha, y=mean_scrs, group=seq_status)) +
  geom_line(aes(color=seq_status)) + 
  geom_line(data=neo_alpha, aes(x=alpha, y=mean_scrs, group=1), inherit.aes = FALSE) +
  scale_x_discrete(name="Binomial Test alpha") +
  scale_y_continuous(name="Mean MPSE score (NeoSeq; n=36)",
                     limits=c(-35,30),
                     breaks=seq(-35, 30, 5),
                     minor_breaks=NULL) +
  theme_bw()

ggplot(data=train_alpha, aes(x=alpha, y=scr, fill=seq_status)) + 
  geom_boxplot(outlier.shape = NA) + 
  scale_x_discrete(name="Binomial Test alpha") +
  scale_y_continuous(name="Mean MPSE score (NeoSeq; n=36)",
                     limits=c(-75,125),
                     breaks=seq(-75, 125, 50),
                     minor_breaks=NULL)

ggplot(data=alpha_summary_long, aes(x=alpha, y=value, group=seq_status)) +
  geom_line(aes(color=seq_status)) + 
  geom_line(data=neo_alpha, aes(x=alpha, y=mean_scrs, group=1), inherit.aes = FALSE) +
  facet_wrap(~ name) +
  scale_x_discrete(name="Binomial Test alpha") +
  scale_y_continuous(name="Mean MPSE score (NeoSeq; n=36)",
                     limits=c(-45,30),
                     breaks=seq(-45, 30, 5),
                     minor_breaks=NULL) +
  theme_bw()

# alphas <- c(0.001, 0.05, 0.25, 0.5, 0.75, 0.9, 1.0)
# mean_scrs <- c(6.96, 9.77, 12.27, 13.58, 14.61, 14.68, 14.75)
# alpha_tbl <- tibble(alpha=binom_alphas, scr=mean_scrs)
# alpha_gg <- ggplot(data=alpha_tbl, aes(x=alpha, y=scr)) +
#   geom_point() +
#   geom_line() +
#   scale_x_continuous(name="Binomial Test alpha",
#                      breaks=binom_alphas,
#                      labels=c("0.001","0.05","0.25","0.5","0.75","0.9","1.0"),
#                      minor_breaks=NULL) +
#   scale_y_continuous(name="Mean MPSE score (NeoSeq; n=36)",
#                      limits=c(0,15),
#                      breaks=seq(0, 15, 2.5),
#                      labels=c("0","2.5","5","7.5","10","12.5","15"),
#                      minor_breaks=NULL) + 
#   theme_bw()



samp_frac <- c(0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
mean_scrs1 <- c(9.00, 11.20, 8.47, 12.94, 13.3, 14.75)
mean_scrs2 <- c(6.22, 8.12, 10.75, 11.68, 15.60, 14.75)
mean_scrs3 <- c(7.62, 5.45, 10.31, 10.43, 14.33, 14.75)
mean_scrs4 <- c(7.36, 8.62, 10.98, 10.25, 12.68, 14.75)
mean_scrs5 <- c(5.75, 8.97, 10.09, 10.06, 14.26, 14.75)
samp_feat_tbl <- tibble(samp_frac=rep(samp_frac, 5), 
                        scr=c(mean_scrs1, mean_scrs2, mean_scrs3, mean_scrs4, mean_scrs5))
samp_feat_gg <- ggplot(data=samp_feat_tbl, aes(x=samp_frac, y=scr)) + 
  geom_point(shape=1, color="red") +
  stat_summary(fun.data=mean_sdl, fun.args=list(mult=1), geom="errorbar", width=0.01) +
  stat_summary(fun=mean, geom="line", aes(group=1)) +
  scale_x_continuous(name="% Features Sampled",
                     breaks=samp_frac,
                     labels=c("50%","60%","70%","80%","90%","100%"),
                     minor_breaks=NULL) +
  scale_y_continuous(name="Mean MPSE score (NeoSeq; n=36)",
                     limits=c(0,16),
                     breaks=seq(0, 15, 2.5),
                     labels=c("0","2.5","5","7.5","10","12.5","15"),
                     minor_breaks=NULL) + 
  theme_bw()



fudge_n <- c(-30, -20, -10, -5, 0, 5, 10, 20, 30)
mean_scrs1 <- c()
mean_scrs2 <- c()
mean_scrs3 <- c()
mean_scrs4 <- c()
mean_scrs5 <- c()
fudge_terms_tbl <- tibble(fudge_n=rep(fudge_n, 5),
                          scr=c(mean_scrs1, mean_scrs2, mean_scrs3, mean_scrs4, mean_scrs5))
fudge_terms_gg <- ggplot(data=fudge_terms_tbl, aes(x=fudge_n, y=scr)) + 
  geom_point(shape=1, color="red") +
  stat_summary(fun.data=mean_sdl, fun.args=list(mult=1), geom="errorbar", width=0.01) +
  stat_summary(fun=mean, geom="line", aes(group=1)) +
  scale_x_continuous(name="# Terms Added/Removed",
                     breaks=fudge_n,
                     labels=c("30","20","10","5","0","5","10","20","30"),
                     minor_breaks=NULL) +
  scale_y_continuous(name="Mean MPSE score (NeoSeq; n=36)",
                     limits=c(0,16),
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


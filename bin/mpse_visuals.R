#!/usr/bin/env R

library(tidyverse)
library(ggwordcloud)


setwd("~/Documents/PhD/Yandell/software/MPSE")
set.seed(42)

rady <- read_delim("analysis/update/training_preds_ba1.0_sf1.0.tsv",
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
  mutate(hpo_cnt = str_count(hpo, "HP")) %>% 
  arrange(desc(pos_log_proba), neg_log_proba) %>% 
  mutate(rank_cumsum = cumsum(as.integer(diagnostic)),
         diag_rate = rank_cumsum / 1:nrow(.),
         list_fraction = 1:nrow(.) / nrow(.))


rady_cases <- rady %>% 
  select(seq_status, diagnostic, incidental, pos_log_proba, neg_log_proba) %>%
  filter(diagnostic=="1" & incidental=="0")
rady_ctrls <- rady %>% 
  select(seq_status, diagnostic, incidental, pos_log_proba, neg_log_proba) %>%
  filter(diagnostic=="0" | incidental=="1") %>% 
  slice_sample(n=208)
rady_sample <- bind_rows(rady_cases, rady_ctrls) %>% 
  arrange(desc(pos_log_proba), neg_log_proba) %>% 
  mutate(rank_cumsum = cumsum(as.integer(diagnostic)),
         diag_rate = rank_cumsum / 1:nrow(.),
         list_fraction = 1:nrow(.) / nrow(.))


rady_sample <- read_delim("analysis/update/training_preds_sample_ba1.0_sf1.0.tsv",
                           delim = "\t",
                           col_types = cols(
                             sex = col_factor(),
                             race = col_factor(),
                             ethnicity = col_factor(),
                             seq_status = col_factor(),
                             incidental = col_factor(),
                             class = col_factor()
                           )) %>% 
  arrange(desc(pos_log_proba), neg_log_proba) %>% 
  mutate(rank_cumsum = cumsum(as.integer(diagnostic)),
         diag_rate = rank_cumsum / 1:nrow(.),
         list_fraction = 1:nrow(.) / nrow(.))


utah_phy <- read_delim("analysis/update/physician_terms_preds.tsv",
                       delim = "\t", 
                       col_types = cols(
                         seq_status = col_factor(),
                         class = col_factor()
                       )) %>% 
  arrange(desc(pos_log_proba), neg_log_proba) %>% 
  mutate(hpo_cnt = str_count(hpo, "HP"),
         rank_cumsum = cumsum(as.integer(diagnostic)),
         diag_rate = rank_cumsum / 1:nrow(.),
         list_fraction = 1:nrow(.) / nrow(.))

utah_man <- read_delim("analysis/update/validation_preds.tsv",
                    delim = "\t", 
                    col_types = cols(
                      seq_status = col_factor(),
                      class = col_factor()
                    )) %>% 
  arrange(desc(pos_log_proba), neg_log_proba) %>% 
  mutate(hpo_cnt = str_count(hpo, "HP"),
         rank_cumsum = cumsum(as.integer(diagnostic)),
         diag_rate = rank_cumsum / 1:nrow(.),
         list_fraction = 1:nrow(.) / nrow(.))


utah_edw <- read_delim("analysis/update/NeoSeq_from_EDW.tsv",
                        delim = "\t", 
                        col_types = cols(
                          seq_status = col_factor(),
                          class = col_factor()
                        )) %>% 
  arrange(desc(pos_log_proba), neg_log_proba) %>% 
  mutate(hpo_cnt = str_count(hpo, "HP"),
         rank_cumsum = cumsum(as.integer(diagnostic)),
         diag_rate = rank_cumsum / 1:nrow(.),
         list_fraction = 1:nrow(.) / nrow(.))


utah_ctrl <- read_delim("analysis/utah_controls_n3000/control_preds.tsv",
                         delim = "\t", 
                         col_types = cols(
                           seq_status = col_factor(),
                           class = col_factor()
                         )) %>% 
  arrange(desc(pos_log_proba), neg_log_proba) %>% 
  mutate(hpo_cnt = str_count(hpo, "HP"),
         rank_cumsum = cumsum(as.integer(diagnostic)),
         diag_rate = rank_cumsum / 1:nrow(.),
         list_fraction = 1:nrow(.) / nrow(.)) %>% 
  filter(!pid %in% utah_edw$pid)


dens_gg_join <- bind_rows("rady_ctrl"=select(rady, scr, seq_status) %>% filter(seq_status=="0"),
                          "rady_case"=select(rady, scr, seq_status) %>% filter(seq_status=="1"),
                          "utah_man"=select(utah_man, scr) %>% mutate(seq_status="2"), 
                          #"utah_edw"=select(utah_edw, scr) %>% mutate(seq_status="3"),
                          "utah_ctrl"=select(utah_ctrl, scr) %>% mutate(seq_status="-1"),
                          .id = "cohort")
dens_gg_facet <- bind_rows("rady"=select(rady, scr, seq_status) %>% mutate(site="RCHSD"),
                          "utah_man"=select(utah_man, scr) %>% mutate(seq_status="1", site="UofU"), 
                          "utah_ctrl"=select(utah_ctrl, scr) %>% mutate(seq_status="0", site="UofU"),
                          .id = "cohort")


card_phen <- read_delim("analysis/test/cardinal_phenotypes.tsv",
                        delim = "\t", 
                        col_types = cols(
                          pid = col_character(),
                          term_id = col_character(),
                          term_name = col_character(),
                          coef = col_double()
                        ))


line_gg_sample <- ggplot(rady_sample, aes(x=list_fraction, y=diag_rate)) +
  geom_line(aes(color="RCHSD"), size=1) + 
  scale_y_continuous(name="Diagnostic rate", 
                     limits=c(0, 1.0),
                     breaks=c(0.00, 0.29, 0.36, 0.41, 0.50, 1.00),
                     labels=c("0%","29%","36%","41%","50%","100%"),
                     minor_breaks=NULL,
                     expand=c(0,0)) + 
  scale_x_continuous(name="Top scoring fraction of probands",
                     limits=c(0, 1.0),
                     breaks=c(0.00, 0.20, 0.50, 0.98),
                     labels=c("0%","20%","50%","100%"),
                     minor_breaks=NULL,
                     expand=c(0,0)) + 
  geom_line(data=utah_man, aes(x=list_fraction, y=diag_rate, color="NeoSeq-Manual"),
            size=1) +
  geom_line(data=utah_edw, aes(x=list_fraction, y=diag_rate, color="NeoSeq-Auto"),
            size=1) +
  scale_color_manual(name="Cohort", breaks=c("NeoSeq-Auto", "NeoSeq-Manual", "RCHSD"), values=c("NeoSeq-Auto"="green", "NeoSeq-Manual"="blue", "RCHSD"="red")) +
  geom_hline(yintercept=0.29, linetype="dotted") + 
  geom_hline(yintercept=0.36, linetype="dotted") + 
  geom_hline(yintercept=0.41, linetype="dotted") + 
  geom_vline(xintercept=0.20, alpha=0.3, size=0.7) + 
  theme_bw() + 
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.position=c(0.82, 0.82),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


line_gg_full <- ggplot(rady, aes(x=list_fraction, y=diag_rate)) +
  geom_smooth(aes(color="RCHSD"), method = lm, formula = y ~ poly(x, 25), se = FALSE, size=1) + 
  scale_y_continuous(name="Diagnostic rate", 
                     limits=c(0, 1.0),
                     breaks=c(0.00, 0.18, 0.35, 0.50, 1.00),
                     labels=c("0%","18%","35%","50%","100%"),
                     minor_breaks=NULL,
                     expand=c(0,0)) + 
  scale_x_continuous(name="Top scoring fraction of probands",
                     limits=c(0, 1.0),
                     breaks=c(0.00, 0.20, 0.50, 0.98),
                     labels=c("0%","20%","50%","100%"),
                     minor_breaks=NULL,
                     expand=c(0,0)) + 
  geom_line(data=utah_man, aes(x=list_fraction, y=diag_rate, color="NeoSeq-Manual"),
            size=1) +
  geom_line(data=utah_edw, aes(x=list_fraction, y=diag_rate, color="NeoSeq-Auto"),
            size=1) +
  scale_color_manual(name="Cohort", breaks=c("NeoSeq-Auto", "NeoSeq-Manual", "RCHSD"), values=c("NeoSeq-Auto"="green", "NeoSeq-Manual"="blue", "RCHSD"="red")) +
  geom_hline(yintercept=0.18, linetype="dashed") + 
  geom_hline(yintercept=0.35, linetype="dotted") + 
  geom_vline(xintercept=0.20, alpha=0.3, size=0.7) + 
  geom_segment(aes(x=0.7, xend=0.6, y=0.47, yend=0.37), 
               inherit.aes=FALSE, size=0.3, 
               arrow=arrow(length=unit(0.01, "npc"))) +
  annotate(geom="text", label="RCHSD Diag. Rate (35%)", x=0.8, y=0.5, size=4.5) +
  annotate(geom="text", label="NSIGHT2 rate of Mendelian", x=0.4, y=0.21, size=4.5) +
  annotate(geom="text", label="disease in NICU (18%)", x=0.4, y=0.15, size=4.5) +
  theme_bw() + 
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=12),
        legend.title=element_text(size=14),
        legend.text=element_text(size=12),
        legend.position=c(0.82, 0.82),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


dens_gg <- ggplot(dens_gg_join, aes(x=scr, fill=seq_status)) + 
  geom_density(alpha=0.5) + 
  geom_segment(aes(x=-75, xend=350, y=0, yend=0), inherit.aes=FALSE, size=0.3, alpha=0.5) +
  scale_y_continuous(name="Density", 
                     limits=c(0, 0.05),
                     minor_breaks=NULL) + 
  scale_x_continuous(name="MPSE Score",
                     limits=c(-75, 350),
                     breaks=seq(-60, 350, 30),
                     minor_breaks=NULL) +
  scale_fill_manual(values = c("#00BFC4","#F8766D","#00BA38","#C77CFF"),
                    guide = guide_legend(title=NULL),
                    labels = c("UofU Not Sequenced","RCHSD Not Sequenced","RCHSD Sequenced","UofU NeoSeq")) +
  theme_bw() + 
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=12),
        legend.text=element_text(size=12),
        legend.position=c(0.62, 0.27))


dens_gg_facet <- ggplot(dens_gg_facet, aes(x=scr, fill=seq_status)) + 
  geom_density(alpha=0.5) + 
  geom_segment(aes(x=-75, xend=350, y=0, yend=0), inherit.aes=FALSE, size=0.3, alpha=0.5) +
  scale_y_continuous(name="Density", 
                     limits=c(0, 0.05),
                     minor_breaks=NULL) + 
  scale_x_continuous(name="MPSE Score",
                     limits=c(-75, 350),
                     breaks=seq(-60, 350, 30),
                     minor_breaks=NULL) +
  scale_fill_manual(values = c("#00BFC4","#F8766D"),
                    guide = guide_legend(title=NULL),
                    labels = c("Not Sequenced", "Sequenced")) +
  theme_bw() + 
  facet_wrap(~site, ncol=1)


wordcloud_FUN <- function(data) {
  ggplot(data, aes(label=term_name, size=coef, color=coef)) + 
    geom_text_wordcloud() + 
    scale_color_gradient(low = "darkgreen", high = "red")
}

wordcloud_gg <- card_phen %>% 
  split(.$pid) %>% 
  map(wordcloud_FUN)


barplot_FUN <- function(data, nterm=10) {
  # see Hadley Wickham's 'Modern Graphics' for subplot description
  if (nrow(data) > 20) {
    data <- bind_rows(top_n(data, nterm, coef), top_n(data, nterm, -coef))
  }
  ggplot(data, aes(x=coef, y=reorder(term_name, coef, sum))) + 
    geom_bar(stat="identity") + 
    scale_x_continuous(name="Term Weight") + 
    scale_y_discrete(name="Cardinal Term")
}
barplot_viewport <- function(data) {
  ggplot(data, aes(x=coef, y=reorder(term_name, coef, sum))) + 
    geom_bar(stat="identity", width = 1.0) + 
    labs(x=NULL, y=NULL) +
    theme(axis.text=element_blank(),
          axis.ticks=element_blank(),
          panel.background=element_blank())
}

barplot_gg <- barplot_FUN(filter(card_phen, pid=="5"), nterm=10)
vp_gg <- barplot_viewport(filter(card_phen, pid=="5"))

pdf("card_pheno_viewport.pdf", width=10, height=6)
subvp <- viewport(width=0.3, height=0.35, x=0.75, y=0.3)
barplot_gg
print(vp_gg, vp=subvp)
dev.off()


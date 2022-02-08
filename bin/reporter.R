#!/usr/bin/env R

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(magrittr)


setwd("~/Documents/PhD/Yandell/software/MPSE")

desc_stats <- function(data) {
	return(table)
}



df <- read_delim("analysis/example/tables/training_predictions.csv", 
                 delim = "\t",
                 col_types = cols(
                   pid = col_character(),
                   sex = col_factor(),
                   race = col_factor(),
                   ethnicity = col_factor(),
                   age = col_double(),
                   seq_status = col_factor(),
                   diagnostic = col_factor(),
                   incidental = col_factor(),
                   neg_proba = col_double(),
                   pos_proba = col_double(),
                   classes = col_double()
                 )) %>% 
  mutate(classes = factor(classes),
         hpo_cnt = str_count(hpo, "HP:"),
         hpo_clean_cnt = str_count(hpo_clean, "HP:"),
         age_floor = floor(age / 5) * 5)


grouped_summary <- function(data, group, var) {
  summary <- data %>% 
    group_by(.data[[group]]) %>% 
    summarise(n = n(),
              min = min(.data[[var]]),
              max = max(.data[[var]]),
              mean = mean(.data[[var]]),
              q1 = quantile(.data[[var]], 0.25),
              median = quantile(.data[[var]], 0.50),
              q3 = quantile(.data[[var]], 0.75))
  return(summary)
}


age_by_seq <- grouped_summary(df, "seq_status", "age_floor")
hpo_by_seq <- grouped_summary(df, "seq_status", "hpo_clean_cnt")

hpo_by_age <- grouped_summary(df, "age_floor", "hpo_clean_cnt")
hpo_by_age %>% 
  select(age_floor, n, mean, median) %>% 
  pivot_longer(cols=c(mean, median), names_to="measure", values_to="value") %>% 
  ggplot(aes(x=log2(age_floor+1), y=value, colour=measure)) +
  geom_line()

df %>% 
  ggplot(aes(x=log2(age_floor+1), y=hpo_clean_cnt)) + 
  geom_smooth() + 
  geom_smooth(aes(colour=seq_status))



df %>% 
  ggplot(aes(x=age_floor, y=hpo_clean_cnt)) + 
  geom_smooth() + 
  geom_point(alpha=0.2, size=0.3)
df %>% 
  filter(age_floor <= 365) %>% 
  ggplot(aes(x=age_floor, y=hpo_clean_cnt)) + 
  geom_smooth() + 
  geom_point(alpha=0.2, size=0.3)


df %>% 
  ggplot(aes(x=hpo_cnt, colour=seq_status)) + 
  geom_histogram()
df %>% 
  ggplot(aes(x=hpo_clean_cnt, colour=seq_status)) + 
  geom_histogram()









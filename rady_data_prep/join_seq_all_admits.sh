#!/usr/bin/env bash

awk 'BEGIN{OFS=","}{printf("%s,%s", $0, NR>1?NR-2 RS:"all_fid" RS)}' ALL_ADMITS_EXP/data_views/rady_hpo_all_admits_remove_repeats.csv > ALL_ADMITS_EXP/data_views/rady_hpo_all_admits_remove_repeats_enumer.csv

awk 'BEGIN{OFS=","}{printf("%s,%s", $0, NR>1?NR-2 RS:"seq_fid" RS)}' SEQ_ADMITS/data_views/rady_hpo_seq_admits_CliniThink.csv > SEQ_ADMITS/data_views/rady_hpo_seq_admits_CliniThink_enumer.csv



# outer join of ALL_ADMITS_EXP and SEQ_ADMITS
csvjoin -c ResearchID --outer ALL_ADMITS_EXP/data_views/rady_hpo_all_admits_remove_repeats_enumer.csv SEQ_ADMITS/data_views/rady_hpo_seq_admits_CliniThink_enumer.csv | \
	csvsort -c CT_HPO_FileID > rady_hpo_all_seq_joined.csv

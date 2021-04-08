#!/usr/bin/env bash


# outer join of ALL_ADMITS_EXP and SEQ_ADMITS
csvjoin -c ResearchID --outer ALL_ADMITS_EXP/data_views/rady_hpo_all_admits_remove_repeats.csv SEQ_ADMITS/data_views/rady_hpo_seq_admits_CliniThink.csv | \
	csvsort -c CT_HPO_FileID \
	> rady_hpo_all_seq_joined.csv

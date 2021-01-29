#!/usr/bin/env bash

# make list of patient CSV files
ls ../../rady_data/RADY_HPO_SEQ_ADMITS | \
	grep 'UR.*csv' | \
	sed -re 's/_/,/;s/.csv//' > seq_file_list

# make list of sorted & concatenated CliniThink HPO terms from patient CSV files
ls ../../rady_data/RADY_HPO_SEQ_ADMITS/*.csv | \
	parallel sort | \
	cut -f1 -d, | \
	tr '\n' ':' | \
	sed -re 's/:Criterion:/\n/g;s/Criterion://' | \
	sed '$ s/.$/\n/' > seq_hpo_list

# concatenate two lists horizontally
paste -d, seq_file_list seq_hpo_list > seq_combined_list.csv

# add header row
sed -i 1i"ResearchID,seq_file_tag,seq_CliniThink_HPO" seq_combined_list.csv

# make output directory
mkdir data_views

# clean up the patient data file
in2csv -f csv ../../rady_data/RADY_HPO_SEQ_ADMITS/Utah_Rady_CT_HPO_Manifest_16APR2020_DeIdentified.txt | \
	csvgrep -c 1 -r "^$" -i > data_views/Utah_Rady_CT_HPO_Manifest_16APR2020_DeIdentified.csv

# join the CliniThink term list to the patient data file
# Manual_HPO must be removed because it contains "," and screws up formatting
csvjoin -c ResearchID --outer data_views/Utah_Rady_CT_HPO_Manifest_16APR2020_DeIdentified.csv seq_combined_list.csv | \
	csvcut -C Manual_HPO,ResearchID2 > data_views/rady_hpo_seq_admits_combined.csv

# recode years as days
grep years data_views/rady_hpo_seq_admits_combined.csv | sed -re 's/ years//' | awk -F, -v OFS=, '$2*=365' > years2days

csvgrep -c Age -m 'years' -i data_views/rady_hpo_seq_admits_combined.csv | \
	sed -re 's/ days//;s/Misisng DoB//' > just_days

cat just_days years2days > all_ages

csvsort -c ResearchID all_ages > data_views/rady_hpo_seq_admits_combined.csv

# remove temp files
rm seq_file_list seq_hpo_list seq_combined_list.csv years2days just_days all_ages




#---- Data Views ----#

# subset unique CliniThink term sets
sed "2d" data_views/rady_hpo_seq_admits_combined.csv | \
	csvgrep -c seq_CliniThink_HPO -r "^$" -i > data_views/rady_hpo_seq_admits_CliniThink.csv

# subset Positive Dx
csvgrep -c Positive -m "False" -i data_views/rady_hpo_seq_admits_combined.csv > data_views/rady_hpo_seq_admits_PositiveDx.csv

# subset Negative Dx
csvgrep -c Negative -m "False" -i data_views/rady_hpo_seq_admits_combined.csv > data_views/rady_hpo_seq_admits_NegativeDx.csv

# subset VUS
csvgrep -c VUS -m "False" -i data_views/rady_hpo_seq_admits_combined.csv > data_views/rady_hpo_seq_admits_VUS.csv

# subset Incidental
csvgrep -c Incidental -m "False" -i data_views/rady_hpo_seq_admits_combined.csv > data_views/rady_hpo_seq_admits_Incidental.csv

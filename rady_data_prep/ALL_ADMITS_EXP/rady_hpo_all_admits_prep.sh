#!/usr/bin/env bash

# make list of patient CSV files
ls ../../rady_data/RADY_HPO_ALL_ADMITS_EXP | \
	grep 'UR20.*csv' | \
	sed -re 's/_/,/;s/.csv//' > all_file_list

# make list of sorted & concatenated CliniThink HPO terms from patient CSV files
ls ../../rady_data/RADY_HPO_ALL_ADMITS_EXP/*.csv | \
	parallel sort | \
	cut -f1 -d, | \
	tr '\n' ':' | \
	sed -re 's/:Criterion:/\n/g;s/Criterion://' | \
	sed '$ s/.$/\n/' > all_hpo_list

# concatenate two lists horizontally
paste -d, all_file_list all_hpo_list > all_combined_list.csv

# add header row
sed -i 1i"CT_HPO_FileID,all_file_tag,all_CliniThink_HPO" all_combined_list.csv

# make output directory
mkdir data_views

# clean up the patient data file
in2csv -f csv ../../rady_data/RADY_HPO_ALL_ADMITS_EXP/CT_HPO_Manifest_Expanded_De_Identified.txt > data_views/CT_HPO_Manifest_Expanded_De_Identified.csv

# replace header
new_head="CT_HPO_FileID,Discharge_Date,Unit,Repeat_cases,Gender,Race,Ethnicity,Age5d,Accession_for_Seq,ResearchID"
sed -i "1s/.*/$new_head/" data_views/CT_HPO_Manifest_Expanded_De_Identified.csv

# join the CliniThink term list to the patient data file
csvjoin -c CT_HPO_FileID --outer data_views/CT_HPO_Manifest_Expanded_De_Identified.csv all_combined_list.csv | \
	csvcut -C CT_HPO_FileID2 | \
	csvsort -c CT_HPO_FileID > data_views/rady_hpo_all_admits_combined.csv

# remove temp files
rm all_file_list all_hpo_list all_combined_list.csv




#---- Data Views ----#

# subset repeat cases
csvgrep -c Repeat_cases -r "^$" -i data_views/rady_hpo_all_admits_combined.csv > data_views/rady_hpo_all_admits_repeat_cases.csv

# subset accessioned for sequencing
csvgrep -c Accession_for_Seq -r "^$" -i data_views/rady_hpo_all_admits_combined.csv > data_views/rady_hpo_all_admits_accessioned.csv

# subset those not accessioned for sequencing
csvgrep -c Accession_for_Seq -r "^$" data_views/rady_hpo_all_admits_combined.csv > data_views/rady_hpo_all_admits_not_accessioned.csv

# remove all duplicate patient rows - keep most recent admit
csvsort -c Repeat_cases,Discharge_Date -r data_views/rady_hpo_all_admits_repeat_cases.csv | \
	uniq -s29 -w2 > remove_dup
csvgrep -c Repeat_cases -r "^$" data_views/rady_hpo_all_admits_combined.csv > no_repeat
csvstack remove_dup no_repeat | \
	csvsort -c CT_HPO_FileID > data_views/rady_hpo_all_admits_remove_repeats.csv
rm remove_dup no_repeat

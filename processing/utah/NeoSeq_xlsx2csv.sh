#!/usr/bin/env bash

ml miniconda3/latest

in2csv NeoSeq_for_Bennet.xlsx --write-sheets -

csvcut -c MRN NeoSeq_for_Bennet_0.csv | \
	grep -E "[0-9]{8}" | \
	awk -v OFS=, '{print $0, "1"}' | \
	sed 1i"mrn,screened" > screened.csv

grep -E "[0-9]{8}" NeoSeq_for_Bennet_1.csv | \
	awk -v OFS=, '{print $0, "1"}' | \
	sed 1i"mrn,sequenced" > sequenced.csv
echo "21271499,1" >> sequenced.csv

csvcut -c1 NeoSeq_for_Bennet_2.csv | \
	grep -E "[0-9]{8}" | \
	awk -v OFS=, '{print $0, "1"}' | \
	sed 1i"mrn,diagnostic" > diagnostic.csv

#csvcut -c2 NeoSeq_for_Bennet_2.csv | \
#	grep -E "[0-9]{8}" | \
#	awk -v OFS=, '{print $0, "1"}' | \
#	sed 1i"mrn,nondiagnostic" > nondiagnostic.csv

csvjoin -c mrn --left -I sequenced.csv diagnostic.csv | \
	csvformat -T > status.tsv

in2csv NeoSeq_for_Bennet_HPOs.xlsx --sheet "Enrolled in NeoSeq" | \
	csvcut -c MRN,HPOs | \
	csvgrep -c MRN -r "\d{8}" | \
	csvformat -T > curated_HPOs.tsv

csvjoin -c MRN,mrn --left -I curated_HPOs.tsv diagnostic.csv | \
	csvformat -T > neoseq_curated_ref.tsv

#csvcut -t -c DW_PID,PAT_ID demographs.tsv | \
#	csvjoin -c PAT_ID,mrn -I - status.tsv | \
#	csvformat -T > neoseq_cross_demographs.tsv

csvcut -t -c DW_PID,PAT_ID NICU_patients.tsv | \
	csvjoin -c PAT_ID,mrn -I - status.tsv | \
	csvformat -T > neoseq_cross_demographs.tsv

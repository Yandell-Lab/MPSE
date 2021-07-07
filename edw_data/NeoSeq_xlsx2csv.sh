#!/usr/bin/env bash

ml miniconda3/latest

in2csv NeoSeq_for_Bennet.xlsx --write-sheets -

#csvcut -c MRN NeoSeq_for_Bennet_0.csv | \
#	grep -E "[0-9]{8}" | \
#	awk -v OFS=, '{print $0, "1"}' | \
#	sed 1i"mrn,screened" > neoseq_screened.csv

grep -E "[0-9]{8}" NeoSeq_for_Bennet_1.csv | \
	awk -v OFS=, '{print $0, "1"}' | \
	sed 1i"mrn,sequenced" > neoseq_sequenced.csv

csvcut -c1 NeoSeq_for_Bennet_2.csv | \
	grep -E "[0-9]{8}" | \
	awk -v OFS=, '{print $0, "1"}' | \
	sed 1i"mrn,diagnostic" > neoseq_diagnostic.csv

#csvcut -c2 NeoSeq_for_Bennet_2.csv | \
#	grep -E "[0-9]{8}" | \
#	awk -v OFS=, '{print $0, "1"}' | \
#	sed 1i"mrn,nondiagnostic" > neoseq_nondiagnostic.csv

csvjoin -c mrn --left -I neoseq_sequenced.csv neoseq_diagnostic.csv > EDW_IDs.csv

csvcut -t -c DW_PID,PAT_ID demographs.tsv | \
	csvjoin -c PAT_ID,mrn -I - EDW_IDs.csv > EDW_NeoSeq_IDs.csv



in2csv NeoSeq_for_Bennet_HPOs.xlsx --sheet "Enrolled in NeoSeq" | \
	csvcut -c MRN,HPOs | \
	csvgrep -c MRN -r "\d{8}" | \
	csvformat -T > tmp_NeoSeq_HPOs.csv

csvjoin -c MRN,mrn --left -I tmp_NeoSeq_HPOs.csv neoseq_diagnostic.csv | \
	csvformat -T > NeoSeq_HPOs.csv



rm NeoSeq_for_Bennet_*.csv
rm neoseq_*.csv
rm EDW_IDs.csv tmp_NeoSeq_HPOs.csv

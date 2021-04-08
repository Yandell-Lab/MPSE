#!/usr/bin/env python

import os
import sys
import argparse
import re
from collections import Counter
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'
from scipy.special import bdtr
from pyhpo.ontology import Ontology


def reader(ftag):
	with open(ftag, 'r') as f:
		dat = pd.read_csv(f, dtype={
			"Age":'Int64', "Age5d":'Int64', 
			"seq_fid":"Int64", "all_fid":"Int64"
			})
	return dat

def grp_count(dat, col_name):
	tgrp = dat.loc[dat[col_name].notnull()].shape[0]
	bgrp = dat.loc[dat[col_name].isnull()].shape[0]
	return tgrp, bgrp

def hpo_regex(string):
	hpo_reg = re.compile(r'HP:\d{7}')
	return hpo_reg.findall(string)

def counter(dat, col_name, grp):
	if grp == "t":
		sub = dat.loc[dat[col_name].notnull()]
		sub["hpo_format"] = sub["seq_HPO_clean"].apply(hpo_regex)
		lst = sub["hpo_format"].tolist()
	else:
		sub = dat.loc[dat[col_name].isnull()]
		sub["hpo_format"] = sub["all_HPO_clean"].apply(hpo_regex)
		lst = sub["hpo_format"].tolist()
	cnt = Counter([item for sublst in lst for item in sublst])
	return cnt

def printer(dat, count, tn, bn):
	count.to_csv(sys.stdout, index=False, header=False, sep="\t")

	print('\t'.join(["#", "NUM_B", str(bn)]))
	print('\t'.join(["#", "NUM_T", str(tn)]))

	def get_line(line):
		if str(line["seq_HPO_clean"]) != "nan":
			grp = "t"
			fid = line["seq_fid"]
			cnt = len(hpo_regex(line["seq_HPO_clean"]))
			age = line["Age"]
		else:
			grp = "b"
			fid = line["all_fid"]
			cnt = len(hpo_regex(line["all_HPO_clean"]))
			age = line["Age5d"]
		msg = "\t".join(["#", "CNT", grp, str(fid), str(cnt), str(age)])
		print(msg)

	keep = ["all_fid","seq_fid",
			"Age5d","Age",
			"all_HPO_clean","seq_HPO_clean"]

	dat[keep].apply(get_line, axis=1)


def main():
	os.chdir("/scratch/ucgd/lustre-work/yandell/u1323262/MPSE")
	
	parser = argparse.ArgumentParser()
	parser.add_argument('-d', '--data', action='store', required=True,
			help='CSV file containing data for targets and background.')
	args = parser.parse_args()

	df = reader(args.data)
	tn, bn = grp_count(df, "seq_HPO_clean")

	tcnt = counter(df, "seq_HPO_clean", "t")
	bcnt = counter(df, "seq_HPO_clean", "b")
	cnts = pd.DataFrame({
		"tcnt":pd.Series(tcnt), 
		"bcnt":pd.Series(bcnt)
		}).fillna(value=0).convert_dtypes()

	cnts["tfrq"] = cnts["tcnt"] / tn
	cnts["bfrq"] = cnts["bcnt"] / bn
	cnts["up"] = np.where(cnts["tfrq"] > cnts["bfrq"], "t", "b")
	
	cnts["binom"] = cnts.apply(
			lambda row: bdtr(row["bcnt"], bn, row["tfrq"]) if row["up"] == "t" \
					else bdtr(row["tcnt"], tn, row["bfrq"]),
					axis=1
					)

	cnts.reset_index(inplace=True)
	cnts.rename(columns={"index":"term"}, inplace=True)
	
	_ = Ontology()
	cnts["name"] = cnts.apply(
			lambda row: Ontology.get_hpo_object(row["term"]).name,
			axis=1
			)

	order_cnts = cnts[["tcnt","bcnt","binom","up","term","name"]]
	printer(df, order_cnts, tn, bn)


if __name__ == '__main__':
	main()

#!/usr/bin/env python

import sys
import re
import numpy as np
import pandas as pd
from pyhpo.ontology import Ontology
from pyhpo.set import HPOSet


all_f = "ALL_ADMITS_EXP/data_views/rady_hpo_all_admits_remove_repeats.csv"
with open(all_f, "r") as f:
	all_df = pd.read_csv(f)

seq_f = "SEQ_ADMITS/data_views/rady_hpo_seq_admits_CliniThink.csv"
with open(seq_f, "r") as f:
	seq_df = pd.read_csv(f)

def hpo_parse(hpo_str):
	hpo_reg = re.compile(r'hp\d{7}')
	hpo_split = [x.replace('hp', 'HP:', 1) for x in hpo_reg.findall(hpo_str)]
	return hpo_split

def subset_child(hpo_lst):
	hpo_set = HPOSet.from_queries(hpo_lst)
	hpo_subset = hpo_set.child_nodes()
	subset_dic = hpo_subset.toJSON()
	out_str = "_".join([x["id"] for x in subset_dic])
	return out_str



def main():
	_ = Ontology()
	all_df["all_HPO_clean"] = all_df.apply(
		lambda row: subset_child(hpo_parse(row["all_CliniThink_HPO"])),
		axis=1
		)
	seq_df["seq_HPO_clean"] = seq_df.apply(
		lambda row: subset_child(hpo_parse(row["seq_CliniThink_HPO"])),
		axis=1
		)
	df = all_df.merge(seq_df, on="ResearchID", how="outer")
	df.sort_values(by=["CT_HPO_FileID"]).to_csv("rady_hpo_all_seq_joined.csv", index=False)

if __name__ == "__main__":
	main()

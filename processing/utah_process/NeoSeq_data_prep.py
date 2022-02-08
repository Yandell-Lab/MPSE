#!/usr/bin/env python3

import sys
import re
import numpy as np
import pandas as pd
from hashlib import sha256
from pyhpo.ontology import Ontology
from pyhpo.set import HPOSet


#----------- HPO from EDW -----------#
#------------------------------------#
edw_tag = "/scratch/ucgd/lustre-work/yandell/u1323262/MPSE/edw_data/EDW_NeoSeq_IDs.csv"
ids_tag = "/scratch/ucgd/lustre-work/yandell/u1323262/FTP/data/edw_data_lemmon/hpo_ids.tsv"

edw = pd.read_csv(edw_tag, sep=",", usecols=["DW_PID","sequenced","diagnostic"])
ids = pd.read_csv(ids_tag, sep="\t", usecols=["mrn","terms"])

def hasher(dw_pid):
	m = sha256()
	m.update(int(dw_pid).to_bytes(4, "big"))
	return m.hexdigest()

edw["mrn"] = edw["DW_PID"].apply(hasher)

df = edw.merge(ids, on="mrn", how="inner").fillna(0)

def hpo_parse(hpo_str):
	hpo_reg = re.compile(r'HP:\d{7}')
	hpo_split = hpo_reg.findall(hpo_str)
	return hpo_split

def subset_child(hpo_lst):
	hpo_set = HPOSet.from_queries(hpo_lst)
	hpo_subset = hpo_set.child_nodes()
	subset_dic = hpo_subset.toJSON()
	out_str = "_".join([x["id"] for x in subset_dic])
	return out_str

_ = Ontology()
df["terms_clean"] = df["terms"].apply(lambda x: subset_child(hpo_parse(x)))

onehot = df.join(df["terms_clean"].str.get_dummies(sep="_"))

seq_samples = onehot.loc[onehot["sequenced"] == 1].drop(columns=["DW_PID","sequenced","terms","terms_clean"])
#diag_samples = onehot.loc[onehot["diagnostic"] == 1]

seq_samples.to_csv("/scratch/ucgd/lustre-work/yandell/u1323262/MPSE/edw_data/EDW_NeoSeq_validation_cases.csv", index=False)



#---------- HPO from NeoSeq ----------#
#-------------------------------------#
neo_tag = "/scratch/ucgd/lustre-work/yandell/u1323262/MPSE/edw_data/NeoSeq_HPOs.csv"
neo = pd.read_csv(neo_tag, sep="\t")

neo["terms_clean"] = neo["HPOs"].apply(lambda x: subset_child(hpo_parse(x)))
neo_onehot = neo.join(neo["terms_clean"].str.get_dummies(sep="_"))

neo_onehot.to_csv("/scratch/ucgd/lustre-work/yandell/u1323262/MPSE/edw_data/NeoSeq_validation_cases.csv", index=False)

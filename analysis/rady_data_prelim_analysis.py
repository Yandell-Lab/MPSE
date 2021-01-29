#!/usr/bin/env python

import os
import sys
import math
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import seaborn as sns


os.chdir("/scratch/ucgd/lustre-work/yandell/u1323262/MPSE")

all_path = "rady_data_prep/ALL_ADMITS_EXP/data_views/rady_hpo_all_admits_remove_repeats.csv"
seq_path = "rady_data_prep/SEQ_ADMITS/data_views/rady_hpo_seq_admits_CliniThink.csv"

all_keep = ["CT_HPO_FileID","Age5d","ResearchID","all_CliniThink_HPO"]
seq_keep = ["ResearchID","Age","seq_CliniThink_HPO"]

all_data = pd.read_csv(all_path, sep=",", index_col="CT_HPO_FileID", usecols=all_keep)
seq_data = pd.read_csv(seq_path, sep=",", index_col="ResearchID", usecols=seq_keep)

all_data = all_data[all_data["ResearchID"].isnull()]
all_data.drop(columns="ResearchID", inplace=True)

all_data.rename(columns={"Age5d":"Age", "all_CliniThink_HPO":"HPO"}, inplace=True)
seq_data.rename(columns={"seq_CliniThink_HPO":"HPO"}, inplace=True)

seq_data = seq_data[seq_data["Age"].notnull()]
seq_data = seq_data.astype({"Age":"int64"})

data = pd.concat([all_data, seq_data], axis=0)

data["hpo_count"] = data["HPO"].map(lambda x: len(x.split(":")))
data["age_log"] = data["Age"].map(lambda x: math.log10(1+x))

ols_fit = smf.ols("hpo_count ~ Age", data=data).fit()
#print(ols_fit.rsquared)

#scat = sns.regplot(x="Age", y="hpo_count", data=data, ci=None)
#fig = scat.get_figure()
#fig.savefig("analysis/rady_age_vs_hpo_count_scatter.png")

scat1 = sns.regplot(x="age_log", y="hpo_count", data=data, ci=None)
fig = scat1.get_figure()
fig.savefig("analysis/rady_age_vs_hpo_count_scatter_log.png")

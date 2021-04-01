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
seq_keep = ["ResearchID","Age","Positive","seq_CliniThink_HPO"]

all_data = pd.read_csv(all_path, sep=",", index_col="CT_HPO_FileID", usecols=all_keep)
seq_data = pd.read_csv(seq_path, sep=",", index_col="ResearchID", usecols=seq_keep)

# remove all sequenced patients
all_data = all_data[all_data["ResearchID"].isnull()]
all_data.drop(columns="ResearchID", inplace=True)

all_data.rename(columns={"Age5d":"Age", "all_CliniThink_HPO":"HPO"}, inplace=True)
seq_data.rename(columns={"seq_CliniThink_HPO":"HPO"}, inplace=True)

# remove 1 observation with missing age
seq_data = seq_data[seq_data["Age"].notnull()]
seq_data = seq_data.astype({"Age":"int64"})

all_data["group"] = "control"
seq_data["group"] = np.where(seq_data["Positive"], "Dx", "case")
seq_data.drop(columns="Positive", inplace=True)

data = pd.concat([all_data, seq_data], axis=0)

# create HPO term counts for each observation
data["hpo_count"] = data["HPO"].map(lambda x: len(x.split(":")))
# create log2 transform of age for plotting
data["age_log"] = data["Age"].map(lambda x: math.log2(1+x))

# remove high ages
age_limit = 12000 #days
data_rm_ol = data[data["Age"]<age_limit]

def scat_plotter(x, y, dat, group=None, f=None):
	fit = smf.ols(y + " ~ " + x, data=dat).fit()
	if f:
		if group:
			p = sns.lmplot(x=x, y=y, data=dat, hue=group, ci=None)
			p.set(xlabel='Age (log2)', ylabel='Distinct HPO Term Counts')
			p.savefig("analysis/" + f)
		else:
			p = sns.regplot(x=x, y=y, data=dat, ci=None)
			p.set(xlabel='Age (log2)', ylabel='Distinct HPO Term Counts')
			fig = p.get_figure()
			fig.savefig("analysis/" + f)
		print("r-squared for {0}: {1}".format(f, str(fit.rsquared)))

#scat_plotter("Age", "hpo_count", data, f="prelim_plots/rady_age_vs_hpo_count_scatter.png")
#scat_plotter("Age", "hpo_count", data_rm_ol, f="prelim_plots/rady_age_vs_hpo_count_scatter_remove_outlier.png")
scat_plotter("age_log", "hpo_count", data, f="prelim_plots/rady_log2_age_vs_hpo_count_scatter.png")
#scat_plotter("age_log", "hpo_count", data_rm_ol, f="prelim_plots/rady_log2_age_vs_hpo_count_scatter_remove_outlier.png")
#scat_plotter("age_log", "hpo_count", data, group="group", f="prelim_plots/rady_log2_age_vs_hpo_count_group_scatter.png")
#scat_plotter("age_log", "hpo_count", data_rm_ol, group="group", f="prelim_plots/rady_log2_age_vs_hpo_count_group_scatter_remove_outlier.png")

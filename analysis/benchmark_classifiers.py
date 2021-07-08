#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.model_selection import train_test_split
#from sklearn.model_selection import RepeatedKFold
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import cross_validate
from sklearn.model_selection import cross_val_predict
from sklearn.naive_bayes import BernoulliNB
from sklearn import metrics


rady_fh = "/scratch/ucgd/lustre-work/yandell/u1323262/MPSE/rady_data_prep/rady_hpo_all_seq_joined.csv"
edw_fh = "/scratch/ucgd/lustre-work/yandell/u1323262/MPSE/edw_data/EDW_NeoSeq_validation_cases.csv"
neo_fh = "/scratch/ucgd/lustre-work/yandell/u1323262/MPSE/edw_data/NeoSeq_validation_cases.csv"

raw = pd.read_csv(rady_fh, dtype={"Positive": np.int32})
edw = pd.read_csv(edw_fh, dtype={"diagnostic": np.int32}, nrows=9)
neo = pd.read_csv(neo_fh)
neo["diagnostic"] = neo["diagnostic"].fillna(0, downcast="infer")

bcols = {
        "CT_HPO_FileID": "pid", 
        "all_HPO_clean": "terms", 
        }
tcols = {
        "ResearchID": "pid", 
        "seq_HPO_clean": "terms"
        }
# split rady data into background & target dataframes
bdat = raw.loc[raw["ResearchID"].isnull()][["CT_HPO_FileID", "all_HPO_clean"]].rename(columns=bcols)
tdat = raw.loc[raw["ResearchID"].notnull()][["ResearchID", "seq_HPO_clean", "Positive"]].rename(columns=tcols)

# one-hot-encode HPO terms
b_onehot = bdat.join(bdat["terms"].str.get_dummies(sep="_"))
t_onehot = tdat.join(tdat["terms"].str.get_dummies(sep="_"))
b_onehot["outcome"] = 0
t_onehot["outcome"] = 1

# recombine background & target data
rad = b_onehot.merge(t_onehot, how="outer", on=None).fillna(0, downcast="infer")

rad_terms = rad.columns[pd.Series(rad.columns).str.startswith("HP:")]
edw_terms = edw.columns[pd.Series(edw.columns).str.startswith("HP:")]
neo_terms = neo.columns[pd.Series(neo.columns).str.startswith("HP:")]
rad_X = rad[rad_terms]
edw_X = edw[edw_terms]
neo_X = neo[neo_terms]

concat_X = pd.concat([rad_X, edw_X, neo_X], keys=["rady","edw","neo"])
concat_X = concat_X.loc[:,concat_X.iloc[0,:].notna()].fillna(0, downcast="infer")
concat_terms = concat_X.columns[pd.Series(concat_X.columns).str.startswith("HP:")]

rad_X = concat_X.loc[["rady"]].reset_index(drop=True).to_numpy()
edw_X = concat_X.loc[["edw"]].reset_index(drop=True).to_numpy()
neo_X = concat_X.loc[["neo"]].reset_index(drop=True).to_numpy()

rad_y = rad["outcome"].to_numpy()

mod = BernoulliNB()
fit = mod.fit(rad_X, rad_y)

edw_pre = fit.predict_proba(edw_X)
neo_pre = fit.predict_proba(neo_X)
#edw_scr = -np.log(edw_pre[:,0])
#neo_scr = -np.log(neo_pre[:,0])

edw_pred = pd.DataFrame(edw_pre, columns=["neg_proba","pos_proba"]).join(edw[["diagnostic"]])
neo_pred = pd.DataFrame(neo_pre, columns=["neg_proba","pos_proba"]).join(neo[["diagnostic"]])
edw_pred.to_csv("/scratch/ucgd/lustre-work/yandell/u1323262/MPSE/analysis/edw_case_predictions.csv", index=False)
neo_pred.to_csv("/scratch/ucgd/lustre-work/yandell/u1323262/MPSE/analysis/neo_case_predictions.csv", index=False)



# cross-validation
loo = LeaveOneOut()

# Bernoulli naive bayes model
clf = BernoulliNB()
#scores = cross_validate(clf, X, y, 
#        cv=loo,
#        scoring=["accuracy"],#,"roc_auc"],
#        return_train_score=True,
#        n_jobs=8)
#
#print("Average model accuracy across validation sets")
#print(scores["test_accuracy"].mean())

# generate class probabilities predictions
predictions = cross_val_predict(clf, X, y, cv=loo, method="predict_proba", n_jobs=8)

# compute AUC & plot ROC
#fpr, tpr, thresholds = metrics.roc_curve(y, predictions[:,1], pos_label=1)
#roc_auc = metrics.auc(fpr, tpr)
#roc_plot = metrics.RocCurveDisplay(fpr=fpr, tpr=tpr, roc_auc=roc_auc)
#roc_plot.plot()
#plt.show()
#plt.clf()

# join probabilities with positive diagnosis indicator
pred_df = pd.DataFrame(predictions, columns=["neg_proba", "pos_proba"]).join(rad[["pid", "Positive", "outcome"]])
# calculate log-probabilities and create rank indicator
pred_df["neg_score"] = -np.log(pred_df["neg_proba"])
pred_df["score_rank"] = pred_df["neg_score"].rank(method="first", ascending=False).astype("int32")

# plot kernel density of log-probabilities
#dist_plot = sns.kdeplot(pred_df["pos_score"], bw_adjust=0.5)
#plt.savefig("pred_log_proba_kde.png")
#plt.clf()

# split into sequenced and non-sequenced subgroups for sampling
pos_samples = pred_df.loc[pred_df["Positive"]==1,]
neg_samples = pred_df.loc[pred_df["Positive"]==0,]

# sample from non-sequenced subgroup to give 18% overall diagnostic rate
diag_rate_sample = pos_samples.merge(neg_samples.sample(n=388, random_state=42), how="outer", on=None)
# re-calculate ranks as a result of sampling
diag_rate_sample["score_rank"] = diag_rate_sample["neg_score"].rank(method="first", ascending=False).astype("int32")
# cumulative sum of positive diagnoses across the rank list
diag_rate_sample["rank_cumsum"] = diag_rate_sample.sort_values(by=["score_rank"])["Positive"].cumsum()

# compute diagnostic rate and rank percentile
diag_rate_sample["diag_rate"] = diag_rate_sample["rank_cumsum"] / diag_rate_sample["score_rank"]
diag_rate_sample["list_fraction"] = diag_rate_sample["score_rank"] / diag_rate_sample.shape[0]


pred_df.to_csv("bernoulli_nb_predictions.csv")
diag_rate_sample.to_csv("bernoulli_nb_diagnostic_rate.csv")

#line_plot = sns.lineplot(data=diag_rate_sample, x="list_fraction", y="diag_rate")
#plt.savefig("diagnostic_rate_lineplot.png")
#plt.clf()

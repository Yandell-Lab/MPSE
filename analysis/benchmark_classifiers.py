#!/usr/bin/env python3

import sys
from time import time
import logging
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.model_selection import train_test_split
#from sklearn.model_selection import RepeatedKFold
from sklearn.model_selection import LeaveOneOut
from sklearn.model_selection import cross_validate
from sklearn.model_selection import cross_val_predict
from sklearn.naive_bayes import BernoulliNB
from sklearn import metrics


# Display progress logs on stdout
logging.basicConfig(
        level=logging.INFO, 
        format="%(asctime)s %(levelname)s %(message)s"
        )

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--data", action="store", required=True)
args = parser.parse_args()

keep_cols = {
        "CT_HPO_FileID": "background", 
        "ResearchID": "target", 
        "all_HPO_clean": "b_terms", 
        "seq_HPO_clean": "t_terms"
        }
raw = pd.read_csv(args.data, usecols=keep_cols.keys())

dat = raw.rename(columns=keep_cols)

bdat = dat.loc[dat["target"].isnull()][["background", "b_terms"]]
tdat = dat.loc[dat["target"].notnull()][["target", "t_terms"]]

b_onehot = bdat["b_terms"].str.get_dummies(sep="_")
t_onehot = tdat["t_terms"].str.get_dummies(sep="_")

b_onehot["outcome"] = 0
t_onehot["outcome"] = 1

df = b_onehot.merge(t_onehot, how="outer", on=None).fillna(0, downcast="infer")

y = df.pop("outcome").to_numpy()
X = df.to_numpy()
term_cols = df.columns

seed=42
#rkf = RepeatedKFold(n_splits=10, n_repeats=5, random_state=seed)
loo = LeaveOneOut()

clf = BernoulliNB()
scores = cross_validate(clf, X, y, 
        cv=loo,
        scoring=["accuracy"],#,"roc_auc"],
        return_train_score=True,
        n_jobs=8)

predictions = cross_val_predict(clf, X, y, cv=loo, method="predict_proba", n_jobs=8)
fpr, tpr, thresholds = metrics.roc_curve(y, predictions[:,1], pos_label=1)
roc_auc = metrics.auc(fpr, tpr)
display = metrics.RocCurveDisplay(fpr=fpr, tpr=tpr, roc_auc=roc_auc)
display.plot()
plt.show()


print(scores["test_accuracy"].mean())



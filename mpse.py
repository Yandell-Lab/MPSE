#!/usr/bin/env python3

import sys
import datetime as dt
import argparse
import random
import csv 
from processing.rady_process.rady_data_prep import ready

from pyhpo.ontology import Ontology
from pyhpo.set import HPOSet

import numpy as np
import pandas as pd

from sklearn.model_selection import LeaveOneOut, cross_validate, cross_val_predict
from sklearn.naive_bayes import BernoulliNB
from sklearn import metrics


def argue():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--training", 
            default="processing/rady_process/rady_training_data.csv", 
            help="Case/control training data in standard format.")
    parser.add_argument("-v", "--validate",
            default="processing/utah_process/utah_validation_data.csv",
            help="Validation data in standard format.")
    parser.add_argument("-o", "--outdir", 
            default="analysis/example/",
            help="Output directory for results & reports.")
    return parser.parse_args()


def get_col_position(data, name):
    index = data[0].index(name)
    return index


def child_terms(hpo):
    hpo_lst = hpo.split(";")
    hpo_set = HPOSet.from_queries(hpo_lst)
    hpo_subset = hpo_set.child_nodes()
    subset_dic = hpo_subset.toJSON()
    out_str = ";".join([x["id"] for x in subset_dic])
    return out_str


def quality_report(data):
    pass


def lst2array(data):
    df = pd.DataFrame(data[1:], columns=data[0])
    one_hot = df["hpo"].str.get_dummies(sep=";")
    #terms = one_hot.columns
    X = one_hot.to_numpy()
    y = df["seq_status"].to_numpy()
    return X,y


def trainer(X, y, mod=BernoulliNB(), cv=LeaveOneOut()):
    scores = cross_validate(mod, X, y,
            cv=cv,
            scoring=["accuracy",],
            return_train_score=True,
            n_jobs=-1)

    predictions = cross_val_predict(mod, X, y, 
            cv=cv, 
            method="predict_log_proba", 
            n_jobs=-1)
    return predictions


def sample_cohort(data, diagnos_rate=0.18):
    cases_idx = get_col_position(data, "diagnostic")
    incident_idx = get_col_position(data, "incidental")

    cases = [x for x in data if x[cases_idx]=="1" and x[incident_idx]=="0"]
    controls = [x for x in data if x[cases_idx]!="1"]

    case_n = len(cases)
    control_n = len(controls)
    n = len(data)-1

    if case_n / n < diagnos_rate:
        control_n = round(case_n / diagnos_rate)
        control_sample = random.sample(controls, control_n)
        sample = [data[0]] + cases + control_sample
    else:
        case_n = round(diagnos_rate * control_n / (1.0 - diagnos_rate))
        case_sample = random.sample(cases, case_n)
        sample = [data[0]] + case_sample + controls
    return sample


def rank_list(data):
    df = pd.DataFrame(data, columns=["neg_log_proba","pos_log_proba"])
    df["score_rank"] = df["neg_log_proba"].rank(method="first", ascending=True)
    return df


def main():
    args = argue()
    #_ = Ontology()

    train = ready(args.training, delim="\t")
    #valid = ready(args.validate, delim="\t")

    #hpo_idx = get_col_position(train, "hpo")
    #for row in train[1:]:
    #    row[hpo_idx] = child_terms(row[hpo_idx])

    train_X, train_y = lst2array(train)

    train_preds = trainer(train_X, train_y)
    print(rank_list(train_preds))


if __name__ == "__main__":
    main()

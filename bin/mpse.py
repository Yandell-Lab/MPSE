#!/usr/bin/env python3

from os import mkdir
import os.path as path
import sys
import datetime as dt
import argparse
import random
import csv 
from processing.rady_process.rady_data_prep import ready, writey

from pyhpo.ontology import Ontology
from pyhpo.set import HPOSet

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from sklearn.model_selection import LeaveOneOut, cross_validate, cross_val_predict
from sklearn.naive_bayes import BernoulliNB
from sklearn import metrics


def argue():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--training", 
            default="processing/rady/rady_training_data.tsv", 
            help="Case/control training data in standard format.")
    parser.add_argument("-v", "--validate",
            default="processing/utah/utah_validation_data.tsv",
            help="Validation data in standard format.")
    parser.add_argument("-o", "--outdir", 
            default="analysis/example/",
            help="Output directory for results & reports.")
    return parser.parse_args()


def check_file_paths(args):
    if not path.isfile(args.training):
        print("Aborting - training input file not specified correctly.")
        sys.exit()

    if path.isdir(args.outdir):
        while True:
            decision = input("Output directory already exists.\nWould you like to proceed? (y/yes or n/no):\n")
            if decision == "y" or decision == "yes":
                break
            if decision == "n" or decision == "no":
                print("Aborting - please specify a different output directory.")
                sys.exit()
            else:
                print("Please provide valid response...\n")
                continue
    else:
        mkdir(args.outdir)
    try:
        mkdir(path.join(args.outdir, "tables"))
        mkdir(path.join(args.outdir, "figures"))
    except:
        return True


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
    y = df["seq_status"].astype("int8").to_numpy()
    return X,y


def trainer(X, y, mod=BernoulliNB(), cv=LeaveOneOut()):
    scores = cross_validate(mod, X, y,
            cv=cv,
            scoring=["accuracy",],
            return_train_score=True,
            n_jobs=-1)

    probas = cross_val_predict(mod, X, y, 
            cv=cv, 
            method="predict_proba", 
            n_jobs=-1)

    log_probas = cross_val_predict(mod, X, y, 
            cv=cv, 
            method="predict_log_proba", 
            n_jobs=-1)

    y_uniq = np.unique(y)
    indices = np.argmax(probas, axis=1)
    classes = np.expand_dims(y_uniq[indices], axis=1)
    return scores, np.hstack((probas, log_probas, classes))


def sample_cohort(data, diagnos_rate=0.18):
    cases_idx = get_col_position(data, "diagnostic")
    incident_idx = get_col_position(data, "incidental")

    cases = [x for x in data[1:] if x[cases_idx]=="1" and x[incident_idx]=="0"]
    controls = [x for x in data[1:] if x[cases_idx]!="1"]

    case_n = len(cases)
    control_n = len(controls)
    n = len(data)-1

    if case_n / n < diagnos_rate:
        control_n = round((1.0 - diagnos_rate) * case_n / diagnos_rate)
        control_samp = random.sample(controls, control_n)
        samp = [data[0]] + cases + control_samp
    else:
        case_n = round(diagnos_rate * control_n / (1.0 - diagnos_rate))
        case_samp = random.sample(cases, case_n)
        samp = [data[0]] + case_samp + controls
    return samp


def rocy(preds, outcome, fpath):
    fpr, tpr, thresholds = metrics.roc_curve(outcome, preds[:,1], pos_label=1)
    roc_auc = metrics.auc(fpr, tpr)
    roc_plot = metrics.RocCurveDisplay(fpr=fpr, tpr=tpr, roc_auc=roc_auc)
    roc_plot.plot()
    plt.savefig(fpath)
    return True


def main():
    args = argue()
    check_file_paths(args)
    _ = Ontology()

    train = ready(args.training, delim="\t")
    valid = ready(args.validate, delim="\t")

    hpo_idx = get_col_position(train, "hpo")
    for row in train[1:]:
        row[hpo_idx] = child_terms(row[hpo_idx])

    train_X, train_y = lst2array(train)

    train_scores, train_preds = trainer(train_X, train_y)
    preds_header = ["neg_proba","pos_proba","neg_log_proba","pos_log_proba","classes"]
    train = [x+y for x,y in zip(
        train, 
        [preds_header] + train_preds.tolist())] 
    train_sample = sample_cohort(train)

    rocy(train_preds, train_y, path.join(args.outdir, "figures/seq_status_ROC.png"))

    writey(train, path.join(args.outdir, "tables/training_predictions.csv"))
    writey(train_sample, path.join(args.outdir, "tables/training_predictions_sample.csv"))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

from os import mkdir
import os.path as path
import sys
import argparse
import random
import csv 

from pyhpo.ontology import Ontology
from pyhpo.set import HPOSet

from fhir.resources.observation import Observation

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from sklearn.model_selection import LeaveOneOut, cross_validate, cross_val_predict
from sklearn.naive_bayes import BernoulliNB
from sklearn import metrics


def argue():
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--training", 
            default="data/test/fake_training_data.tsv", 
            help="Case/control training data in standard format.")
    parser.add_argument("-v", "--validate",
            default="data/test/fake_validation_data.tsv",
            help="Validation data in standard format.")
    parser.add_argument("-F", "--FHIR", 
			action="store_true",
            help="Return results as FHIR Observation Resource (JSON).")
    parser.add_argument("-o", "--outdir", 
            default="analysis/test",
            help="Output directory for results & reports.")
    return parser.parse_args()


def ready(ftag, delim="\t", drop_header=False):
	with open(ftag, newline="") as f:
		reader = csv.reader(f, delimiter=delim)
		if drop_header:
			out = [row for row in reader][1:]
		else:
			out = [row for row in reader]
	return out 


def writey(data, ftag, header=None, delim="\t"):
	with open(ftag, "w", newline="") as f:
		writer = csv.writer(f, delimiter=delim)
		if header is not None:
			writer.writerow(header)
		writer.writerows(data)
	return True


def get_col_pos(data, col_names):
	idx_dic = {name: data[0].index(name) for name in col_names}
	return idx_dic


def child_terms(hpo):
    hpo_lst = hpo.split(";")
    hpo_set = HPOSet.from_queries(hpo_lst)
    hpo_subset = hpo_set.child_nodes()
    subset_dic = hpo_subset.toJSON()
    out_str = ";".join([x["id"] for x in subset_dic])
    return out_str


def compliant(data, dataset_name, col_idx):
	# check identifiers are unique
	ids = [x[0] for x in data[1:]]
	if len(ids) != len(set(ids)):
		sys.exit("{0} pids are not unique\nAborting...".format(dataset_name))

	# check value sets for seq_status, diagnostic, incidental
	# fill "" with 0
	check_cols = ["seq_status","diagnostic","incidental"]
	for row in data[1:]:
		for col in check_cols:
			if row[col_idx[col]] == "":
				row[cold_idx[col]] = "0"
			elif row[col_idx[col]] not in ["0","1"]:
				sys.exit("{0}: non-compliant {1} value\nAborting...".format(dataset_name, col))
	
	# order HPO list
	# call child_terms()
	for row in data[1:]:
		row[col_idx["hpo"]] = ";".join(sorted(row[col_idx["hpo"]].split(";")))
		row[col_idx["hpo"]] = child_terms(row[col_idx["hpo"]])

	return data


def lst2array(training, validation=None):
	if validation:
		train = pd.DataFrame(training[1:], columns=training[0])
		valid = pd.DataFrame(validation[1:], columns=validation[0])

		train_onehot = train["hpo"].str.get_dummies(sep=";")
		valid_onehot = valid["hpo"].str.get_dummies(sep=";")
		concat = pd.concat([train_onehot, valid_onehot], keys=["train","valid"])
		concat = concat.loc[:,concat.iloc[0,:].notna()].fillna(0, downcast="infer")

		train_X = concat.loc[["train"]].reset_index(drop=True).to_numpy()
		valid_X = concat.loc[["valid"]].reset_index(drop=True).to_numpy()
		train_y = train["seq_status"].astype("int8").to_numpy()
		valid_y = valid["seq_status"].astype("int8").to_numpy()
	else:
		df = pd.DataFrame(training[1:], columns=training[0])
		one_hot = df["hpo"].str.get_dummies(sep=";")
		train_X = one_hot.to_numpy()
		train_y = df["seq_status"].astype("int8").to_numpy()
		valid_X = None
		valid_y = None
	return (train_X, train_y, valid_X, valid_y)


def training(X, y, mod=BernoulliNB(), cv=LeaveOneOut()):
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


def validation(train_X, train_y, valid_X):
	mod = BernoulliNB()
	fit = mod.fit(train_X, train_y)

	probas = fit.predict_proba(valid_X)
	log_probas = fit.predict_log_proba(valid_X)

	y_uniq = np.unique(train_y)
	indices = np.argmax(probas, axis=1)
	classes = np.expand_dims(y_uniq[indices], axis=1)
	return np.hstack((probas, log_probas, classes))


def sample_cohort(data, col_idx, diagnos_rate=0.18):
    cases = [x for x in data[1:] if x[col_idx["diagnostic"]]=="1" and x[col_idx["incidental"]]=="0"]
    controls = [x for x in data[1:] if x[col_idx["diagnostic"]]!="1"]

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


def build_resources(data, col_idx):
	resources = []
	for pt in data:
		injest = {
				"resourceType": "Observation",
				"identifier": [{"value": pt[col_idx["pid"]]}],
				"status": "final",
				"code": {
					"coding": [{"system": "???", "code": "???", "display": "???"}], 
					"text": "???"
					}
				}
		obs = Observation.parse_obj(injest)
		resources.append(obs)
	return resources


def main():
	args = argue()
	_ = Ontology()
	
	train = ready(args.training)
	valid = ready(args.validate)
	
	data_col_names = ["pid","seq_status","diagnostic","incidental","hpo"]
	train_col_idx = get_col_pos(train, data_col_names)
	valid_col_idx = get_col_pos(valid, data_col_names)

	train = compliant(train, "train_data", train_col_idx)
	valid = compliant(valid, "valid_data", valid_col_idx)

	train_X, train_y, valid_X, valid_y = lst2array(train, valid)
	
	train_scores, train_preds = training(train_X, train_y)
	valid_preds = validation(train_X, train_y, valid_X)
	
	preds_header = ["neg_proba","pos_proba","neg_log_proba","pos_log_proba","classes"]
	train_out = [x+y for x,y in zip(
	    train, 
	    [preds_header] + train_preds.tolist())] 
	train_sample = sample_cohort(train_out, train_col_idx)

	valid_out = [x+y for x,y in zip(
		valid,
		[preds_header] + valid_preds.tolist())]
	
	if not path.isdir(args.outdir):
		mkdir(args.outdir)
		mkdir(path.join(args.outdir, "tables"))
		mkdir(path.join(args.outdir, "figures"))

	rocy(train_preds, train_y, path.join(args.outdir, "figures/seq_status_ROC.png"))
	
	writey(train_out, path.join(args.outdir, "tables/training_predictions.tsv"))
	writey(train_sample, path.join(args.outdir, "tables/training_predictions_sample.tsv"))
	writey(valid_out, path.join(args.outdir, "tables/validation_predictions.tsv"))

	if args.FHIR:
		resources = build_resources()
		for obs in resources:
			json.dumps(obs.json(indent=True))


if __name__ == "__main__":
    main()

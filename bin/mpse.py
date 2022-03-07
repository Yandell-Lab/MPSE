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
from scipy.stats import binom_test
from matplotlib import pyplot as plt

from sklearn.model_selection import LeaveOneOut, cross_validate, cross_val_predict
from sklearn.naive_bayes import BernoulliNB
from sklearn import metrics
from joblib import dump, load


def argue():
	parser = argparse.ArgumentParser()
	parser.add_argument("-t", "--training",
			default="data/test/fake_training_data.tsv", 
			help="Case/control training data in standard format.")
	parser.add_argument("-m", "--model",
			required=False,
			help="Serialized model (pickle object) to load from disc.")
	parser.add_argument("-p", "--prospective",
	        required=False,
	        help="Prospective data in standard format.")
	parser.add_argument("-s", "--select", 
			action="store_true",
	        help="Perform feature selection using binomial significance test.")
	parser.add_argument("-F", "--FHIR", 
			action="store_true",
	        help="Return results as FHIR Observation Resource (JSON).")
	parser.add_argument("-P", "--Pickle", 
			action="store_true",
	        help="Pickle model object to file 'data/test/model_obj.pickle'")
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


def compliant(data, dataset_name, col_idx, check_cols=None):
	# check identifiers are unique
	ids = [x[col_idx["pid"]] for x in data[1:]]
	if len(ids) != len(set(ids)):
		sys.exit("{0} pids are not unique\nAborting...".format(dataset_name))

	# check value sets for seq_status, diagnostic, incidental
	# fill "" with 0
	if check_cols:
		for row in data[1:]:
			for col in check_cols:
				if row[col_idx[col]] == "":
					row[col_idx[col]] = "0"
				elif row[col_idx[col]] not in ["0","1"]:
					sys.exit("{0}: non-compliant {1} value\nAborting...".format(dataset_name, col))
	
	# order HPO list
	# call child_terms()
	for row in data[1:]:
		row[col_idx["hpo"]] = ";".join(sorted(row[col_idx["hpo"]].split(";")))
		row[col_idx["hpo"]] = child_terms(row[col_idx["hpo"]])
	return data


def onehot_encode(data):
	df = pd.DataFrame(data[1:], columns=data[0])
	onehot = df["hpo"].str.get_dummies(sep=";")
	return onehot


def select_features(data, col_idx, alpha=0.05, direction_bias="two-sided"):
	onehot = onehot_encode(data)
	col_names = onehot.columns

	keep_idx = []
	drop_idx = []
	for i in range(onehot.shape[1]):
		x_b, x_t = 0,0
		n_b, n_t = 0,0
		for j in range(onehot.shape[0]):
			if data[1:][j][col_idx["seq_status"]] == "0":
				x_b += onehot.iloc[j,i]
				n_b += 1
			else:
				x_t += onehot.iloc[j,i]
				n_t += 1
		p = x_b / n_b
		test = binom_test(x=x_t, n=n_t, p=p, alternative=direction_bias)
		if test <= alpha:
			keep_idx.append(i)
		else:
			drop_idx.append(i)
	
	keep_terms = col_names[keep_idx]
	drop_terms = col_names[drop_idx]
	trimmed = onehot.iloc[:,keep_idx]
	return trimmed, keep_terms, drop_terms


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
    scrs = np.log(probas[:,1] / probas[:,0])
    return scores, np.hstack((probas, log_probas, classes, scrs[:, np.newaxis]))


def score_probands(mod, valid_X):
	probas = mod.predict_proba(valid_X)
	log_probas = mod.predict_log_proba(valid_X)

	y_uniq = np.unique(mod.classes_)
	indices = np.argmax(probas, axis=1)
	classes = np.expand_dims(y_uniq[indices], axis=1)
	scrs = np.log(probas[:,1] / probas[:,0])
	return np.hstack((probas, log_probas, classes, scrs[:, np.newaxis]))


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

	if not path.isdir(args.outdir):
		mkdir(args.outdir)
		mkdir(path.join(args.outdir, "tables"))
	
	if args.model:
		mod = load(args.model)
		prosp = ready(args.prospective)
		prosp_col_idx = get_col_pos(prosp, ["pid","hpo"])
		prosp = compliant(prosp, "prosp_data", prosp_col_idx)

		keep_terms = mod.feature_names_in_
		df_concat = [onehot_encode(prosp), pd.DataFrame(columns=keep_terms)]
		prosp_X = pd.concat(df_concat)[keep_terms].fillna(0).astype("int8")

		prosp_preds = score_probands(mod, prosp_X)
		preds_header = ["neg_proba","pos_proba","neg_log_proba","pos_log_proba","class","scr"]
		prosp_out = [x+y for x,y in zip(prosp, [preds_header] + prosp_preds.tolist())]
		writey(prosp_out, path.join(args.outdir, "tables/prospective_predictions.tsv"))
	else:
		train = ready(args.training)
		col_pos_names = ["pid","seq_status","diagnostic","incidental","hpo"]
		train_col_idx = get_col_pos(train, col_pos_names)
		train = compliant(train, 
				"train_data", 
				train_col_idx, 
				check_cols=["seq_status","diagnostic","incidental"])

		if args.select:
			train_X, keep_terms, drop_terms = select_features(train, train_col_idx)
		else:
			train_X = onehot_encode(train)
			keep_terms = train_X.columns
		train_y = np.array([x[train_col_idx["seq_status"]] for x in train[1:]])
		if args.Pickle:
			fit = BernoulliNB().fit(train_X, train_y)
			dump(fit, "data/test/model_obj.pickle")

		train_scores, train_preds = training(train_X, train_y)
		preds_header = ["neg_proba","pos_proba","neg_log_proba","pos_log_proba","class","scr"]
		train_out = [x+y for x,y in zip(train, [preds_header] + train_preds.tolist())] 
		train_sample = sample_cohort(train_out, train_col_idx)

		writey(train_out, path.join(args.outdir, "tables/training_predictions.tsv"))
		writey(train_sample, path.join(args.outdir, "tables/training_predictions_sample.tsv"))

		if args.prospective:
			prosp = ready(args.prospective)
			prosp_col_idx = get_col_pos(prosp, ["pid","hpo"])
			prosp = compliant(prosp, "prosp_data", prosp_col_idx)

			df_concat = [onehot_encode(prosp), pd.DataFrame(columns=keep_terms)]
			prosp_X = pd.concat(df_concat)[keep_terms].fillna(0).astype("int8")

			mod = BernoulliNB().fit(train_X, train_y)
			prosp_preds = score_probands(mod, prosp_X)
			prosp_out = [x+y for x,y in zip(prosp, [preds_header] + prosp_preds.tolist())]
			writey(prosp_out, path.join(args.outdir, "tables/prospective_predictions.tsv"))

	if args.FHIR:
		resources = build_resources()
		for obs in resources:
			json.dumps(obs.json(indent=True))


if __name__ == "__main__":
    main()

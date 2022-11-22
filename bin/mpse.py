#!/usr/bin/env python3

from os import mkdir
import os.path as path
import sys
import argparse
import random
import json
import csv 
import re

from pyhpo.ontology import Ontology
from pyhpo.set import HPOSet
import simple_icd_10_cm as cm
#from fhir.resources.observation import Observation

from datetime import date
from datetime import datetime as dt
from datetime import timezone as tz

from collections import defaultdict

import numpy as np
import pandas as pd
from scipy.stats import binom_test
from matplotlib import pyplot as plt

from sklearn.model_selection import LeaveOneOut, cross_validate, cross_val_predict
from sklearn.naive_bayes import BernoulliNB
from sklearn import metrics
from joblib import dump, load

csv.field_size_limit(sys.maxsize)


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
	parser.add_argument("--timestamps", 
			action="store_true",
			help="Input features have associated timestamps.")
	parser.add_argument("-a", "--alpha", 
			type=float,
			default=1.0,
			help="""Binomial test alpha used for feature selection. 
			Default: no feature selection.""")
	parser.add_argument("--sample_features", 
			type=float,
			default=1.0,
			help="""Fraction of features to sample 
			from training data for resilience testing. 
			Default: no sampling.""")
	parser.add_argument("--fudge_terms", 
			type=int,
			default=0,
			help="""Number of random terms to add/subtract 
			from prospective observations for resilience testing. 
			Default: no fudging.""")
	parser.add_argument("-R", "--Rady",
			action="store_true",
			help="Prospective data comes from Rady Clinithink process.")
	parser.add_argument("-C", "--Cardinal",
			action="store_true",
			help="Return cardinal phenotypes for prospective cases.")
	parser.add_argument("-P", "--Pickle", 
			action="store_true",
	        help="Dump pickled model object to file '{outdir}/trained_model.pickle'")
#	parser.add_argument("-F", "--FHIR", 
#			action="store_true",
#	        help="Return results as FHIR Observation Resource (JSON).")
	parser.add_argument("--json",
			action="store_true",
	        help="Return results as JSON object.")
	parser.add_argument("-o", "--outdir", 
	        default="analysis/test",
	        help="Output directory for results & reports.")
	return parser


def check_args(parser):
	args = parser.parse_args()
	low, high = 0.0, 1.0
	if args.alpha <= low or args.alpha > high:
		print("Arg 'ALPHA' defined on interval ({0},{1}]\nAborting".format(low, high))
		parser.print_usage()
		sys.exit()

	low, high = 0.5, 1.0
	if args.sample_features < low or args.sample_features > high:
		print("Arg 'SAMPLE_FEATURES' defined on interval [{0},{1}]\nAborting".format(low, high))
		parser.print_usage()
		sys.exit()

	low, high = -50, 50
	if args.fudge_terms < low or args.fudge_terms > high:
		print("Arg 'FUDGE_TERMS' defined on interval [{0},{1}]\nAborting".format(low, high))
		parser.print_usage()
		sys.exit()
	return True


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


def hpo_parse(hpo_str):
	hpo_reg = re.compile(r"hp\d{7}")
	srch = hpo_reg.search(hpo_str).group()
	return srch.replace("hp", "HP:", 1)


def extract_timestamps(data, col_idx):
	extract = [data[0] + ["manifest_date"]]
	for row in data[1:]:
		tmsp = row[col_idx["codes"]].split(";")
		tmsp = [(date.fromisoformat(x.split("|")[1]), x.split("|")[0]) for x in tmsp]

		d = defaultdict(set)
		for k, v in tmsp:
			d[k].add(v)

		dsort = sorted(d.items())
		previous = set()
		for k, v in dsort:
			new_row = row + [k.isoformat()]
			new_cde = v | previous
			new_row[col_idx["codes"]] = ";".join(list(new_cde))
			extract.append(new_row)
			previous = new_cde
	return extract


def child_terms(hpo_lst):
    hpo_set = HPOSet.from_queries(hpo_lst)
    hpo_subset = hpo_set.child_nodes()
    subset_dic = hpo_subset.toJSON()
    out_str = sorted([x["id"] for x in subset_dic])
    return out_str


def clean_codes(codes, keep_others):
	"""
	Basic pre-processing of input codes in preparation for modeling.
	Code types that are currently handled...
	- Human Phenotype Ontology (HPO)
	- ICD-10-CM
	Other code types are optionally included without processing.
	"""
	valid_hpo = Ontology.to_dataframe().index.tolist()
	hpo = []
	icd= []
	other = []
	for cde in codes:
		if cde in valid_hpo:
			hpo.append(cde)
		elif cm.is_valid_item(cde):
			icd.append(cm.add_dot(cde))
		else:
			other.append(cde)
	hpo_clean = child_terms(hpo)
	icd_clean = sorted(icd)
	other_clean = sorted(other) if keep_others else []
	return hpo_clean + icd_clean + other_clean


def compliant(data, dataset_name, col_idx, check_cols=None, keep_all_codes=False):
	# check identifiers are unique
	ids = [x[col_idx["pid"]] for x in data[1:]]
	if len(ids) != len(set(ids)):
		msg = "Warning: the dataset '{0}' PIDs are not unique. Please check this is expected..."
		print(msg.format(dataset_name), file=sys.stderr)

	# check value sets for seq_status, diagnostic, incidental
	# fill "" with 0
	if check_cols:
		for row in data[1:]:
			for col in check_cols:
				if row[col_idx[col]] == "":
					row[col_idx[col]] = "0"
				elif row[col_idx[col]] not in ["0","1"]:
					msg = "{0}: non-compliant value for column '{1}'\nAborting..."
					sys.exit(msg.format(dataset_name, col))
	
	# order code list and remove duplicate codes
	# clean codes
	data[0].append("codes_clean")
	for row in data[1:]:
		dirty = sorted(set(row[col_idx["codes"]].split(";")))
		clean = clean_codes(dirty, keep_others=keep_all_codes)
		row[col_idx["codes"]] = ";".join(dirty)
		row.append(";".join(clean))
	return data


def onehot_encode(data):
	df = pd.DataFrame(data[1:], columns=data[0])
	onehot = df["codes_clean"].str.get_dummies(sep=";")
	return onehot


def select_features(data, col_idx, alpha, direction_bias="two-sided"):
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


def sample_features(data, frac):
	onehot = onehot_encode(data)
	ncol = onehot.shape[1]
	sample_n = int(ncol * frac)
	keep_idx = np.random.choice(range(ncol), size=sample_n, replace=False)
	sort_idx = np.sort(keep_idx)
	sampled = onehot.iloc[:,sort_idx]
	return sampled, sampled.columns, sample_n


def fudge_terms(data, col_idx, sample_set, fudge_n):
	if fudge_n > 0:
		for row in data[1:]:
			choose = np.random.choice(sample_set, size=fudge_n, replace=False)
			add = ";".join(choose)
			new = row[col_idx["codes_clean"]] + ";" + add
			row[col_idx["codes_clean"]] = ";".join(sorted(new.split(";")))
	else:
		for row in data[1:]:
			new = row[col_idx["codes_clean"]].split(";")
			choose = np.random.choice(range(len(new)), size=abs(fudge_n), replace=False)
			for i in np.sort(choose)[::-1]:
				del new[i]
			row[col_idx["codes_clean"]] = ";".join(new)
	return data


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

	y_uniq = mod.classes_
	indices = np.argmax(probas, axis=1)
	classes = np.expand_dims(y_uniq[indices], axis=1)
	scrs = np.log(probas[:,1] / probas[:,0])
	return np.hstack((probas, log_probas, classes, scrs[:, np.newaxis]))


def rocy(preds, outcome):
	fpr, tpr, thresholds = metrics.roc_curve(outcome.astype("int8")[:, np.newaxis], preds[:,1], pos_label=1)
	roc_auc = metrics.auc(fpr, tpr)
	print(roc_auc, file=sys.stderr)
	return True


def get_cardinal(pid, data, feature_probs):
	term_names = data.columns
	valid_hpo = Ontology.to_dataframe().index.tolist()
	cards = []
	for idx, row in data.iterrows():#.tolist():
		boolpt = list(map(bool, row))
		card = (term_names[boolpt].tolist(), feature_probs[boolpt].tolist())
		for feat in range(len(card[0])):
			if card[0][feat] in valid_hpo:
				cards.append([
					pid[idx], 
					"HPO",
					card[0][feat], 
					Ontology.get_hpo_object(card[0][feat]).name, 
					card[1][feat]])
			elif cm.is_valid_item(card[0][feat]):
				cards.append([
					pid[idx],
					"ICD-10-CM",
					card[0][feat],
					cm.get_description(card[0][feat]),
					card[1][feat]])
			else:
				cards.append([
					pid[idx],
					"Unspecified",
					card[0][feat],
					"Unknown",
					card[1][feat]])
	return cards


def null_dist(fit, data, col_idx, keep_terms, preds_header):
	# Create a 'null' distribution of randomly generated code lists
	term_cnts = [len(x[col_idx["codes_clean"]].split(";")) for x in data[1:]]
	null_samp = [np.random.choice(keep_terms, size=n, replace=False) for n in term_cnts]
	null_samp = [["codes_clean"]] + [[";".join(x)] for x in null_samp]
	df_concat = [onehot_encode(null_samp), pd.DataFrame(columns=keep_terms)]
	null_X = pd.concat(df_concat)[keep_terms].fillna(0).astype("int8")
	null_preds = score_probands(fit, null_X)
	null_out = [x+y for x,y in zip(null_samp, [preds_header] + null_preds.tolist())]
	#writey(null_out, path.join(args.outdir, "null_preds.tsv"))
	return null_out


#def build_resources(data, col_idx):
#	stamp = dt.now() #timezone??
#	resources = []
#	for pt in data[1:]:
#		injest = {
#				"resourceType": "Observation",
#				"identifier": [{"value": pt[col_idx["pid"]]}],
#				"status": "final",
#				"code": {
#					#"coding": [{"system": "???", "code": "???", "display": "???"}], 
#					"text": "MPSE score: {0}".format(pt[col_idx["scr"]])
#					},
#				"effectiveDateTime": stamp.isoformat()
#				}
#		obs = Observation.parse_obj(injest)
#		resources.append(obs)
#		#resources.append(injest)
#	return resources


def build_JSON(data, col_idx, cards):
	stamp = dt.now()

	card_dict = {}
	for row in cards[1:]:
		card = {"term_id": row[1], "term_name": row[2], "coef": row[3]}
		if row[0] in card_dict:
			card_dict[row[0]].append(card)
		else:
			card_dict[row[0]] = [card]

	pt_lst = []
	for row in data[1:]:
		pt_id = row[col_idx["pid"]]
		pt = {
				"identifier": pt_id,
				"mpse_score": row[col_idx["scr"]],
				"cardinal_phenotypes": card_dict[pt_id]
				}
		pt_lst.append(pt)

	JSON = {
			"timestamp": stamp.isoformat(),
			"MPSE_version": "1.0",
			"MPSE_manifest": pt_lst
			}
	return JSON


def main():
	parser = argue()
	args = parser.parse_args()
	check_args(parser)
	_ = Ontology()

	if not path.isdir(args.outdir):
		mkdir(args.outdir)
	
	if args.model:
		mod = load(args.model)
		keep_terms = mod.feature_names_in_

		if args.Rady:
			raw = ready(args.prospective, delim=",", drop_header=True)
			cde_lst = [hpo_parse(x[0]) for x in raw]
			prosp = [["pid","codes"], [path.basename(args.prospective), ";".join(cde_lst)]]
		else:
			prosp = ready(args.prospective)

		prosp_col_idx = get_col_pos(prosp, ["pid","codes"])

		if args.timestamps:
			prosp = extract_timestamps(prosp, prosp_col_idx)

		if args.fudge_terms != 0:
			prosp = fudge_terms(prosp, prosp_col_idx, keep_terms, args.fudge_terms)
		prosp = compliant(prosp, "prosp_data", prosp_col_idx, False)

		df_concat = [onehot_encode(prosp), pd.DataFrame(columns=keep_terms)]
		prosp_X = pd.concat(df_concat)[keep_terms].fillna(0).astype("int8")

		prosp_preds = score_probands(mod, prosp_X)
		preds_header = ["neg_proba","pos_proba","neg_log_proba","pos_log_proba","class","scr"]
		prosp_out = [x+y for x,y in zip(prosp, [preds_header] + prosp_preds.tolist())]

		prosp_writer = csv.writer(sys.stdout, delimiter="\t")
		prosp_writer.writerows(prosp_out)

		if args.Cardinal:
			coefs = mod.feature_log_prob_
			prosp_pid = [x[0] for x in prosp[1:]]
			cards = get_cardinal(prosp_pid, prosp_X, coefs[1] - coefs[0])
			writey(cards, 
					path.join(args.outdir, "cardinal_phenotypes.tsv"), 
					header=["pid","domain","term_id","term_name","coef"])

	else:
		train = ready(args.training)
		col_pos_names = ["pid","seq_status","diagnostic","codes"]
		train_col_idx = get_col_pos(train, col_pos_names)
		train = compliant(train, 
				"train_data", 
				train_col_idx, 
				["seq_status","diagnostic"],
				False)

		if args.alpha != 1.0:
			train_X, keep_terms, drop_terms = select_features(train, train_col_idx, args.alpha)
		elif args.sample_features != 1.0:
			train_X, keep_terms, sample_n = sample_features(train, args.sample_features)
		else:
			train_X = onehot_encode(train)
			keep_terms = train_X.columns
		train_y = np.array([x[train_col_idx["seq_status"]] for x in train[1:]])

		train_scores, train_preds = training(train_X, train_y)
		#rocy(train_preds, train_y)
		preds_header = ["neg_proba","pos_proba","neg_log_proba","pos_log_proba","class","scr"]
		train_out = [x+y for x,y in zip(train, [preds_header] + train_preds.tolist())] 

		writey(train_out, path.join(args.outdir, 
			"training_preds_ba{0}_sf{1}.tsv".format(args.alpha, args.sample_features)))

		fit = BernoulliNB().fit(train_X, train_y)

		if args.Pickle:
			dump(fit, path.join(args.outdir, "trained_model.pickle"))

		if args.prospective:
			if args.Rady:
				raw = ready(args.prospective, delim=",", drop_header=True)
				cde_lst = [hpo_parse(x[0]) for x in raw]
				prosp = [["pid","codes"], [path.basename(args.prospective), ";".join(cde_lst)]]
			else:
				prosp = ready(args.prospective)

			prosp_col_idx = get_col_pos(prosp, ["pid","codes"])

			if args.timestamps:
				prosp = extract_timestamps(prosp, prosp_col_idx)

			if args.fudge_terms != 0:
				prosp = fudge_terms(prosp, prosp_col_idx, keep_terms, args.fudge_terms)
			prosp = compliant(prosp, "prosp_data", prosp_col_idx, False)

			df_concat = [onehot_encode(prosp), pd.DataFrame(columns=keep_terms)]
			prosp_X = pd.concat(df_concat)[keep_terms].fillna(0).astype("int8")

			prosp_preds = score_probands(fit, prosp_X)
			prosp_out = [x+y for x,y in zip(prosp, [preds_header] + prosp_preds.tolist())]

			prosp_writer = csv.writer(sys.stdout, delimiter="\t")
			prosp_writer.writerows(prosp_out)

			if args.Cardinal:
				coefs = fit.feature_log_prob_
				prosp_pid = [x[0] for x in prosp[1:]]
				cards = get_cardinal(prosp_pid, prosp_X, coefs[1] - coefs[0])
				writey(cards, 
						path.join(args.outdir, "cardinal_phenotypes.tsv"), 
						header=["pid","domain","term_id","term_name","coef"])

#	if args.FHIR:
#		resources = build_resources(prosp_out, get_col_pos(prosp_out, ["pid","scr"]))
#		for obs in resources:
#			fname = path.join(args.outdir, "{0}_FHIR.json".format(obs["identifier"][0]["value"]))
#			with open(fname, "w") as json_out:
#				json.dump(obs, json_out)

	if args.prospective and args.json:
		JSON = build_JSON(prosp_out, get_col_pos(prosp_out, ["pid","scr"]), cards)
		fname = path.join(args.outdir, "prospective_preds.json")
		with open(fname, "w") as json_out:
			json.dump(JSON, json_out, indent=4)
	


if __name__ == "__main__":
    main()

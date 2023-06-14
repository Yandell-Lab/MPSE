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
from joblib import dump, load

import numpy as np
import pandas as pd
from scipy.stats import binom_test
from matplotlib import pyplot as plt

from sklearn.model_selection import KFold, StratifiedKFold, LeaveOneOut
from sklearn.model_selection import cross_validate, cross_val_predict
from sklearn import metrics

from sklearn.naive_bayes import BernoulliNB
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC

csv.field_size_limit(sys.maxsize)


def argue():
    """Parses command line arguments.

    Returns:
        ArgumentParser: An ArgumentParser object that can be used to define and parse command line arguments.
    """
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
    parser.add_argument("--compare_models", 
            action="store_true",
            help="Train & test input data using multiple classifiers.")
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
#   parser.add_argument("-F", "--FHIR", 
#           action="store_true",
#           help="Return results as FHIR Observation Resource (JSON).")
    parser.add_argument("--json",
            action="store_true",
            help="Return results as JSON object.")
    parser.add_argument("-o", "--outdir", 
            default="analysis/test",
            help="Output directory for results & reports.")
    return parser


def check_args(parser):
    """Performs validity checks for the testing parameters 'alpha', 'sample_features', and 'fudge_terms'.

    Args:
        parser (ArgumentParser): An ArgumentParser object that defines and parses command line arguments.

    Returns:
        bool: True if all arguments pass the validity checks.
    """
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
    """Reads a delimited file and returns its contents as a list of rows.

    Args:
        ftag (str): The path or filename of the delimited file to be read.
        delim (str, optional): The delimiter used in the delimited file. Defaults to "\t".
        drop_header (bool, optional): Whether to drop the header row from the delimited file. Defaults to False.

    Returns:
        list: A list containing the rows of the delimited file, where each row is a list of its fields.
    """
    with open(ftag, newline="") as f:
        reader = csv.reader(f, delimiter=delim)
        if drop_header:
            out = [row for row in reader][1:]
        else:
            out = [row for row in reader]
    return out 


def writey(data, ftag, header=None, delim="\t"):
    """Writes data to a delimited file.

    Args:
        data (iterable): The data to be written to the delimited file.
        ftag (str): The path or filename of the delimited file to be written.
        header (iterable or None, optional): The header row to be written to the delimited file. Defaults to None.
        delim (str, optional): The delimiter to be used in the delimited file. Defaults to "\t".

    Returns:
        bool: True if the data is successfully written to the delimited file.
    """
    with open(ftag, "w", newline="") as f:
        writer = csv.writer(f, delimiter=delim)
        if header is not None:
            writer.writerow(header)
        writer.writerows(data)
    return True


def get_column_positions(data, col_names):
    """Returns a dictionary mapping column names to their corresponding indices in the data.

    Args:
        data (list): A list of lists representing tabular data.
        col_names (list): A list of column names.

    Returns:
        dict: A dictionary mapping column names to their corresponding indices in the data.
    """
    idx_dic = {name: data[0].index(name) for name in col_names}
    return idx_dic


def parse_hpo(hpo_str):
    """Parses and formats an HPO code from a string.

    Args:
        hpo_str (str): A string containing an HPO code with variable formatting.

    Returns:
        str: The parsed and formatted HPO code in the format 'HP:#######'.
    """
    pattern = r'(?i)([hp]{2}):?(\d{7})'
    replacement = r'HP:\2'
    return re.sub(pattern, replacement, hpo_str)


def extract_timestamps(data, col_idx):
    """Extracts timestamp-code pairs ('HP:#######|yyyy-mm-dd') and returns a list for each unique timestamp that includes all codes from that and previous timestamps.

    Args:
        data (list): A list of lists representing tabular data.
        col_idx (dict): A dictionary mapping column names to their corresponding indices.

    Returns:
        list: A modified extract of the data with added "manifest_date" column.
    """
    extract = [data[0] + ["manifest_date"]]
    for row in data[1:]:
        if row[col_idx["codes"]] != "":
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
        else:
            extract.append(row + [""])
    return extract


def remove_parent_terms(hpo_lst):
    """Retrieves the most specific HPO term of each subtree.

    Args:
        hpo_lst (list): A list of HPO terms.

    Returns:
        list: A sorted list of child terms' IDs.
    """
    hpo_set = HPOSet.from_queries(hpo_lst)
    hpo_subset = hpo_set.child_nodes()
    subset_dic = hpo_subset.toJSON()
    out_str = sorted([x["id"] for x in subset_dic])
    return out_str


def clean_codes(codes, keep_others=False):
    """Basic pre-processing of input codes in preparation for modeling.

    Args:
        codes (iterable): An iterable of codes.
        keep_others (bool, optional): Whether to keep non-HPO/ICD-10-CM codes. Defaults to False.

    Returns:
        list: A filtered list of codes.
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
    hpo_clean = remove_parent_terms(hpo)
    icd_clean = sorted(icd)
    other_clean = sorted(other) if keep_others else []
    return hpo_clean + icd_clean + other_clean


def make_compliant(data, dataset_name, col_idx, check_cols=None, keep_all_codes=False):
    """Performs data compliance checks and modifications on the provided dataset.

    Args:
        data (list): A list of lists representing tabular data.
        dataset_name (str): The name of the dataset.
        col_idx (dict): A dictionary mapping column names to their corresponding indices.
        check_cols (list, optional): A list of column names to perform value set checks on. Defaults to None.
        keep_all_codes (bool, optional): Whether to keep all codes or remove unrecognized codes. Defaults to False.

    Returns:
        list: The modified dataset after performing compliance checks and modifications.
    """
    # remove rows with no codes
    data = [row for row in data if row[col_idx["codes"]] != ""]

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
                    msg = "{0}: non-compliant value for column '{1}'; must be '0' or '1'"
                    raise ValueError(msg.format(dataset_name, col))
                    #sys.exit(msg.format(dataset_name, col))
    
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
    """Performs one-hot encoding on the 'codes_clean' column of the provided data.

    Args:
        data (list): A list of lists representing tabular data.

    Returns:
        pandas.DataFrame: The one-hot encoded representation of the 'codes_clean' column.
    """
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
    """Performs training and evaluation of a classifier model on the provided data.

    Args:
        X (array-like): The feature matrix.
        y (array-like): The target variable.
        mod (object, optional): The classifier model to use. Defaults to BernoulliNB().
        cv (object, optional): The cross-validation strategy to use. Defaults to LeaveOneOut().

    Returns:
        tuple: A tuple containing the evaluation scores and the concatenated results.
    """
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


def run_multimodels(X_train, y_train, X_test, y_test=None):
    """Runs multiple classification models and evaluates their performance using cross-validation.

    Args:
        X_train (array-like): Training data features.
        y_train (array-like): Training data labels.
        X_test (array-like): Test data features.
        y_test (array-like): Test data labels.

    Returns:
        pandas.DataFrame: Concatenated results of the cross-validation evaluation for each model.
    """
    dfs = []
    models = [ 
        ('LogReg', LogisticRegression()),
        ('DT', DecisionTreeClassifier()),
        ('RF', RandomForestClassifier()),
        ('KNN', KNeighborsClassifier()),
        ('SVM', SVC()),
        ('BNB', BernoulliNB()),
        ('GBM', GradientBoostingClassifier())
    ]
    results = []
    names = []
    #scoring = ['accuracy', 'precision_weighted', 'recall_weighted', 'f1_weighted', 'roc_auc']
    #scoring = ['accuracy', 'precision', 'recall', 'f1', 'roc_auc']
    scoring = {"accuracy": "accuracy",
            "precision": metrics.make_scorer(metrics.precision_score, pos_label="1"),
            "recall": metrics.make_scorer(metrics.recall_score, pos_label="1"),
            "f1": metrics.make_scorer(metrics.f1_score, pos_label="1"),
            "auc": "roc_auc"
            }
    target_names = ["0", "1"]
    for name, model in models:
        #kfold = KFold(n_splits=5, shuffle=True, random_state=1234)
        skfold = StratifiedKFold(n_splits=3, shuffle=True, random_state=1234)
        #loocv = LeaveOneOut()
        cv_results = cross_validate(model, X_train, y_train, cv=skfold, scoring=scoring)
        clf = model.fit(X_train, y_train)
        y_pred = clf.predict(X_test)
        print(name)
        print(metrics.classification_report(y_test, y_pred, target_names=target_names))
        results.append(cv_results)
        names.append(name)
        this_df = pd.DataFrame(cv_results)
        this_df["model"] = name
        dfs.append(this_df)
    final = pd.concat(dfs, ignore_index=True)
    return final


def score_probands(mod, valid_X):
    """Scores probands using the trained model and returns the predicted probabilities and scores.

    Args:
        mod (object): The trained classifier model.
        valid_X (array-like): The feature matrix of the probands to be scored.

    Returns:
        array-like: The concatenated results including predicted probabilities, log probabilities, predicted classes, and scores.
    """
    probas = mod.predict_proba(valid_X)
    log_probas = mod.predict_log_proba(valid_X)

    y_uniq = mod.classes_
    indices = np.argmax(probas, axis=1)
    classes = np.expand_dims(y_uniq[indices], axis=1)
    scrs = np.log(probas[:,1] / probas[:,0])
    return np.hstack((probas, log_probas, classes, scrs[:, np.newaxis]))


def process_prospective(mod, keep_terms, header, args):
    """Processes and scores prospective data using the provided model.

    Args:
        mod (object): The model for scoring the prospective data.
        keep_terms (list): A list of terms to keep in the processed data.
        header (list): The header row for the output data.
        args (ArgumentParser): The command line arguments.

    Returns:
        tuple: A tuple containing the processed prospective data, the processed prospective features, and the output data.
    """
    if args.Rady:
        raw = ready(args.prospective, delim=",", drop_header=True)
        cde_lst = [parse_hpo(x[0]) for x in raw]
        prosp = [["pid","codes"], [path.basename(args.prospective), ";".join(cde_lst)]]
    else:
        prosp = ready(args.prospective)

    prosp_col_idx = get_column_positions(prosp, ["pid","codes"])

    if args.timestamps:
        prosp = extract_timestamps(prosp, prosp_col_idx)

    if args.fudge_terms != 0:
        prosp = fudge_terms(prosp, prosp_col_idx, keep_terms, args.fudge_terms)
    prosp = make_compliant(prosp, "prosp_data", prosp_col_idx, False)

    df_concat = [onehot_encode(prosp), pd.DataFrame(columns=keep_terms)]
    prosp_X = pd.concat(df_concat)[keep_terms].fillna(0).astype("int8")

    prosp_preds = score_probands(mod, prosp_X)
    prosp_out = [x+y for x,y in zip(prosp, [header] + prosp_preds.tolist())]
    return prosp, prosp_X, prosp_out


def rocy(preds, outcome):
    """Computes the Receiver Operating Characteristic (ROC) curve and the Area Under the Curve (AUC).

    Args:
        preds (array-like): The predicted probabilities of the positive class.
        outcome (array-like): The true outcome labels.

    Returns:
        bool: True if the ROC AUC is successfully computed.
    """
    fpr, tpr, thresholds = metrics.roc_curve(outcome.astype("int8")[:, np.newaxis], preds[:,1], pos_label=1)
    roc_auc = metrics.auc(fpr, tpr)
    print(roc_auc, file=sys.stderr)
    return True


def get_cardinal(pid, data, feature_probs):
    """Retrieves the cardinality information for the given patient IDs from the data and feature probabilities.

    Args:
        pid (list): List of patient IDs.
        data (pandas.DataFrame): The data containing the features.
        feature_probs (list): The feature probabilities (aka coefficients).

    Returns:
        list: A list of lists representing the cardinality information for each feature.
              Each sublist contains the patient ID, feature type, feature code, feature name, and feature probability.
    """
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
#   stamp = dt.now() #timezone??
#   resources = []
#   for pt in data[1:]:
#       injest = {
#               "resourceType": "Observation",
#               "identifier": [{"value": pt[col_idx["pid"]]}],
#               "status": "final",
#               "code": {
#                   #"coding": [{"system": "???", "code": "???", "display": "???"}], 
#                   "text": "MPSE score: {0}".format(pt[col_idx["scr"]])
#                   },
#               "effectiveDateTime": stamp.isoformat()
#               }
#       obs = Observation.parse_obj(injest)
#       resources.append(obs)
#       #resources.append(injest)
#   return resources


def build_JSON(data, col_idx, cards):
    """Builds a JSON representation of the data and cardinality information.

    Args:
        data (list): A list of lists representing tabular data.
        col_idx (dict): A dictionary mapping column names to their corresponding indices.
        cards (list): A list of lists representing the cardinality information.

    Returns:
        dict: A JSON object containing the timestamp, MPSE version, and MPSE manifest.
              The MPSE manifest includes patient information, MPSE score, and cardinal phenotypes.
    """
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
    preds_header = ["neg_proba","pos_proba","neg_log_proba","pos_log_proba","class","scr"]

    if not path.isdir(args.outdir):
        mkdir(args.outdir)
    
    if args.model:
        mod = load(args.model)
        keep_terms = mod.feature_names_in_

        prosp, prosp_X, prosp_out = process_prospective(mod, keep_terms, preds_header, args)
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
        train_col_idx = get_column_positions(train, col_pos_names)
        train = make_compliant(train, 
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
        train_out = [x+y for x,y in zip(train, [preds_header] + train_preds.tolist())] 

        writey(train_out, path.join(args.outdir, 
            "training_preds_ba{0}_sf{1}.tsv".format(args.alpha, args.sample_features)))

        fit = BernoulliNB().fit(train_X, train_y)

        if args.Pickle:
            dump(fit, path.join(args.outdir, "trained_model.pickle"))

        if args.prospective:
            if args.compare_models:
                prosp = ready(args.prospective)
                prosp_col_idx = get_column_positions(prosp, ["pid","codes","seq_status"])
                prosp = make_compliant(prosp, "prosp_data", prosp_col_idx, False)
                df_concat = [onehot_encode(prosp), pd.DataFrame(columns=keep_terms)]
                prosp_X = pd.concat(df_concat)[keep_terms].fillna(0).astype("int8")
                prosp_y = np.array([x[prosp_col_idx["seq_status"]] for x in prosp[1:]])
                model_comparison = run_multimodels(train_X, train_y, prosp_X, prosp_y)
                print(model_comparison)
                sys.exit()
            else:
                prosp, prosp_X, prosp_out = process_prospective(fit, keep_terms, preds_header, args)
                prosp_writer = csv.writer(sys.stdout, delimiter="\t")
                prosp_writer.writerows(prosp_out)

                if args.Cardinal:
                    coefs = fit.feature_log_prob_
                    prosp_pid = [x[0] for x in prosp[1:]]
                    cards = get_cardinal(prosp_pid, prosp_X, coefs[1] - coefs[0])
                    writey(cards, 
                            path.join(args.outdir, "cardinal_phenotypes.tsv"), 
                            header=["pid","domain","term_id","term_name","coef"])

#   if args.FHIR:
#       resources = build_resources(prosp_out, get_column_positions(prosp_out, ["pid","scr"]))
#       for obs in resources:
#           fname = path.join(args.outdir, "{0}_FHIR.json".format(obs["identifier"][0]["value"]))
#           with open(fname, "w") as json_out:
#               json.dump(obs, json_out)

    if args.prospective and args.json:
        JSON = build_JSON(prosp_out, get_column_positions(prosp_out, ["pid","scr"]), cards)
        fname = path.join(args.outdir, "prospective_preds.json")
        with open(fname, "w") as json_out:
            json.dump(JSON, json_out, indent=4)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3

import re
import sys
import csv 
import json
import argparse
from os import mkdir
import os.path as path

import simple_icd_10_cm as cm
from pyhpo.set import BasicHPOSet
from pyhpo.ontology import Ontology

from datetime import date
from datetime import datetime as dt

from joblib import dump, load
from collections import defaultdict

import numpy as np
import pandas as pd

from sklearn import metrics
from sklearn.model_selection import cross_validate, cross_val_predict
from sklearn.model_selection import KFold, StratifiedKFold, LeaveOneOut

from sklearn.svm import SVC
from sklearn.naive_bayes import BernoulliNB
from sklearn.tree import DecisionTreeClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier

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
    parser.add_argument("--keep_all_codes",
            action="store_true",
            help="Keep non-vocabulary codes.")
    parser.add_argument("--timestamps", 
            action="store_true",
            help="Input features have associated timestamps.")
    parser.add_argument("-R", "--Rady",
            action="store_true",
            help="Prospective data comes from Rady Clinithink process.")
    parser.add_argument("-C", "--Cardinal",
            action="store_true",
            help="Return cardinal phenotypes for prospective cases.")
    parser.add_argument("-P", "--Pickle", 
            action="store_true",
            help="Dump pickled model object to file '{outdir}/trained_model.pickle'")
    parser.add_argument("--json",
            action="store_true",
            help="Return results as JSON object.")
    parser.add_argument("-o", "--outdir", 
            default="analysis/test",
            help="Output directory for results & reports.")
    return parser


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
    # BasicHPOSet does the following:
    # -removes parent terms
    # -removes modifier terms
    # -replaces obsolete terms
    hpo_set = BasicHPOSet.from_queries(hpo_lst)
    hpo_dic = hpo_set.toJSON()
    hpo_str = sorted([x["id"] for x in hpo_dic])
    return hpo_str


def clean_codes(codes, valid_hpo, keep_others=False):
    """Basic pre-processing of input codes in preparation for modeling.

    Args:
        codes (iterable): An iterable of codes.
        keep_others (bool, optional): Whether to keep non-HPO/ICD-10-CM codes. Defaults to False.

    Returns:
        list: A filtered list of codes.
    """
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


def annotate_codes(cde, valid_hpo):
    if cde in valid_hpo:
        return Ontology.get_hpo_object(cde).name
    elif cm.is_valid_item(cde):
        return cm.get_description(cde)
    else:
        return "other"


def make_compliant(data, valid_hpo, dataset_name, col_idx, check_cols=None, keep_all_codes=False):
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
        clean = clean_codes(dirty, valid_hpo, keep_others=keep_all_codes)
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
        #('SVM', SVC()),
        ('SVM', SVC(probability=True)),
        ('GBM', GradientBoostingClassifier()),
        ('MLP', MLPClassifier()),
        ('BNB', BernoulliNB())
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
        y_score = clf.predict_proba(X_test)[:,1]
        print(name)
        print(metrics.classification_report(y_test, y_pred, target_names=target_names))
        print(metrics.roc_auc_score(y_test, y_score))
        print(metrics.confusion_matrix(y_test, y_pred))
        print()
        results.append(cv_results)
        names.append(name)
        this_df = pd.DataFrame(cv_results)
        this_df["model"] = name
        this_df.drop(columns=["fit_time","score_time"], inplace=True)
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

    prosp = make_compliant(prosp, valid_hpo, "prosp_data", prosp_col_idx, args.keep_all_codes)

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
    _ = Ontology()
    valid_hpo = Ontology.to_dataframe().index.tolist()
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
                valid_hpo,
                "train_data", 
                train_col_idx, 
                ["seq_status","diagnostic"],
                args.keep_all_codes)

        train_X = onehot_encode(train)
        keep_terms = train_X.columns
        train_y = np.array([x[train_col_idx["seq_status"]] for x in train[1:]])

        train_scores, train_preds = training(train_X, train_y)
        # rocy(train_preds, train_y)
        train_out = [x+y for x,y in zip(train, [preds_header] + train_preds.tolist())] 
        writey(train_out, path.join(args.outdir, "training_preds.tsv"))

        fit = BernoulliNB().fit(train_X, train_y)
        mod_attr = pd.concat([pd.DataFrame(fit.feature_names_in_).rename(columns={0:"codes"}),
                              pd.DataFrame(np.transpose(fit.feature_count_)).rename(columns={0:"ctrl_cnt", 1:"case_cnt"}),
                              pd.DataFrame(np.transpose(fit.feature_log_prob_)).rename(columns={0:"ctrl_log_prob", 1:"case_log_prob"})], axis = 1)
        mod_attr["coef"] = mod_attr["case_log_prob"] - mod_attr["ctrl_log_prob"]
        mod_attr["ctrl_n"] = fit.class_count_[0]
        mod_attr["case_n"] = fit.class_count_[1]
        mod_attr["code_descrip"] = mod_attr["codes"].apply(annotate_codes, valid_hpo=valid_hpo)
        mod_attr.to_csv(path.join(args.outdir, "feature_coefficients.tsv"), sep="\t", index=False)

        if args.Pickle:
            dump(fit, path.join(args.outdir, "trained_model.pickle"))

        if args.prospective:
            if args.compare_models:
                prosp = ready(args.prospective)
                prosp_col_idx = get_column_positions(prosp, ["pid","codes","seq_status"])
                prosp = make_compliant(prosp, valid_hpo, "prosp_data", prosp_col_idx, args.keep_all_codes)
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

    if args.prospective and args.json:
        JSON = build_JSON(prosp_out, get_column_positions(prosp_out, ["pid","scr"]), cards)
        fname = path.join(args.outdir, "prospective_preds.json")
        with open(fname, "w") as json_out:
            json.dump(JSON, json_out, indent=4)


if __name__ == "__main__":
    main()

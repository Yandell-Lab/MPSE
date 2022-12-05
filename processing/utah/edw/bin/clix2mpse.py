#!/usr/bin/env python3

import sys
import csv
import argparse
from datetime import datetime as dt


def argue():
	parser = argparse.ArgumentParser()
	parser.add_argument("-d", "--data", required=True)
	parser.add_argument("-M", "--Master", default="docs/NeoSeq_MASTER_all.tsv")
	parser.add_argument("--timestamps", action="store_true", 
			help="Add timestamps to HPO terms.")
	parser.add_argument("outfile", nargs="?", type=argparse.FileType("w"), default=sys.stdout,
			help="Provide optional name of output file. Default (-) writes to stdout.")
	args = parser.parse_args()
	return args


def ready(ftag, delim=",", drop_header=False):
	with open(ftag, newline="") as f:
		reader = csv.reader(f, delimiter=delim)
		if drop_header:
			out = [row for row in reader][1:]
		else:
			out = [row for row in reader]
	return out


def get_col_pos(data, col_names):
	header = [x.lower() for x in data[0]]
	idx_dic = {name: header.index(name) for name in col_names}
	return idx_dic


def generate_lookup(data, idx):
	lookup = {}
	for i in range(len(data)):
		lookup[data[i][idx]] = i
	return lookup


def main():
	args = argue()
	raw = ready(args.data, delim="\t")

	col_names = ["patient_id","dob","observation_datetime","project","abstractions"]
	data_col_pos = get_col_pos(raw, col_names)
	pid_idx = data_col_pos["patient_id"]
	dob_idx = data_col_pos["dob"]
	obs_idx = data_col_pos["observation_datetime"]
	prj_idx = data_col_pos["project"]
	abs_idx = data_col_pos["abstractions"]

	data = [row for row in raw if row[abs_idx] != ""]

	if args.timestamps:
		for row in data[1:]:
			odt = dt.fromisoformat(row[obs_idx])
			hpo = row[abs_idx].split(";")
			tsp = ["|".join([x, odt.strftime("%Y-%m-%d")]) for x in hpo]
			row[abs_idx] = ";".join(tsp)
	
	tmp = {}
	for row in data[1:]:
		if row[pid_idx] not in tmp.keys():
			tmp[row[pid_idx]] = {"project": row[prj_idx], 
					"dob": row[dob_idx], 
					"abstractions": [row[abs_idx]]}
		else:
			tmp[row[pid_idx]]["abstractions"].append(row[abs_idx])

	col_names = ["patient_id","dob","project","abstractions"]
	data = [col_names]
	for k, v in tmp.items():
		data.append([k, v["dob"], v["project"], ";".join(v["abstractions"])])
	data_col_pos = get_col_pos(data, col_names)
	pid_idx = data_col_pos["patient_id"]
	dob_idx = data_col_pos["dob"]
	prj_idx = data_col_pos["project"]
	abs_idx = data_col_pos["abstractions"]

	master = ready(args.Master, delim="\t")
	col_names = ["mrn","neoseq_id","diagnostic","edw_exists"]
	master_col_pos = get_col_pos(master, col_names)
	master = [x for x in master if x[master_col_pos["edw_exists"]] == "1"]
	neo_lookup = generate_lookup(master, master_col_pos["mrn"])

	mpse = []
	for row in data[1:]:
		mrn = row[pid_idx]
		if mrn in neo_lookup.keys():
			out = [mrn,
					master[neo_lookup[mrn]][master_col_pos["neoseq_id"]],
					row[dob_idx],
					"1",
					master[neo_lookup[mrn]][master_col_pos["diagnostic"]],
					row[abs_idx]]
		else:
			out = [mrn,
					row[prj_idx],
					row[dob_idx],
					"0",
					"0",
					row[abs_idx]]
		mpse.append(out)

	header = ["pid","project","dob","seq_status","diagnostic","codes"]
	args.outfile.write("\t".join(header) + "\n")
	for row in mpse:
		args.outfile.write("\t".join(row) + "\n")


if __name__ == "__main__":
	main()

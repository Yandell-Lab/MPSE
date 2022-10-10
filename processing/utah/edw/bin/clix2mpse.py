#!/usr/bin/env python3

import sys
import csv
import argparse
#from datetime import date
#from datetime import datetime as dt


def argue():
	parser = argparse.ArgumentParser()
	parser.add_argument("-d", "--data", required=True)
	parser.add_argument("-g", "--grouped", action="store_true")
	parser.add_argument("-N", "--NeoSeq", action="store_true")
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
	idx_dic = {name: data[0].index(name) for name in col_names}
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
	data = [row for row in raw if row[data_col_pos["abstractions"]] != ""]

	if args.NeoSeq:
		master = ready(args.Master, delim="\t")
		col_names = ["mrn","neoseq_id","diagnostic"]
		master_col_pos = get_col_pos(master, col_names)
		neo_lookup = generate_lookup(master[1:], master_col_pos["mrn"])

	mpse = []
	for row in data[1:]:
		mrn = row[data_col_pos["patient_id"]]
		dob = row[data_col_pos["dob"]],
		odt = row[data_col_pos["observation_datetime"]],
		hpo = row[data_col_pos["abstractions"]]
		if args.NeoSeq:
			out = [mrn,
					master[neo_lookup[mrn]][master_col_pos["neoseq_id"]],
					dob,
					odt,
					"1",
					master[neo_lookup[mrn]][master_col_pos["diagnostic"]],
					hpo]
		else:
			out = [mrn,
					row[data_col_pos["project"]],
					dob,
					odt,
					"0",
					"0",
					hpo]
		mpse.append(out)

	header = ["pid","project","dob","obs_datetime","seq_status","diagnostic","hpo"]
	args.outfile.write("\t".join(header) + "\n")
	for row in mpse:
		args.outfile.write("\t".join(row) + "\n")


if __name__ == "__main__":
	main()

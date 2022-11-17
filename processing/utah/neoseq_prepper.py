#!/usr/bin/env python3

import sys
import csv
import re
import argparse
from hashlib import sha256
from pyhpo.ontology import Ontology
from pyhpo.set import HPOSet


def argue():
	parser = argparse.ArgumentParser()
	parser.add_argument("-n", "--neoseq", 
			default="data/utah/processed/neoseq_cross_demographs.tsv",
			help="File containing NeoSeq DW_PIDs with sequenced/diagnostic status.")
	parser.add_argument("-H", "--hpo",
			default="../FTP/data/edw_data_lemmon/hpo_ids.tsv",
			help="File of EDW HPO term lists.")
	parser.add_argument("outfile", nargs="?", type=argparse.FileType("w"), default=sys.stdout,
			help="Provide optional name of output file. Default (-) writes to stdout.")
	args = parser.parse_args()
	return args

def ready(ftag, delim="\t", drop_header=False):
	with open(ftag, newline="") as f:
		reader = csv.reader(f, delimiter=delim)
		if drop_header:
			out = [row for row in reader][1:]
		else:
			out = [row for row in reader]
	return out

def hasher(dw_pid):
	m = sha256()
	m.update(int(dw_pid).to_bytes(4, "big"))
	return m.hexdigest()

def generate_lookup(f, c_idx):
	lookup = {}
	for i in range(len(f)):
		lookup[f[i][c_idx]] = i
	return lookup

def inner_join(f1, f2, c1_idx, c2_idx):
	f1_lookup = generate_lookup(f1, c1_idx)
	f2_lookup = generate_lookup(f2, c2_idx)
	f2_ncol = len(f2[0])
	
	f = []
	header = f1[0] + f2[0][:c2_idx] + f2[0][c2_idx+1:]
	f.append(header)
	for row in f1[1:]:
		join_id = row[c1_idx]
		if join_id in f2_lookup.keys():
			match = f2[f2_lookup[row[c1_idx]]]
			add = row + match[:c2_idx] + match[c2_idx+1:]
			f.append(add)
	return f

def hpo_parse(hpo_str):
	hpo_reg = re.compile(r"HP:\d{7}")
	hpo_split = hpo_reg.findall(hpo_str)
	return hpo_split

def main():
	args = argue()
	_ = Ontology()

	neo = ready(args.neoseq, drop_header=True)
	hpo = ready(args.hpo, drop_header=True)
	hpo = [[x[0], x[5]] for x in hpo]
	
	for row in neo:
		row.append(hasher(row[0]))
	
	joined = inner_join(neo, hpo, 4, 0)

	for row in joined:
		row[5] = ";".join(hpo_parse(row[5]))
		row.append(row[0])
	
	header = ["DW_PID","PAT_ID","seq_status","diagnostic","hashed_id","codes","pid"]
	args.outfile.write("\t".join(header) + "\n")
	for row in joined:
		args.outfile.write("\t".join(row) + "\n")


if __name__ == "__main__":
	main()

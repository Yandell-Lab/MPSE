#!/usr/bin/env python3

import sys
import re
import csv
import json
import argparse
import pypandoc

csv.field_size_limit(sys.maxsize)


def argue():
	parser = argparse.ArgumentParser()
	parser.add_argument("-d", "--data", required=True)
	parser.add_argument("-t", "--type", choices=["clix","clinphen"], required=True,
			help="Input data source type.")
	parser.add_argument("outfile", nargs="?", type=argparse.FileType("w", encoding="ascii"), default=sys.stdout, help="Provide optional name of output file. Default (-) writes to stdout.")
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


def clean_text(txt, rtf_regex):
	if rtf_regex.match(txt) is not None:
		txt = pypandoc.convert_text(txt, "plain", format="rtf")
		txt = txt.replace("|", "")
	txt = re.sub(r'[^\x00-\x7f]', ' ', txt)
	txt = re.sub(r'\n', ' ', txt)
	txt = re.sub(r' +', ' ', txt)
	return txt


def main():
	args = argue()
	data = ready(args.data, delim="\t")
	rtf_regex = re.compile(r'^{\\rtf\d')
	col_names = ["pat_id","birth_date","rpt_id","proc_dt","rpt_type_cd","text"]
	col_pos = get_col_pos(data, col_names)

	if args.type == "clix":
		docs = []
		for row in data[1:]:
			txt = row[col_pos["text"]]
			if txt == "":
				continue
			elif len(txt) > 50000:
				continue
			else:
				txt = clean_text(txt, rtf_regex)
				meta = {"patient_id": row[col_pos["pat_id"]],
						"document_id": row[col_pos["rpt_id"]],
						"patient_dob": row[col_pos["birth_date"]],
						"observation_datetime": row[col_pos["proc_dt"]],
						"project": "neoseq4mpse",
						"visit_id": row[col_pos["rpt_type_cd"]],
						"author": "bennet_peterson"}
				docs.append({"data": txt, "metadata": meta})
		out = {"documents": docs}
		args.outfile.write(json.dumps(out, indent=4))

	elif args.type == "clinphen":
		header = ["pid","document_id","patient_dob","observation_datetime","visit_id","text"]
		args.outfile.write("\t".join(header) + "\n")
		for row in data[1:]:
			txt = row[col_pos["text"]]
			if txt == "":
				continue
			elif len(txt) > 50000:
				continue
			else:
				txt = clean_text(txt, rtf_regex)
				out = [row[col_pos["pat_id"]],
					row[col_pos["rpt_id"]],
					row[col_pos["birth_date"]],
					row[col_pos["proc_dt"]],
					row[col_pos["rpt_type_cd"]],
					txt]
				args.outfile.write("\t".join(out) + "\n")


if __name__=="__main__":
	main()

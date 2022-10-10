#!/usr/bin/env python3

import os
import sys
import re
import json
import argparse
from jinja2 import Template
import cx_Oracle
import pypandoc
from common import get_db


def argue():
	parser = argparse.ArgumentParser()
	parser.add_argument("-s", "--sql", required=True)
	parser.add_argument("-p", "--pat_id", type=ascii, required=True)
	parser.add_argument("outfile", nargs="?", type=argparse.FileType("w", encoding="ascii"), default=sys.stdout, help="Provide optional name of output file. Default (-) writes to stdout.")
	args = parser.parse_args()
	return args


def serialize_cell(cell):
	if cell is None:
		return ""
	elif isinstance(cell, cx_Oracle.Object):
		text = " ".join(cell.aslist())
	elif isinstance(cell, bytes):
		text = str(int.from_bytes(cell[0:8], byteorder="little"))
	else:
		text = str(cell)
	return " ".join(text.split())


def format_note(cursor):
	rtf_regex = re.compile(r'^{\\rtf\d')
	docs = []
	for row in cursor:
		serial = [x for x in map(serialize_cell, row)]
		if serial[4] == "":
			continue
		elif len(serial[4]) > 50000:
			continue
		elif rtf_regex.match(serial[4]) is not None:
			serial[4] = pypandoc.convert_text(serial[4], "plain", format="rtf")
			serial[4] = serial[4].replace("|", "")
		serial[4] = re.sub(r'[^\x00-\x7f]', ' ', serial[4])
		serial[4] = re.sub(r'\n', ' ', serial[4])
		serial[4] = re.sub(r' +', ' ', serial[4])
		meta = {"patient_id": serial[0],
				"document_id": serial[0],
				"patient_dob": serial[1],
				"observation_datetime": serial[2],
				"project": "neoseq4mpse",
				"visit_id": serial[3],
				"author": "bennet_peterson"}
		docs.append({"data": serial[4], "metadata": meta})
	return {"documents": docs}


def main():
	args = argue()
	sql_file = args.sql
	tsv_file = sql_file.replace('sql','tsv')
	template = Template(open(sql_file).read())

	username = os.environ.get("EDW_USER")
	password = os.environ.get("EDW_PASS")
	query = template.render(PAT_ID=args.pat_id)

	try:
		cursor = get_db(username, password).cursor()
		cursor.arraysize = 1000
		cursor.execute(query)
		#headers = map(lambda metadata: metadata[0], cursor.description)
		output = format_note(cursor)
		args.outfile.write(json.dumps(output, indent=4))

	except Exception as e:
		print("Error getting " + tsv_file)
		raise e


if __name__=="__main__":
	main()

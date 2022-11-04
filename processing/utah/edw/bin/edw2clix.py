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
		if serial[5] == "":
			continue
		elif len(serial[5]) > 50000:
			continue
		elif rtf_regex.match(serial[5]) is not None:
			serial[5] = pypandoc.convert_text(serial[5], "plain", format="rtf")
			serial[5] = serial[5].replace("|", "")
		serial[5] = re.sub(r'[^\x00-\x7f]', ' ', serial[5])
		serial[5] = re.sub(r'\n', ' ', serial[5])
		serial[5] = re.sub(r' +', ' ', serial[5])
		meta = {"patient_id": serial[0],
				"document_id": serial[2],
				"project": "edw4mpse",
				"visit_id": serial[4],
				"author": "bennet_peterson"}
		docs.append({"data": serial[5], "metadata": meta})
	return {"documents": docs}


def main():
	args = argue()
	sql_file = args.sql
	tsv_file = sql_file.replace('sql','tsv')
	template = Template(open(sql_file).read())

	username = os.environ.get("EDW_USER")
	password = os.environ.get("EDW_PASS")

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

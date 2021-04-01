#!/usr/bin/env python3

import sys
import argparse
import csv


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data', action='store', required=True,
		help='CSV file containing data for targets and background.')
args = parser.parse_args()


def reader(ftag):
	llist = []
	with open(ftag, 'r', newline='') as f:
		csvreader = csv.reader(f, delimiter=',')
		for row in csvreader:
			llist.append(row)
	return llist

def hpo_parse(hpo_lst):
	hpo_split = hpo_lst.split(':')
	hpo_split = [x.replace('hp', 'HP:', 1) for x in hpo_split]
	hpo_split = [x.replace('_', '\t', 1) for x in hpo_split]
	hpo_str = '\n'.join(hpo_split)
	return hpo_str


def printer(llist):
	with open('RADY_HPO_SEQ_ADMITS.fasta.txt', 'w') as targ, open('RADY_HPO_ALL_ADMITS.fasta.txt', 'w') as back:
		for row in llist[1:]:
			if row[13]=='':
				b0 = '>' + row[12]
				b1 = 'RID:' + row[0]
				b2 = 'NAME:' + row[0] + '_' + row[10]
				b3 = 'AGE:' + row[7]
				b  = '\t'.join([b0, b1, b2, b3])
				bseq = hpo_parse(row[11])
				back.write(b + '\n')
				back.write(bseq + '\n')
			else:
				t0 = '>' + row[22]
				t1 = 'RID:' + row[13]
				t2 = 'NAME:' + row[13] + '_' + row[20]
				t3 = 'AGE:' + row[14]
				t4 = 'DIAG:' + ("0" if row[18]=="True" else \
						"1" if row[15]=="True" else "0")
				t  = '\t'.join([t0, t1, t2, t3, t4])
				tseq = hpo_parse(row[21])
				targ.write(t + '\n')
				targ.write(tseq + '\n')

def main():
	printer(reader(args.data))

if __name__ == '__main__':
	main()

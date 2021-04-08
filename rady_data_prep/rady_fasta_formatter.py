#!/usr/bin/env python3

import sys
import argparse
import pandas as pd
from pyhpo.ontology import Ontology


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data', action='store', required=True,
		help='CSV file containing data for targets and background.')
args = parser.parse_args()


def reader(ftag):
	with open(ftag, 'r') as f:
		df = pd.read_csv(f, 
				dtype={
					"all_fid":"Int64",
					"seq_fid":"Int64",
					"all_file_tag":"Int64", 
					"seq_file_tag":"Int64",
					"Age5d":"Int64",
					"Age":"Int64"}, 
				keep_default_na=False)
	return df

def hpo_parse(hpo_lst):
	hpo_split = hpo_lst.split('_')
	hpo_name = [x + "\t" + Ontology.get_hpo_object(x).name for x in hpo_split]
	hpo_str = '\n'.join(hpo_name)
	return hpo_str

def printer(d):
	with open('RADY_HPO_SEQ_ADMITS.fasta.txt', 'w') as targ, open('RADY_HPO_ALL_ADMITS.fasta.txt', 'w') as back:
		for row in d.itertuples(index=False):
			if row.ResearchID=="":
				b0 = ">" + str(row.all_fid)
				b1 = "RID:" + str(row.CT_HPO_FileID)
				b2 = 'NAME:' + str(row.CT_HPO_FileID) + '_' + str(row.all_file_tag)
				b3 = 'AGE:' + str(row.Age5d)
				b  = '\t'.join([b0, b1, b2, b3])
				bseq = hpo_parse(row.all_HPO_clean)
				back.write(b + '\n')
				back.write(bseq + '\n')
			else:
				t0 = '>' + str(row.seq_fid)
				t1 = 'RID:' + str(row.ResearchID)
				t2 = 'NAME:' + str(row.ResearchID) + '_' + str(row.seq_file_tag)
				t3 = 'AGE:' + str(row.Age)
				t4 = 'DIAG:' + ("0" if row.Incidental=="True" else \
						"1" if row.Positive=="True" else "0")
				t  = '\t'.join([t0, t1, t2, t3, t4])
				tseq = hpo_parse(row.seq_HPO_clean)
				targ.write(t + '\n')
				targ.write(tseq + '\n')


def main():
	_ = Ontology()
	printer(reader(args.data))

if __name__ == '__main__':
	main()

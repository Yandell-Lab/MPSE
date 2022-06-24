#!/usr/bin/env python

import sys
import re
import csv


def ready(ftag, delim=",", drop_header=False):
    with open(ftag, newline="") as f:
        reader = csv.reader(f, delimiter=delim)
        if drop_header:
            out = [row for row in reader][1:]
        else:
            out = [row for row in reader]
    return out

def writey(data, ftag, header=None, delim="\t"):
    with open(ftag, "w", newline="") as f:
        writer = csv.writer(f, delimiter=delim)
        if header is not None:
            writer.writerow(header)
        writer.writerows(data)
    return True

def hpo_parse(hpo_str):
	hpo_reg = re.compile(r"hp\d{7}")
	hpo_split = [x.replace("hp", "HP:", 1) for x in hpo_reg.findall(hpo_str)]
	return ";".join(hpo_split)


def main():
	all_f = "ALL_ADMITS_EXP/data_views/rady_hpo_all_admits_remove_repeats.csv"
	seq_f = "SEQ_ADMITS/data_views/updated/rady_hpo_seq_admits_CliniThink.csv"
	
	all_d = ready(all_f, drop_header=True)
	seq_d = ready(seq_f, drop_header=True)
	
	all_sub = [[x[0],x[4],x[5],x[6],x[7],"0","0","0",hpo_parse(x[11])] for x in all_d if not x[9]]
	#seq_sub = [[x[0],"","","",x[1],"1",x[2],x[5],hpo_parse(x[8])] for x in seq_d]
	seq_sub = []
	for x in seq_d:
		if x[0] == "UR225":
			# UR225 was mistakenly coded as non-diagnostic
			seq_sub.append([x[0], "", "", "", x[1], "1", "1", x[5], hpo_parse(x[8])])
		else:
			seq_sub.append([x[0], "", "", "", x[1], "1", x[2], x[5], hpo_parse(x[8])])

	header = ["pid","sex","race","ethnicity","age","seq_status","diagnostic","incidental","hpo"]
	writey(all_sub + seq_sub, "rady_training_data.tsv", header=header)


if __name__ == "__main__":
	main()

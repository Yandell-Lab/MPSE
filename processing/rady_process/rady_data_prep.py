#!/usr/bin/env python

import sys
import re
import csv


def ready(ftag, delim=",", header=True):
    with open(ftag, newline="") as f:
        reader = csv.reader(f, delimiter=delim)
        if header:
            out = [row for row in reader][1:]
        else:
            out = [row for row in reader]
    return out


def hpo_parse(hpo_str):
	hpo_reg = re.compile(r"hp\d{7}")
	hpo_split = [x.replace("hp", "HP:", 1) for x in hpo_reg.findall(hpo_str)]
	return ";".join(hpo_split)


def main():
    all_f = "ALL_ADMITS_EXP/data_views/rady_hpo_all_admits_remove_repeats.csv"
    seq_f = "SEQ_ADMITS/data_views/rady_hpo_seq_admits_CliniThink.csv"

    all_d = ready(all_f)
    seq_d = ready(seq_f)

    all_sub = [[x[0],x[4],x[5],x[6],x[7],"0","","",hpo_parse(x[11])] for x in all_d if not x[9]]
    seq_sub = [[x[0],"","","",x[1],"1",x[2],x[5],hpo_parse(x[8])] for x in seq_d]

    header = ["pid","sex","race","ethnicity","age","seq_status","diagnostic","incidental","hpo"]

    with open("rady_training_data.csv", "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(header)
        writer.writerows(all_sub + seq_sub)

if __name__ == "__main__":
	main()

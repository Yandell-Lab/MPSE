#!/usr/env/bin bash


cd ./rady_data_prep/ALL_ADMITS_EXP
./rady_hpo_all_admits_prep.sh

cd ../SEQ_ADMITS
./rady_hpo_seq_admits_prep.sh

cd ..
./join_seq_all_admits.sh

./rady_fasta_formatter.py -d ./rady_hpo_all_seq_joined.csv

cd ../analysis
./table_maker.py -d ../rady_data_prep/rady_hpo_all_seq_joined.csv > foo.clean

cd ..
../GRAVITY/bin/mpse -s -t ./analysis/foo.clean -q ./rady_data_prep/RADY_HPO_SEQ_ADMITS.fasta.txt -m -j t > tar
../GRAVITY/bin/mpse -s -t ./analysis/foo.clean -q ./rady_data_prep/RADY_HPO_ALL_ADMITS.fasta.txt -m -j b > bar

../GRAVITY/bin/mpse_stats -t tar -b bar -p 0.5

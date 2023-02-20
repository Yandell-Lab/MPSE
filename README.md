# Mendelian Phenotype Search Engine
MPSE employs HPO-based phenotype descriptions derived from patient EHRs to compute a score. This score can be used to determine the likelihood that a Mendelian condition is contributing to a patient’s clinical presentation, and thus, can be used for the prioritization of patients for WGS.


## CLI Arguments
The following core parameters are available while executing mpse.py
- -t,  --training ... Case/control training data in standard format
- -m,  --model ... Serialized model (pickle object) to load from disc
- -p,  --prospective ... Prospective data in standard format
- -C,  --Cardinal ... Return cardinal phenotypes for prospective data
- -P,  --Pickle ... Dump pickled model object to file
- -o,  --outdir ... Output directory for results and reports


## Example Commands
To get the full usage statement, use the -h/--help flag:  
> $ ./bin/mpse.py -h

The most basic use is training a model and returning scores for the training cohort:  
> $ ./bin/mpse.py -t data/test/fake_training_data.tsv  

By default, output files are written to './analysis/test/'. This command will create a single file named ./analysis/test/training_preds.tsv

Prospective cases can be scored using the trained model:  
> $ ./bin/mpse.py -t data/test/fake_training_data.tsv -p data/test/fake_prospective_data.tsv  

The prospective results are sent to standard output

Setting the -P/--Pickle flag will write the model object to disc in pickle format:  
> $ ./bin/mpse.py -t data/test/fake_training_data.tsv -P  

A file named trained_model.pickle will appear in the output directory

This pickle object can then be used in place of a training data file to score prospective cases:  
> $ ./bin/mpse.py -m analysis/test/trained_model.pickle -p data/test/fake_prospective_data.tsv

The -C/--Cardinal flag generates a file containing cardinal phenotypes for each prospective case:  
> $ ./bin/mpse.py -m analysis/test/trained_model.pickle -p data/test/fake_prospective_data.tsv -C  

Cardinal phenotypes are written to a file named cardinal_phenotypes.tsv


## Input File Format
MPSE demands certain characteristics of input data. The data must be tab-delimited with rows corresponding to patients and columns corresponding to features/variables. The first row must be a header line with feature names. The following features must be included (in no particular order) with specified names.

- pid - unique observation identifier
- seq_status - sequencing status (0=not sequenced, 1=sequenced)
- diagnostic - diagnosis status (0=not diagnostic, 1=diagnostic)
- incidental - incidental finding (0=no, 1=yes)
- codes - semicolon-delimited list of hpo terms
    - ex. "HP:0000001;HP:0000002;HP:0000003"

There is an exception for prospective data, when *seq_status*, *diagnostic*, and *incidental* features may be unavailable or not applicable. Any additional features will be ignored by MPSE.


## Output File Format
Output files are made by appending 7 new fields to the original input files. These fields are:
- codes_clean - list of hpo terms with all parent terms removed (terminal terms retained)
- neg_proba - predicted probability data comes from negative class (not sequenced)
- pos_proba - predicted probability data comes from positive class (sequenced)
- neg_log_proba - natural log of neg_proba
- pos_log_proba	- natural log of pos_proba
- class - model's class prediction
- scr - MPSE score = log(pos_proba / neg_proba)

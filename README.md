# Mendelian Phenotype Search Engine


## CLI Arguments
The following arguments are available while executing mpse.py
- -t, --training        Case/control training data in standard format (required)
- -v, --validate        Validation data in standard format (optional)
- -o, --outdir          Output directory for results and reports (optional)


## Input File Format
MPSE demands certain characteristics of input data. The data must be delimited with rows corresponding to unique observations and columns corresponding to features. The first row must be a header line with feature names. The following features must be included (in no particular order) with specified names.

- pid - unique observation identifier
- sex
- race
- ethnicity
- age
- seq_status - sequencing status (0=not sequenced, 1=sequenced)
- diagnostic - diagnosis status (0=not diagnostic, 1=diagnostic)
- incidental - incidental finding (0=no, 1=yes)
- hpo - delimited list of hpo terms
    - ex. "HP:0000001;HP:0000002;HP:0000003"

There is an exception for prospective data, when *seq_status*, *diagnostic*, and *incidental* features are unavailable. Any additional features will be ignored by MPSE.

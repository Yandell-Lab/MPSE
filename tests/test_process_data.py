#!/usr/bin/env python3

import os
import sys
import pytest

parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, parent_dir)
from bin.mpse import *


@pytest.fixture(scope="module")
def example_codes():
    # 806-0 LOINC
    # RX10359383 RxNorm
    # I70.501 ICD-10-CM (with dot)
    # HP:0001182 HPO (tapered finger)
    # HP:0001166 HPO (arachnodactyly)
    # HP:9999999 HPO (invalid code)
    # C84Z0 ICD-10-CM (without dot)
    return ["806-0","RX10359383","I70.501","HP:0001182","HP:0001166","HP:9999999","C84Z0"]


@pytest.fixture(scope="module")
def hponto():
    return Ontology()


def test_get_column_positions():
    data = [["col1","col2","col3"],
            ["aaaa","bbbb","cccc"]]
    col_names = ["col1","col3"]
    assert get_column_positions(data, col_names) == {"col1": 0, "col3": 2}


@pytest.mark.parametrize("hpo_code", ["hp0000001","hp:0000001","HP0000001","HP:0000001"])
def test_parse_hpo(hpo_code):
    assert parse_hpo(hpo_code) == "HP:0000001"


def test_extract_timestamps():
    data = [["pid","codes","abc"],
            ["id1","HP:0000001|2000-01-01;HP:0000002|2000-01-02","a"],
            ["id2","HP:0000003|2000-01-03","b"],
            ["id3","","c"]]
    col_idx = {"pid": 0, "codes": 1}
    # extract = [["pid","codes","abc","manifest_date"],
    #            ["id1","HP:0000001","a","2000-01-01"],
    #            ["id1","HP:0000002;HP:0000001","a","2000-01-02"], ## order of codes is non-deterministic
    #            ["id2","HP:0000003","b","2000-01-03"],
    #            ["id3","","c",""]]
    extract = extract_timestamps(data, col_idx)
    assert len(extract) == 5
    assert extract[0] == ["pid","codes","abc","manifest_date"]
    assert extract[1][1] == "HP:0000001"
    assert len(extract[2][1]) == 21
    assert extract[2][3] == "2000-01-02"
    assert extract[3][2] == "b"
    assert extract[4][3] == ""


def test_remove_parent_terms(hponto):
    terms = ["HP:0001182","HP:0001166","HP:0100807","HP:0001167"]
    assert remove_parent_terms(terms) == ["HP:0001166","HP:0001182"]


def test_clean_codes_with_keep_others(hponto, example_codes):
    valid_hpo = Ontology.to_dataframe().index.tolist()
    expected = ["HP:0001166", "HP:0001182", "C84.Z0", "I70.501",
                "806-0", "HP:9999999", "RX10359383"]
    assert clean_codes(example_codes, valid_hpo, True) == expected


def test_clean_codes_no_keep_others(hponto, example_codes):
    valid_hpo = Ontology.to_dataframe().index.tolist()
    expected = ["HP:0001166", "HP:0001182", "C84.Z0", "I70.501"]
    assert clean_codes(example_codes, valid_hpo, False) == expected


def test_annotate_codes(hponto, example_codes):
    valid_hpo = Ontology.to_dataframe().index.tolist()
    annos = []
    for cde in example_codes:
        annos.append(annotate_codes(cde, valid_hpo, add_vars=[]))
    expected = ["other",
                "other",
                "Unspecified atherosclerosis of nonautologous biological bypass graft(s) of the extremities, right leg",
                "Tapered finger",
                "Arachnodactyly",
                "other",
                "Other mature T/NK-cell lymphomas, unspecified site"]
    assert annos == expected


def test_compliant_data(hponto, example_codes):
    valid_hpo = Ontology.to_dataframe().index.tolist()
    data = [["pid","diagnostic","codes","abc"],
            ["id1","1",";".join(example_codes),"a"],
            ["id2","","HP:0001167;HP:0001167","b"]]
    dataset_name = "my_compliant_data"
    col_idx = {"pid": 0, "diagnostic": 1, "codes": 2}
    expected = [["pid","diagnostic","codes","abc","codes_clean"],
                ["id1","1",";".join(sorted(set(example_codes))),"a","HP:0001166;HP:0001182;C84.Z0;I70.501"],
                ["id2","0","HP:0001167","b","HP:0001167"]]
    assert make_compliant(data, valid_hpo, dataset_name, col_idx, ["diagnostic"], False) == expected


def test_non_compliant_data(example_codes):
    valid_hpo = Ontology.to_dataframe().index.tolist()
    data = [["pid","diagnostic","codes","abc"],
            ["id1","*",";".join(example_codes),"a"],
            ["id2","","HP:0001167;HP:0001167","b"]]
    dataset_name = "my_non_compliant_data"
    col_idx = {"pid": 0, "diagnostic": 1, "codes": 2}
    with pytest.raises(ValueError):
        make_compliant(data, valid_hpo, dataset_name, col_idx, ["diagnostic"], False)


def test_onehot_encode():
    data = [["pid","abc","codes_clean","cov"],
            ["id1","aaa","HP:0000001",1],
            ["id2","bbb","HP:0000001;HP:0000002",1],
            ["id3","ccc","HP:0000002;HP:0000003",1]]
    expected = pd.DataFrame([[1,0,0],[1,1,0],[0,1,1]],
                            columns=["HP:0000001","HP:0000002","HP:0000003"])
    expected_cov = pd.DataFrame([[1,0,0,1],[1,1,0,1],[0,1,1,1]],
                                columns=["HP:0000001","HP:0000002","HP:0000003","cov"])
    assert onehot_encode(data, add_vars=[]).equals(expected)
    assert onehot_encode(data, add_vars=["cov"]).equals(expected_cov)


if __name__ == "__main__":
    pytest.main()

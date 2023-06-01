#!/usr/bin/env python3

import os
import sys
import pytest

parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, parent_dir)
from bin.mpse import *


@pytest.fixture
def example_codes():
    # 806-0 LOINC
    # RX10359383 RxNorm
    # I70.501 ICD-10-CM (with dot)
    # C84Z0 ICD-10-CM (without dot)
    # HP:9999999 HPO (invalid code)
    return ["806-0","RX10359383","I70.501","HP:0001182","HP:0001166","HP:9999999","C84Z0"]


@pytest.fixture
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
    extract = [["pid","codes","abc","manifest_date"],
               ["id1","HP:0000001","a","2000-01-01"],
               ["id1","HP:0000002;HP:0000001","a","2000-01-02"],
               ["id2","HP:0000003","b","2000-01-03"],
               ["id3","","c",""]]
    assert extract_timestamps(data, col_idx) == extract


def test_remove_parent_terms(hponto):
    terms = ["HP:0001182","HP:0001166","HP:0100807","HP:0001167"]
    assert remove_parent_terms(terms) == ["HP:0001166","HP:0001182"]


def test_clean_codes(hponto, example_codes):
    assert clean_codes(example_codes, True) == ["HP:0001166",
                                                "HP:0001182",
                                                "C84.Z0",
                                                "I70.501",
                                                "806-0",
                                                "HP:9999999",
                                                "RX10359383"]
    assert clean_codes(example_codes, False) == ["HP:0001166",
                                                 "HP:0001182",
                                                 "C84.Z0",
                                                 "I70.501"]


def test_compliant_data(hponto, example_codes):
    data = [["pid","diagnostic","codes","abc"],
            ["id1","1",";".join(example_codes),"a"],
            ["id2","","HP:0001167;HP:0001167","b"]]
    dataset_name = "my_compliant_data"
    col_idx = {"pid": 0, "diagnostic": 1, "codes": 2}
    expected = [["pid","diagnostic","codes","abc","codes_clean"],
                ["id1","1",";".join(sorted(set(example_codes))),"a","HP:0001166;HP:0001182;C84.Z0;I70.501"],
                ["id2","0","HP:0001167","b","HP:0001167"]]
    assert make_compliant(data, dataset_name, col_idx, ["diagnostic"], False) == expected


def test_non_compliant_data(example_codes):
    data = [["pid","diagnostic","codes","abc"],
            ["id1","*",";".join(example_codes),"a"],
            ["id2","","HP:0001167;HP:0001167","b"]]
    dataset_name = "my_non_compliant_data"
    col_idx = {"pid": 0, "diagnostic": 1, "codes": 2}
    with pytest.raises(ValueError):
        make_compliant(data, dataset_name, col_idx, ["diagnostic"], False)


def test_onehot_encode():
    data = [["pid","abc","codes_clean"],
            ["id1","aaa","HP:0000001"],
            ["id2","bbb","HP:0000001;HP:0000002"],
            ["id3","ccc","HP:0000002;HP:0000003"]]
    expected = pd.DataFrame([[1,0,0],[1,1,0],[0,1,1]], columns=["HP:0000001","HP:0000002","HP:0000003"])
    assert onehot_encode(data).equals(expected)


if __name__ == "__main__":
    pytest.main()

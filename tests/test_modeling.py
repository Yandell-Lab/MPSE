#!/usr/bin/env python3

import os
import sys
import pytest
from sklearn.datasets import make_classification

parent_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, parent_dir)
from bin.mpse import *


@pytest.fixture(scope="module")
def training_data():
    # Generate dummy data for testing
    X, y = make_classification(n_samples=100, n_features=10, random_state=42)
    return X, y


def test_training(training_data):
    scores, results = training(training_data[0], training_data[1])

    assert isinstance(scores, dict)
    assert isinstance(results, np.ndarray)
    assert results.shape[0] == training_data[0].shape[0]
    assert results.shape[1] == 6  # Number of columns in the concatenated results
    assert np.array_equal(np.unique(training_data[1]), np.unique(results[:, -2])) # Classes column contains correct unique classes


def test_score_probands(training_data):
    model = BernoulliNB()
    model.fit(training_data[0], training_data[1])

    valid_X, _ = make_classification(n_samples=10, n_features=10, random_state=42)
    results = score_probands(model, valid_X)

    assert isinstance(results, np.ndarray)
    assert results.shape[0] == valid_X.shape[0]
    assert results.shape[1] == 6  # Number of columns in the concatenated results


def test_process_prospective(training_data, tmp_path):
    _ = Ontology()
    valid_hpo = Ontology.to_dataframe().index.tolist()
    model = BernoulliNB()
    model.fit(training_data[0], training_data[1])

    keep_terms = [chr(ord("A") + i) for i in range(10)]
    header = ["neg_proba","pos_proba","neg_log_proba","pos_log_proba","class","scr"]
    prospective_data = "pid\tcodes\nid1\tA;B;D;H\nid2\tA;B;G;I\nid3\tA;B;C;G;H;I;K;L"
    # prospective_data = [["pid","codes"],
    #                     ["id1","A;B;D;H"],
    #                     ["id2","A;B;G;I"],
    #                     ["id3","A;B;C;G;H;I;K;L"]]

    temp_file = tmp_path / "test_prospective.csv"
    temp_file.write_text(prospective_data)

    parser = argparse.ArgumentParser()
    parser.add_argument("-R", "--Rady", action="store_true", default=False)
    parser.add_argument("-p", "--prospective")
    parser.add_argument("--timestamps", action="store_true", default=False)
    parser.add_argument("--vars", action="append", default=[])
    parser.add_argument("--keep_all_codes", action="store_true", default=False)
    args = parser.parse_args(["--prospective", str(temp_file), "--keep_all_codes"])

    prosp, prosp_X, prosp_out = process_prospective(model, valid_hpo, keep_terms, header, args)

    assert isinstance(prosp, list)
    assert len(prosp) == 4
    assert prosp[0] == ["pid","codes","codes_clean"]
    assert len(prosp[1]) == 3

    assert isinstance(prosp_X, pd.DataFrame)
    assert prosp_X.shape[0] == len(prosp) - 1
    assert prosp_X.shape[1] == len(keep_terms)

    assert prosp_out[0] == ["pid","codes","codes_clean"] + header
    assert isinstance(prosp_out, list)
    assert len(prosp_out) == len(prosp)


if __name__ == "__main__":
    pytest.main()

"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_databases_from_extann.py -x -v --log-cli-level=DEBUG --capture=tee-sys
coverage report --include=howard/* -m
"""

import os
import argparse
import logging as log
import pytest
import json
from os.path import join as osj
from tempfile import TemporaryDirectory
from howard.functions.commons import remove_if_exists, identical
from howard.functions.from_extann import from_extann, read_json
from test_needed import tests_folder


def edit_config(configfile, key, new, tmp_dir):
    data = read_json(configfile)
    data[key] = new
    with open(osj(tmp_dir, "param.tmp.json"), "w+") as js:
        json.dump(data, js)
    return osj(tmp_dir, "param.tmp.json")


@pytest.mark.parametrize("mode", ["all", "longest", "chosen"])
def test_from_extann(mode):
    """
    Test from extann func while switching on each mode
    """
    with TemporaryDirectory(dir="/tmp") as tmp_dir:
        input_ex = osj(tests_folder, "databases", "extann", "extann.tsv")
        output_ex = osj(tmp_dir, f"OMIM.output.{mode}.bed")
        param = osj(tests_folder, "data", "example.param.extann.json")
        param_ex = edit_config(param, "mode_extann", mode, tmp_dir)
        expected_ex = osj(
            tests_folder, "databases", "extann", f"extann.expected.{mode}.bed"
        )

        # remove
        remove_if_exists([output_ex])
        remove_if_exists([output_ex + ".hdr"])

        args = argparse.Namespace(
            input_extann=input_ex,
            output_extann=output_ex,
            param_extann=param_ex,
            refgene=osj(
                tests_folder,
                "databases",
                "extann",
                "ncbiRefSeq.chunk.bed.gz",
            ),
            transcripts=osj(tests_folder, "data", "transcripts.from_extann.tsv"),
            mode_extann=mode,
        )
        from_extann(args)

        assert os.path.exists(output_ex)
        assert os.path.exists(output_ex + ".hdr")
        assert identical([expected_ex, output_ex], "##")

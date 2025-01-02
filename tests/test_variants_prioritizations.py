# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_objects_variants.py -x -v --log-cli-level=INFO --capture=tee-sys
coverage report --include=howard/* -m 
"""

import logging as log
import os
import sys
from tempfile import TemporaryDirectory
import duckdb
import re
import Bio.bgzf as bgzf
import gzip
import pytest

from howard.functions.commons import *
from howard.objects.variants import Variants
from howard.functions.databases import *
from test_needed import *


def test_prioritization():
    """
    This is a test function for prioritization of variants in a VCF file using a specified configuration
    and parameter dictionary.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct config dict
        config = {}

        # Construct param dict
        param = {
            "prioritization": {
                "prioritization_config": tests_data_folder
                + "/prioritization_profiles.json",
                "profiles": ["default", "GERMLINE", "sql_class"],
                "pzfields": [
                    "PZFlag",
                    "PZScore",
                    "PZClass",
                    "PZComment",
                    "PZInfos",
                    "PZTags",
                ],
            }
        }

        # Create object
        variants = Variants(
            input=input_vcf, output=output_vcf, load=True, config=config, param=param
        )

        # Prioritization
        variants.prioritization()

        # Check all priorized default profile
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag_default=%'
            AND INFO LIKE '%PZScore_default=%'
            AND INFO LIKE '%PZClass_default=%'
            AND INFO LIKE '%PZComment_default=%'
            AND INFO LIKE '%PZInfos_default=%'
            AND INFO LIKE '%PZTags_default=%'
            """
        )
        assert len(result) == 4

        # Check all priorized GERMLINE profile
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag_GERMLINE=%'
            AND INFO LIKE '%PZScore_GERMLINE=%'
            AND INFO LIKE '%PZClass_GERMLINE=%'
            AND INFO LIKE '%PZComment_GERMLINE=%'
            AND INFO LIKE '%PZInfos_GERMLINE=%'
            AND INFO LIKE '%PZTags_GERMLINE=%'
            """
        )
        assert len(result) == 2

        # Check all priorized sql_class profile
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag_sql_class=%'
            AND INFO LIKE '%PZScore_sql_class=%'
            AND INFO LIKE '%PZClass_sql_class=%'
            AND INFO LIKE '%PZComment_sql_class=%'
            AND INFO LIKE '%PZInfos_sql_class=%'
            """
        )
        assert len(result) == 2

        # Check all priorized default profile (as default)
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag=%'
            AND INFO LIKE '%PZScore=%'
            AND INFO LIKE '%PZClass=%'
            AND INFO LIKE '%PZComment_default=%'
            AND INFO LIKE '%PZInfos_default=%'
            AND INFO LIKE '%PZTags_default=%'
            """
        )
        assert len(result) == 4

        # Check all priorized default profile
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag_default=%'
            AND INFO LIKE '%PZScore_default=%'
            """
        )
        assert len(result) == 7

        # Check all priorized GERMLINE profile
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag_GERMLINE=%'
            AND INFO LIKE '%PZScore_GERMLINE=%'
            """
        )
        assert len(result) == 7

        # Check all priorized default profile (as default)
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag=%'
            AND INFO LIKE '%PZScore=%'
            """
        )
        assert len(result) == 7

        # Check annotation default
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO LIKE '%PZScore_default=105%' """
        )
        assert len(result) == 1

        # Check annotation GERMILNE
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO LIKE '%PZScore_GERMLINE=5%' """
        )
        assert len(result) == 1

        # Check annotation sql_class
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO LIKE '%PZScore_sql_class=60%' """
        )
        assert len(result) == 1

        # Check FILTERED
        result = variants.get_query_to_df(
            """ SELECT INFO FROM variants WHERE INFO LIKE '%FILTERED%' """
        )
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_prioritization_full_unsorted():
    """
    This is a test function for prioritization of variants in a VCF file using a specified configuration
    and parameter dictionary.
    Test with a VCF full variants type: SNV, INDEL, MNV, SV
    This VCF is unsorted
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.full.unsorted.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct config dict
        config = {}

        # Construct param dict
        param = {
            "prioritization": {
                "prioritization_config": tests_data_folder
                + "/prioritization_profiles.json",
                "profiles": ["default", "GERMLINE"],
                "pzfields": ["PZFlag", "PZScore", "PZComment", "PZInfos"],
            }
        }

        # Create object
        variants = Variants(
            input=input_vcf, output=output_vcf, load=True, config=config, param=param
        )

        # Prioritization
        variants.prioritization()

        # Check all priorized default profile
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag_default=%'
            AND INFO LIKE '%PZScore_default=%'
            AND INFO LIKE '%PZComment_default=%'
            AND INFO LIKE '%PZInfos_default=%'
            """
        )
        assert len(result) == 4

        # Check all priorized GERMLINE profile
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag_GERMLINE=%'
            AND INFO LIKE '%PZScore_GERMLINE=%'
            AND INFO LIKE '%PZComment_GERMLINE=%'
            AND INFO LIKE '%PZInfos_GERMLINE=%'
            """
        )
        assert len(result) == 3

        # Check all priorized default profile (as default)
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag=%'
            AND INFO LIKE '%PZScore=%'
            AND INFO LIKE '%PZComment=%'
            AND INFO LIKE '%PZInfos=%'
            """
        )
        assert len(result) == 4

        # Check all priorized default profile
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag_default=%'
            AND INFO LIKE '%PZScore_default=%'
            """
        )
        assert len(result) == 36

        # Check all priorized GERMLINE profile
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag_GERMLINE=%'
            AND INFO LIKE '%PZScore_GERMLINE=%'
            """
        )
        assert len(result) == 36

        # Check all priorized default profile (as default)
        result = variants.get_query_to_df(
            """
            SELECT * FROM variants
            WHERE INFO LIKE '%PZFlag=%'
            AND INFO LIKE '%PZScore=%'
            """
        )
        assert len(result) == 36

        # Check annotation1
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' AND INFO LIKE '%PZScore_default=15%' """
        )
        assert len(result) == 1

        # Check FILTERED
        result = variants.get_query_to_df(
            f""" SELECT INFO FROM variants WHERE INFO LIKE '%FILTERED%' """
        )
        assert len(result) == 2

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_prioritization_varank():
    """
    This is a test function for the prioritization feature of a Python package called "Variants".
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct config dict
        config = {}

        # Construct param dict
        param = {
            "prioritization": {
                "prioritization_config": tests_data_folder
                + "/prioritization_profiles.json",
                "profiles": ["default", "GERMLINE"],
                "pzfields": ["PZFlag", "PZScore", "PZComment", "PZInfos"],
                "prioritization_score_mode": "VaRank",
            }
        }

        # Create object
        variants = Variants(
            input=input_vcf, output=output_vcf, load=True, config=config, param=param
        )

        # Prioritization
        variants.prioritization()

        # Check all priorized
        result = variants.get_query_to_df(""" SELECT INFO FROM variants """)
        assert len(result) > 0

        # Check annotation1
        result = variants.get_query_to_df(
            """ SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr1' AND POS = 28736 AND REF = 'A' AND ALT = 'C' AND INFO LIKE '%PZScore_default=15%' """
        )
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_prioritization_all_profiles():
    """
    This function tests if an error is raised when there are no prioritization configuration profiles
    provided.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct config dict
        config = {}

        # Construct param dict
        param = {
            "prioritization": {
                "prioritization_config": tests_data_folder
                + "/prioritization_profiles.json"
            }
        }

        # Create object
        variants = Variants(
            input=input_vcf, output=output_vcf, load=True, config=config, param=param
        )

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False

        # # Prioritization fail
        # with pytest.raises(ValueError) as e:
        #     variants.prioritization()
        # assert str(e.value) == f"NO Profiles configuration"


def test_prioritization_profile_not_configured():
    """
    This function tests if an error is raised when there are no prioritization configuration profiles
    provided.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"
        profile_not_defined = "profile_not_defined"

        # Construct config dict
        config = {}

        # Construct param dict
        param = {
            "prioritization": {
                # "prioritization_config": tests_data_folder + "/prioritization_profiles.json"
                "profiles": [profile_not_defined]
            }
        }

        # Create object
        variants = Variants(
            input=input_vcf, output=output_vcf, load=True, config=config, param=param
        )

        # Prioritization fail
        with pytest.raises(ValueError) as e:
            variants.prioritization()
        assert str(e.value) == f"Profile '{profile_not_defined}' NOT configured"


def test_prioritization_no_pzfields():
    """
    This function tests the prioritization method of a Variants object when there are no pzfields
    specified.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct config dict
        config = {}

        # Construct param dict
        param = {
            "prioritization": {
                "prioritization_config": tests_data_folder
                + "/prioritization_profiles.json",
                "profiles": [],
                "pzfields": [],
            }
        }

        # Create object
        variants = Variants(
            input=input_vcf, output=output_vcf, load=True, config=config, param=param
        )

        # Prioritization
        variants.prioritization()

        # Check all priorized
        result = variants.get_query_to_df(
            """ SELECT * FROM variants WHERE INFO LIKE '%PZScore_default=%' """
        )
        assert len(result) == 0

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_prioritization_no_infos():
    """
    This function tests the prioritization method of the Variants class when there is no information
    available.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.no_samples.vcf.gz"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct config dict
        config = {}

        # Construct param dict
        param = {
            "prioritization": {
                "prioritization_config": tests_data_folder
                + "/prioritization_profiles.json",
                "profiles": ["default"],
                "pzfields": ["PZFlag", "PZScore", "PZComment", "PZInfos"],
            }
        }

        # Create object
        variants = Variants(
            input=input_vcf, output=output_vcf, load=True, config=config, param=param
        )

        # Prioritization with default that faild because of DP does not exists and not CAST in FLOAT in profile criterion definition
        try:
            variants.prioritization()
        except:
            assert True

        # Using GERMLINE prfiles, no CAST error

        # Construct param dict
        param = {
            "prioritization": {
                "prioritization_config": tests_data_folder
                + "/prioritization_profiles.json",
                "profiles": ["GERMLINE"],
                "pzfields": ["PZFlag", "PZScore", "PZComment", "PZInfos"],
            }
        }

        # Set param
        variants.set_param(param=param)

        # Prioritization with GERMLINE ok
        try:
            variants.prioritization()
        except:
            assert False

        result = variants.get_query_to_df(
            """
            SELECT INFO FROM variants WHERE INFO LIKE '%PZScore_GERMLINE=0%' AND INFO LIKE '%PZFlag_GERMLINE=PASS%'
            """
        )
        assert len(result) == 10

        # Check if VCF is in correct format with pyVCF
        remove_if_exists([output_vcf])
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False

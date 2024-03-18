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





def test_annotation_exomiser():
    """
    The function `test_annotation_exomiser()` tests the annotation functionality of the Exomiser tool on
    different example files with and without HPO terms and with a pedigree.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Output file
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Number of threads
        tests_config["threads"] = 8


        ### Example file with some HPO
        # Download Exomiser databases if needed (take a while)

        # Init Param
        input_vcf = tests_data_folder + "/example.vcf.gz"
        tests_config["threads"] = 8
        param = {
            "assembly": "hg19",
            "annotation":
                {
                    "exomiser": {
                        "hpo": ['HP:0001156', '0001363', '0011304', '0010055'],
                        "transcript_source": "refseq",
                        "release": "2109"
                    }
                }
            }
        
        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, config=tests_config, param=param, load=True)

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # Export
        variants.export_output()

        # query annotated variant - check number of variants
        result = variants.get_query_to_df(""" SELECT * FROM variants """)
        log.debug(result["INFO"])
        assert len(result) == 7

        # query annotated variant - check number of annotated variants
        result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%Exomiser%' """)
        assert len(result) == 6


        ### Example file WHITHOUT HPO

        # Init Param
        input_vcf = tests_data_folder + "/example.vcf.gz"
        tests_config["threads"] = 8
        param = {
            "annotation":
                {
                    "exomiser": {
                        "hpo": []
                    }
                }
            }

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, config=tests_config, param=param, load=True)

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # Export
        variants.export_output()

        # query annotated variant - check number of variants
        result = variants.get_query_to_df(""" SELECT * FROM variants """)
        log.debug(result["INFO"])
        assert len(result) == 7

        # query annotated variant - check number of annotated variants
        result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%Exomiser%' """)
        assert len(result) == 6


        ### Example file with pedigree

        # Init Param
        input_vcf =tests_data_folder + "/example.vcf.gz"
        tests_config["threads"] = 8
        param = {
            "annotation":
                {
                    "exomiser": {
                        "phenopacket": {
                            "id": "sample1-analysis",
                            "proband": {
                                "subject": {
                                    "id": "sample1",
                                    "sex": "FEMALE"
                                },
                                "phenotypicFeatures": [
                                    {
                                        "type": {
                                            "id": "HP:0001159",
                                            "label": "Syndactyly"
                                        }
                                    },
                                    {
                                        "type": {
                                            "id": "HP:0000486",
                                            "label": "Strabismus"
                                        }
                                    },
                                    {
                                        "type": {
                                            "id": "HP:0000327",
                                            "label": "Hypoplasia of the maxilla"
                                        }
                                    },
                                    {
                                        "type": {
                                            "id": "HP:0000520",
                                            "label": "Proptosis"
                                        }
                                    },
                                    {
                                        "type": {
                                            "id": "HP:0000316",
                                            "label": "Hypertelorism"
                                        }
                                    },
                                    {
                                        "type": {
                                            "id": "HP:0000244",
                                            "label": "Brachyturricephaly"
                                        }
                                    }
                                ]
                            },
                            "pedigree": {
                                "persons": [
                                    {
                                        "individualId": "sample1",
                                        "paternalId": "sample3",
                                        "maternalId": "sample4",
                                        "sex": "FEMALE",
                                        "affectedStatus": "AFFECTED"
                                    },
                                    {
                                        "individualId": "sample2",
                                        "paternalId": "sample3",
                                        "maternalId": "sample4",
                                        "sex": "MALE",
                                        "affectedStatus": "UNAFFECTED"
                                    },
                                    {
                                        "individualId": "sample3",
                                        "sex": "MALE",
                                        "affectedStatus": "UNAFFECTED"
                                    },
                                    {
                                        "individualId": "sample4",
                                        "sex": "FEMALE",
                                        "affectedStatus": "UNAFFECTED"
                                    }
                                ]
                            }
                        }
                    }
                }
        }

        # Create object
        variants = Variants(conn=None, input=input_vcf, output=output_vcf, config=tests_config, param=param, load=True)

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # Export
        variants.export_output()

        # query annotated variant - check number of variants
        result = variants.get_query_to_df(""" SELECT * FROM variants """)
        log.debug(result["INFO"])
        assert len(result) == 7

        # query annotated variant - check number of annotated variants
        result = variants.get_query_to_df(""" SELECT * FROM variants WHERE INFO LIKE '%Exomiser%' """)
        assert len(result) == 6


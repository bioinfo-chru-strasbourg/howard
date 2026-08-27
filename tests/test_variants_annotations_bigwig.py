# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_variants_annotations_bigwig.py -x -v --log-cli-level=INFO --capture=tee-sys
coverage report --include=howard/* -m 
"""

import logging as log
import os
from tempfile import TemporaryDirectory
import vcf  # type: ignore

from howard.objects.variants import Variants
from howard.functions.commons import remove_if_exists

from test_needed import tests_folder, tests_data_folder, tests_annotations_folder


def test_annotation_bigwig_annotations_fields_and_options():
    """
    This function tests the annotation of a VCF file using bigwig annotations.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        input_vcf = tests_data_folder + "/example.full.vcf.gz"
        annotation_bigwig = os.path.join(tests_annotations_folder, "gerp.bw")
        annotation_bigwig2 = os.path.join(tests_annotations_folder, "gerp2.bw")
        annotation_bigbed_regions = os.path.join(
            tests_annotations_folder, "annotation_regions.bb"
        )
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "annotation": {
                "bigwig": {
                    "annotations": {
                        annotation_bigwig: {
                            "annotation_fields": {
                                "INFO": None
                            }
                        },
                        annotation_bigwig2: {
                            "annotation_fields": {
                                "gerp": "gerp_renamed"
                            },
                            "options": {
                                "method": {
                                    "gerp_renamed": "max",
                                }
                            }
                        },
                        annotation_bigbed_regions: {
                            "annotation_fields": {
                                "score1": "score_one",
                                "score2": "score2",
                                "strand": "strand",
                                "itemRgb": "itemRgb",
                            },
                            "options": {
                                "method": {
                                    # "score_one": "max",
                                    "score2": "min",
                                    "strand": "join",
                                    "itemRgb": "uniq",
                                }
                            }
                        },
                    },
                    "options": {
                        # "default": {
                        "method": {
                            "Integer": "mean",
                            "Float": "mean",
                            "String": "uniq",
                        }
                    },
                }
            }
        }

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant with BigWig - SNP
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;gerp=4.86;gerp_renamed=4.86'"""
        )
        assert len(result) == 1

        # query annotated variant with BigWig - MNP - gerp mean and gerp_renamed max
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'GC' AND ALT = 'TT' AND INFO = 'gerp=5.39;gerp_renamed=5.92'"""
        )
        assert len(result) == 1

        # query annotated variant with BigBed
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr9' AND POS = 37020622 AND INFO LIKE '%score_one=878522;score2=176%'"""
        )
        assert len(result) == 1

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False


def test_annotation_bigwig_annotations():
    """
    This function tests the annotation of a VCF file using bigwig annotations.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        input_vcf = tests_data_folder + "/example.full.vcf.gz"
        annotation_bigwig = os.path.join(tests_annotations_folder, "gerp.bw")
        annotation_bigwig2 = os.path.join(tests_annotations_folder, "gerp2.bw")
        annotation_bigbed_regions = os.path.join(
            tests_annotations_folder, "annotation_regions.bb"
        )
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "annotation": {
                "bigwig": {
                    "annotations": {
                        annotation_bigwig: {
                            "annotation_fields": {
                                "INFO": None
                            }
                        },
                        annotation_bigwig2: {
                            "annotation_fields": {
                                "gerp": "gerp_renamed"
                            }
                        },
                        annotation_bigbed_regions: {
                            "annotation_fields": {
                                "score1": "score_one",
                                "score2": "score2",
                                "strand": "strand",
                                "itemRgb": "itemRgb",
                            }
                        },
                    }
                }
            }
        }

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant with BigWig - SNP
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;gerp=4.86;gerp_renamed=4.86'"""
        )
        assert len(result) == 1, f"Query annotated variant with BigWig - SNP - Expected 1, but got {len(result)}. Result: {result.to_string()}"

        # query annotated variant with BigWig - MNP - no method specified, so gerp and gerp_renamed should be the same
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'GC' AND ALT = 'TT' AND INFO = 'gerp=5.39;gerp_renamed=5.39'"""
        )
        assert len(result) == 1, f"Query annotated variant with BigWig - MNP - no method specified, so gerp and gerp_renamed should be the same - Expected 1, but got {len(result)}. Result: {result.to_string()}"

        # query annotated variant with BigBed
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr9' AND POS = 37020622 AND INFO LIKE '%score_one=878522;score2=176%'"""
        )
        assert len(result) == 1, f"Query annotated variant with BigBed - Expected 1, but got {len(result)}. Result: {result.to_string()}"

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False, f"VCF file {output_vcf} is not in correct format. Please check the output VCF file."



def test_annotation_bigwig_annotations_renamed_header():
    """
    This function tests the annotation of a VCF file using bigwig annotations.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        input_vcf = tests_data_folder + "/example.full.vcf.gz"
        annotation_bigwig = os.path.join(tests_annotations_folder, "gerp.bw")
        annotation_bigwig2 = os.path.join(tests_annotations_folder, "gerp2.bw")
        annotation_bigbed_regions = os.path.join(
            tests_annotations_folder, "annotation_regions.bb"
        )
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "annotation": {
                "bigwig": {
                    "annotations": {
                        annotation_bigwig: {
                            "annotation_fields": {
                                "INFO": None
                            }
                        },
                        annotation_bigwig2: {
                            "annotation_fields": {
                                "gerp": "gerp_renamed"
                            }
                        },
                        annotation_bigbed_regions: {
                            "annotation_fields": {
                                "score1": "score_one",
                                "score2": "score2",
                                "strand": "strand",
                                "itemRgb": "itemRgb",
                            },
                            "options": {
                                "header_fields": {
                                    "score_one": {
                                        "Number": "1",
                                        "Type": "String",
                                        "Description": "Score one renamed field"
                                    }
                                }
                            }
                        },
                    }
                }
            }
        }

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant with BigWig - SNP
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;gerp=4.86;gerp_renamed=4.86'"""
        )
        assert len(result) == 1, f"Query annotated variant with BigWig - SNP - Expected 1, but got {len(result)}. Result: {result.to_string()}"

        # query annotated variant with BigWig - MNP - no method specified, so gerp and gerp_renamed should be the same
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'GC' AND ALT = 'TT' AND INFO = 'gerp=5.39;gerp_renamed=5.39'"""
        )
        assert len(result) == 1, f"Query annotated variant with BigWig - MNP - no method specified, so gerp and gerp_renamed should be the same - Expected 1, but got {len(result)}. Result: {result.to_string()}"

        # query annotated variant with BigBed
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr9' AND POS = 37020622 AND INFO LIKE '%score_one=878522;score2=176%'"""
        )
        assert len(result) == 1, f"Query annotated variant with BigBed - Expected 1, but got {len(result)}. Result: {result.to_string()}"

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False, f"VCF file {output_vcf} is not in correct format. Please check the output VCF file."

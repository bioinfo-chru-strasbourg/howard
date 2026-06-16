# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_objects_variants.py -x -v --log-cli-level=INFO --capture=tee-sys
coverage report --include=howard/* -m 
"""

import os
import logging as log
from tempfile import TemporaryDirectory
import vcf  # type: ignore

from howard.objects.variants import Variants
from howard.functions.commons import remove_if_exists

from test_needed import tests_data_folder, tests_folder, tests_config, tests_annotations_folder



def test_annotation_annovar():
    """
    This function tests the annotation of variants using Annovar and checks if the output VCF file is in
    the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        annotation_annovar = "nci60"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "annotation": {
                "annovar": {"annotations": {annotation_annovar: {"annotation_fields": {"INFO": None}, "options": {"operation": "f"}}}}
            }
        }

        # Create object
        variants = Variants(
            conn=None,
            input=input_vcf,
            output=output_vcf,
            config=tests_config,
            param=param,
            load=True,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci60=0.66'"""
        )
        assert len(result) == 1, f"Check annotated variant with nci60, Expected 1 annotated variant, but got {len(result)}"

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False, "VCF file is not in correct format"


def test_annotation_annovar_complete_format():
    """
    This function tests the annotation of variants using Annovar and checks if the output VCF file is in
    the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        annotation_annovar = "NCI60"
        code_annovar = "nci60"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "annotation": {
                "annovar": {
                    "annotations": {
                        annotation_annovar: {
                            "annotation_fields": {"nci60": "NCI60"},
                            "options": {
                                "protocol": code_annovar,
                                "operation": "f",
                                "genebase": "",
                                "arguments": "",
                                "xref": "",
                                "options": ""
                            }
                        }
                    }
                }
            }
        }

        # Create object
        variants = Variants(
            conn=None,
            input=input_vcf,
            output=output_vcf,
            config=tests_config,
            param=param,
            load=True,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;NCI60=0.66'"""
        )
        assert len(result) == 1, f"Check annotated variant with NCI60, Expected 1 annotated variant, but got {len(result)}"

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False, "VCF file is not in correct format"


def test_annotation_annovar_complete_format_refgene():
    """
    This function tests the annotation of variants using Annovar and checks if the output VCF file is in
    the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        annotation_annovar = "NCI60"
        code_annovar = "nci60"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "annotation": {
                "annovar": {
                    "annotations": {
                        "REFGENE": {
                            "annotation_fields": {
                                "Func_refGene": "location",
                                "ExonicFunc_refGene": "outcome",
                                "Gene_refGene": "symbol",
                                "AAChange_refGene": "hgvs"
                            },
                            "options": {
                                "protocol": "refGene",
                                "operation": "g",
                                "genebase": "-splicing_threshold 3 -indel_splicing_threshold 3",
                                "arguments": "-hgvs",
                                "xref": "",
                                "options": ""
                            }
                            }
                    }
                }
            }
        }

        # Create object
        variants = Variants(
            conn=None,
            input=input_vcf,
            output=output_vcf,
            config=tests_config,
            param=param,
            load=True,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query full variant
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE INFO LIKE '%location=%'"""
        )
        assert len(result) == 7, f"Check full variant annotated with location, Expected 7 annotated variants, but got {len(result)}"

        # query annotated variant
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO LIKE 'DP=125;%symbol=EGFR%'"""
        )
        assert len(result) == 1, f"Check annotated variant with symbol, Expected 1 annotated variant, but got {len(result)}"

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False, "VCF file is not in correct format"



def test_annotation_annovar_complete_format_refgene_xref():
    """
    This function tests the annotation of variants using Annovar and checks if the output VCF file is in
    the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = os.path.join(tests_data_folder, "example.vcf.gz")
        annotation_annovar = "REFGENE"
        code_annovar = "refGene"
        xref = os.path.join(tests_annotations_folder, "gene_xref.txt")
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "annotation": {
                "annovar": {
                    "annotations": {
                        annotation_annovar: {
                            "annotation_fields": {
                                "Func_refGene": "location",
                                "ExonicFunc_refGene": "outcome",
                                "Gene_refGene": "symbol",
                                "AAChange_refGene": "hgvs",
                                "pLi_refGene": "pLiXRef",
                                "Gene_Full_Name_refGene": "geneFullNameXRef"
                            },
                            "options": {
                                "protocol": code_annovar,
                                "operation": "gx",
                                "genebase": "-splicing_threshold 3 -indel_splicing_threshold 3",
                                "arguments": "-hgvs",
                                "xref": xref,
                                "options": ""
                            }
                            }
                    }
                }
            }
        }

        # Create object
        variants = Variants(
            conn=None,
            input=input_vcf,
            output=output_vcf,
            config=tests_config,
            param=param,
            load=True,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query full variant
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE INFO LIKE '%location=%'"""
        )
        assert len(result) == 7, f"Check full variant annotated with location, Expected 7 annotated variants, but got {len(result)}"

        # query annotated variant
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO LIKE 'DP=125;%symbol=EGFR%geneFullNameXRef=alpha-1-B_glycoprotein%'"""
        )
        assert len(result) == 1, f"Check annotated variant with symbol and geneFullNameXRef, Expected 1 annotated variant, but got {len(result)}"

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False, "VCF file is not in correct format"


def test_annotation_annovar_previous_format():
    """
    This function tests the annotation of variants using Annovar and checks if the output VCF file is in
    the correct format.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        annotation_annovar = "nci60"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "annotation": {
                "annovar": {"annotations": {annotation_annovar: {"INFO": None}}}
            }
        }

        # Create object
        variants = Variants(
            conn=None,
            input=input_vcf,
            output=output_vcf,
            config=tests_config,
            param=param,
            load=True,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci60=0.66'"""
        )
        assert len(result) == 1, f"Check annotated variant with nci60, Expected 1 annotated variant, but got {len(result)}"

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False, "VCF file is not in correct format"



def test_annotation_annovar_full_unsorted_previous_format():
    """
    The function tests the annotation of variants using Annovar and checks if the output VCF file is in
    the correct format.
    Test with a VCF full variants type: SNV, INDEL, MNV, SV
    This VCF is unsorted
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.full.unsorted.vcf.gz"
        annotation_annovar = "nci60"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "annotation": {
                "annovar": {"annotations": {annotation_annovar: {"INFO": None}}}
            }
        }

        # Create object
        variants = Variants(
            conn=None,
            input=input_vcf,
            output=output_vcf,
            config=tests_config,
            param=param,
            load=True,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci60=0.66'"""
        )
        assert len(result) == 1, f"Check annotated variant with nci60, Expected 1 annotated variant, but got {len(result)}"

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False, "VCF file is not in correct format"


def test_annotation_annovar_no_samples_previous_format():
    """
    This function tests the annotation of a VCF file using Annovar when there are no samples present.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.no_samples.vcf.gz"
        annotation_annovar = "nci60"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {
            "annotation": {
                "annovar": {"annotations": {annotation_annovar: {"INFO": None}}}
            }
        }

        # Create object
        variants = Variants(
            conn=None,
            input=input_vcf,
            output=output_vcf,
            config=tests_config,
            param=param,
            load=True,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df(
            """ SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr12' AND POS = 68724951 AND REF = 'G' AND ALT = 'T' AND INFO = 'nci60=0.77' """
        )
        assert len(result) == 1, f"Check annotated variant with nci60, Expected 1 annotated variant, but got {len(result)}"

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False, "VCF file is not in correct format"


def test_annotation_annovar_sqlite_previous_format():
    """
    This function tests the annotation of variants using Annovar and SQLite database.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        annotation_annovar = "nci60"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct config dict
        config = tests_config.copy()
        config["connexion_format"] = "sqlite"

        # Construct param dict
        param = {
            "annotation": {
                "annovar": {"annotations": {annotation_annovar: {"INFO": None}}}
            }
        }

        # Create object
        variants = Variants(
            conn=None,
            input=input_vcf,
            output=output_vcf,
            config=config,
            param=param,
            load=True,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci60=0.66'"""
        )
        assert len(result) == 1, f"Check annotated variant with nci60, Expected 1 annotated variant, but got {len(result)}"

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False, "VCF file is not in correct format"


def test_annotation_quick_annovar_previous_format():
    """
    This function tests the annotation of a VCF file using Annovar.
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        input_vcf = tests_data_folder + "/example.vcf.gz"
        annotation_annovar = "nci60"
        output_vcf = f"{tmp_dir}/output.vcf.gz"

        # Construct param dict
        param = {"annotations": f"annovar:{annotation_annovar}"}

        # Create object
        variants = Variants(
            conn=None,
            input=input_vcf,
            output=output_vcf,
            config=tests_config,
            param=param,
            load=True,
        )

        # Remove if output file exists
        remove_if_exists([output_vcf])

        # Annotation
        variants.annotation()

        # query annotated variant
        result = variants.get_query_to_df(
            """SELECT 1 AS count FROM variants WHERE "#CHROM" = 'chr7' AND POS = 55249063 AND REF = 'G' AND ALT = 'A' AND INFO = 'DP=125;nci60=0.66'"""
        )
        assert len(result) == 1, f"Check annotated variant with nci60, Expected 1 annotated variant, but got {len(result)}"

        # Check if VCF is in correct format with pyVCF
        variants.export_output()
        try:
            vcf.Reader(filename=output_vcf)
        except:
            assert False, "VCF file is not in correct format"

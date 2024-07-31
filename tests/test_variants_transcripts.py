# -*- coding: utf-8 -*-
"""
Tests

Usage:
pytest tests/

Coverage:
coverage run -m pytest tests/test_variants_transcripts.py -x -vv --log-cli-level=DEBUG --capture=tee-sys
coverage report --include=howard/* -m
"""

import logging as log
from tempfile import TemporaryDirectory
import pytest

from howard.objects.variants import Variants
from test_needed import tests_folder, tests_config


@pytest.mark.parametrize(
    "input_vcf",
    [
        "tests/data/example.ann.transcripts.vcf.gz",
        "tests/data/example.ann.vcf.gz",
        "tests/data/example.dbnsfp.transcripts.vcf.gz",
        "tests/data/example.dbnsfp.no_transcripts.vcf.gz",
    ],
)
def test_devel_create_transcript_view(input_vcf):
    """
    The function `test_devel_create_transcript_view` creates a transcript view from a VCF file using
    specified parameters and checks the resulting table for data.

    :param input_vcf: It seems like the `input_vcf` parameter is missing in the provided code snippet.
    Could you please provide the value or path that should be assigned to the `input_vcf` variable in
    the `test_devel_create_transcript_view` function?
    """

    with TemporaryDirectory(dir=tests_folder) as tmp_dir:

        # Init files
        output_vcf = f"{tmp_dir}/output.vcf"

        # Construct param dict
        param = {
            "transcripts": {
                "table": "transcripts",
                "struct": {
                    "from_column_format": [  # format List, e.g. snpEff
                        {
                            "transcripts_column": "ANN",
                            "transcripts_infos_column": "Feature_ID",
                        }
                    ],
                    "from_columns_map": [  # format List, e.g. dbNSFP columns
                        {
                            "transcripts_column": "Ensembl_transcriptid",
                            "transcripts_infos_columns": [
                                "genename",
                                "Ensembl_geneid",
                                "LIST_S2_score",
                                "LIST_S2_pred",
                            ],
                        },
                        {
                            "transcripts_column": "Ensembl_transcriptid",
                            "transcripts_infos_columns": [
                                "genename",
                                "VARITY_R_score",
                                "Aloft_pred",
                            ],
                        },
                    ],
                },
            }
        }

        # Create object
        variants = Variants(
            conn=None, input=input_vcf, output=output_vcf, param=param, load=True
        )

        # Create transcript view
        transcripts_table = variants.create_transcript_view()

        # Check table exists
        assert transcripts_table is not None

        # Check table content
        query_check = f"""
            SELECT * FROM {transcripts_table}
            ORDER BY "#CHROM", POS, REF, ALT, transcript
        """
        check = variants.get_query_to_df(query=query_check)
        log.debug(f"check={check}")
        assert len(check) > 0

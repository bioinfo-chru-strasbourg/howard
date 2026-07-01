"""
Tests for annotation header fields override feature:
- variants_annotation.get_annotation_header_fields_override
- variants_annotation.build_info_with_header_override
- annotation_parquet "reheader" logic (natural Number/Type/Description
  fallback computed from the parquet DB's own header, then overridden)

Usage:
pytest tests/test_variants_annotations_header_override.py -x -v --log-cli-level=DEBUG --capture=tee-sys

Coverage:
coverage run -m pytest tests/test_variants_annotations_header_override.py -x -v
coverage report --include=howard/* -m
"""
import pytest
from howard.objects.variants import Variants


def _create_minimal_vcf(tmp_path) -> str:
    """
    Create a minimal header-only VCF file, sufficient to instantiate
    a Variants object (same pattern as Variants(input=db_hdr_file)
    used elsewhere in annotation.py for header-only parsing).
    """
    vcf_content = (
        "##fileformat=VCFv4.2\n"
        "##contig=<ID=chr1,length=249250621>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )
    vcf_file = tmp_path / "minimal.vcf"
    vcf_file.write_text(vcf_content)
    return str(vcf_file)


def _get_variants_object(tmp_path, param: dict = None) -> Variants:
    """
    Helper to build a minimal Variants object, optionally with param set.
    """
    input_vcf = _create_minimal_vcf(tmp_path)
    variants = Variants(input=input_vcf)
    if param is not None:
        variants.set_param(param)
    return variants


###############################################################################
# get_annotation_header_fields_override
###############################################################################


def test_get_annotation_header_fields_override_default(tmp_path):
    """
    No param set -> empty dict
    """
    variants = _get_variants_object(tmp_path, param={})
    assert variants.get_annotation_header_fields_override() == {}


def test_get_annotation_header_fields_override_missing_sections(tmp_path):
    """
    Partial param structure (missing 'options' or 'header_fields') -> empty dict
    """
    variants = _get_variants_object(tmp_path, param={"annotation": {}})
    assert variants.get_annotation_header_fields_override() == {}

    variants = _get_variants_object(tmp_path, param={"annotation": {"options": {}}})
    assert variants.get_annotation_header_fields_override() == {}


def test_get_annotation_header_fields_override_full(tmp_path):
    """
    Full 'header_fields' override dict provided -> returned as-is
    """
    header_fields = {
        "my_field": {
            "number": ".",
            "type": "Float",
            "description": "my field description",
        },
        "another_field": {"type": "Integer"},
    }
    param = {"annotation": {"options": {"header_fields": header_fields}}}
    variants = _get_variants_object(tmp_path, param=param)
    assert variants.get_annotation_header_fields_override() == header_fields


###############################################################################
# build_info_with_header_override
###############################################################################


def test_build_info_no_override_none(tmp_path):
    """
    header_fields_override=None -> original values kept unchanged
    """
    variants = _get_variants_object(tmp_path)
    info = variants.build_info_with_header_override(
        field_name="my_field",
        field_number=".",
        field_type="String",
        field_description="my description",
        field_source="my_source",
        field_version="my_version",
        header_fields_override=None,
    )
    assert info.id == "my_field"
    assert info.num == "."
    assert info.type == "String"
    assert info.desc == "my description"
    assert info.source == "my_source"
    assert info.version == "my_version"
    assert info.type_code == variants.code_type_map["String"]


def test_build_info_no_override_empty_dict(tmp_path):
    """
    header_fields_override={} -> original values kept unchanged
    """
    variants = _get_variants_object(tmp_path)
    info = variants.build_info_with_header_override(
        field_name="my_field",
        field_number="1",
        field_type="Integer",
        field_description="my description",
        field_source="my_source",
        field_version="my_version",
        header_fields_override={},
    )
    assert info.num == "1"
    assert info.type == "Integer"
    assert info.desc == "my description"


def test_build_info_override_field_not_matched(tmp_path):
    """
    Override dict provided but does NOT contain this field -> original kept
    """
    variants = _get_variants_object(tmp_path)
    info = variants.build_info_with_header_override(
        field_name="my_field",
        field_number=".",
        field_type="String",
        field_description="my description",
        field_source="my_source",
        field_version="my_version",
        header_fields_override={"other_field": {"type": "Float"}},
    )
    assert info.type == "String"
    assert info.num == "."
    assert info.desc == "my description"


def test_build_info_full_override(tmp_path):
    """
    Override provides number, type and description -> all overridden
    """
    variants = _get_variants_object(tmp_path)
    header_fields_override = {
        "my_field": {
            "number": "1",
            "type": "Float",
            "description": "overridden description",
        }
    }
    info = variants.build_info_with_header_override(
        field_name="my_field",
        field_number=".",
        field_type="String",
        field_description="original description",
        field_source="my_source",
        field_version="my_version",
        header_fields_override=header_fields_override,
    )
    assert info.num == "1"
    assert info.type == "Float"
    assert info.desc == "overridden description"
    assert info.type_code == variants.code_type_map["Float"]


def test_build_info_partial_override_type_only(tmp_path):
    """
    Override provides only 'type' -> number/description untouched
    """
    variants = _get_variants_object(tmp_path)
    info = variants.build_info_with_header_override(
        field_name="my_field",
        field_number=".",
        field_type="String",
        field_description="original description",
        field_source="my_source",
        field_version="my_version",
        header_fields_override={"my_field": {"type": "Integer"}},
    )
    assert info.type == "Integer"
    assert info.num == "."
    assert info.desc == "original description"


def test_build_info_partial_override_number_only(tmp_path):
    """
    Override provides only 'number' -> type/description untouched
    """
    variants = _get_variants_object(tmp_path)
    info = variants.build_info_with_header_override(
        field_name="my_field",
        field_number=".",
        field_type="String",
        field_description="original description",
        field_source="my_source",
        field_version="my_version",
        header_fields_override={"my_field": {"number": "A"}},
    )
    assert info.num == "A"
    assert info.type == "String"


def test_build_info_partial_override_description_only(tmp_path):
    """
    Override provides only 'description' -> number/type untouched.
    This is the exact pattern used in production param files
    (see param.DIAGGEN.GOMV1_GERMLINE.json) where only 'description' is set.
    """
    variants = _get_variants_object(tmp_path)
    info = variants.build_info_with_header_override(
        field_name="my_field",
        field_number=".",
        field_type="String",
        field_description="original description",
        field_source="my_source",
        field_version="my_version",
        header_fields_override={"my_field": {"description": "new description"}},
    )
    assert info.desc == "new description"
    assert info.num == "."
    assert info.type == "String"


def test_build_info_invalid_type_raises(tmp_path):
    """
    Override 'type' not in CODE_TYPE_MAP -> raises ValueError
    """
    variants = _get_variants_object(tmp_path)
    with pytest.raises(ValueError):
        variants.build_info_with_header_override(
            field_name="my_field",
            field_number=".",
            field_type="String",
            field_description="description",
            field_source="source",
            field_version="version",
            header_fields_override={"my_field": {"type": "NotAValidType"}},
        )


def test_build_info_multiple_fields_isolation(tmp_path):
    """
    Override dict with multiple fields -> only the matching field is impacted
    """
    variants = _get_variants_object(tmp_path)
    header_fields_override = {
        "field_a": {"type": "Float"},
        "field_b": {"type": "Integer"},
    }

    info_a = variants.build_info_with_header_override(
        field_name="field_a",
        field_number=".",
        field_type="String",
        field_description="desc a",
        field_source="src",
        field_version="ver",
        header_fields_override=header_fields_override,
    )
    info_b = variants.build_info_with_header_override(
        field_name="field_b",
        field_number=".",
        field_type="String",
        field_description="desc b",
        field_source="src",
        field_version="ver",
        header_fields_override=header_fields_override,
    )
    info_c = variants.build_info_with_header_override(
        field_name="field_c",
        field_number=".",
        field_type="String",
        field_description="desc c",
        field_source="src",
        field_version="ver",
        header_fields_override=header_fields_override,
    )

    assert info_a.type == "Float"
    assert info_b.type == "Integer"
    assert info_c.type == "String"


###############################################################################
# Parquet "reheader" logic (annotation_parquet)
#
# annotation_parquet computes NATURAL Number/Type/Description values from the
# parquet database's own header (parquet_hdr_vcf_header_infos[annotation_field])
# BEFORE calling build_info_with_header_override(), following this fallback:
#
#   if parquet_type in ["regions"]:
#       number = "."
#   else:
#       number = source_info.num or "."
#   type_ = source_info.type or "String"
#   description = source_info.desc or f"{annotation_field} description"
#
# The tests below replicate this natural computation with a lightweight stub
# for parquet_hdr_vcf_header_infos[annotation_field] (only .num/.type/.desc
# are read by that logic), then feed the result into
# build_info_with_header_override(), exactly as annotation_parquet does,
# using field_name=annotation_fields_new_name (the RENAMED field), which is
# the key used to look up the global header_fields override.
###############################################################################


class _FakeParquetSourceInfo:
    """
    Minimal stand-in for parquet_hdr_vcf_header_infos[annotation_field]
    (only .num/.type/.desc are read by annotation_parquet's natural
    fallback computation).
    """

    def __init__(self, num=None, type_=None, desc=None):
        self.num = num
        self.type = type_
        self.desc = desc


def _parquet_natural_header_values(
    source_info: _FakeParquetSourceInfo,
    annotation_field: str,
    parquet_type: str = "variants",
):
    """
    Replicates annotation_parquet's natural (pre-override) computation of
    Number/Type/Description for a parquet-sourced INFO field.
    """
    if parquet_type in ["regions"]:
        number = "."
    else:
        number = source_info.num or "."
    type_ = source_info.type or "String"
    description = source_info.desc or f"{annotation_field} description"
    return number, type_, description


def test_parquet_reheader_no_override_variants_type(tmp_path):
    """
    parquet_type='variants' (normal), no override -> natural values kept as-is
    """
    variants = _get_variants_object(tmp_path, param={})
    source_info = _FakeParquetSourceInfo(num="1", type_="Float", desc="orig desc")

    number, type_, description = _parquet_natural_header_values(
        source_info, annotation_field="my_field", parquet_type="variants"
    )
    info = variants.build_info_with_header_override(
        field_name="my_field",
        field_number=number,
        field_type=type_,
        field_description=description,
        field_source="parquet_db",
        field_version="unknown",
        header_fields_override=variants.get_annotation_header_fields_override(),
    )

    assert info.num == "1"
    assert info.type == "Float"
    assert info.desc == "orig desc"


def test_parquet_reheader_regions_forces_number_dot(tmp_path):
    """
    parquet_type='regions' -> Number is forced to '.' regardless of source
    Number, when no override is set
    """
    variants = _get_variants_object(tmp_path, param={})
    source_info = _FakeParquetSourceInfo(num="1", type_="Float", desc="orig desc")

    number, type_, description = _parquet_natural_header_values(
        source_info, annotation_field="my_field", parquet_type="regions"
    )
    info = variants.build_info_with_header_override(
        field_name="my_field",
        field_number=number,
        field_type=type_,
        field_description=description,
        field_source="parquet_db",
        field_version="unknown",
        header_fields_override=variants.get_annotation_header_fields_override(),
    )

    assert info.num == "."
    assert info.type == "Float"


def test_parquet_reheader_missing_type_and_desc_fallback(tmp_path):
    """
    Source header has no Type/Description -> fallback to 'String' and
    '{field} description'
    """
    variants = _get_variants_object(tmp_path, param={})
    source_info = _FakeParquetSourceInfo(num=None, type_=None, desc=None)

    number, type_, description = _parquet_natural_header_values(
        source_info, annotation_field="my_field", parquet_type="variants"
    )
    info = variants.build_info_with_header_override(
        field_name="my_field",
        field_number=number,
        field_type=type_,
        field_description=description,
        field_source="parquet_db",
        field_version="unknown",
        header_fields_override=variants.get_annotation_header_fields_override(),
    )

    assert info.num == "."
    assert info.type == "String"
    assert info.desc == "my_field description"


def test_parquet_reheader_regions_override_number_wins(tmp_path):
    """
    parquet_type='regions' forces Number='.', but an explicit override
    for Number must still take final precedence
    """
    param = {
        "annotation": {"options": {"header_fields": {"my_field": {"number": "A"}}}}
    }
    variants = _get_variants_object(tmp_path, param=param)
    source_info = _FakeParquetSourceInfo(num="1", type_="Float", desc="orig desc")

    number, type_, description = _parquet_natural_header_values(
        source_info, annotation_field="my_field", parquet_type="regions"
    )
    info = variants.build_info_with_header_override(
        field_name="my_field",
        field_number=number,
        field_type=type_,
        field_description=description,
        field_source="parquet_db",
        field_version="unknown",
        header_fields_override=variants.get_annotation_header_fields_override(),
    )

    assert info.num == "A"  # override wins over the "regions -> '.'" rule


def test_parquet_reheader_rename_then_override_by_new_name(tmp_path):
    """
    Critical regression test: when a parquet field is renamed via
    'annotation_fields' (e.g. ALLELEFREQ -> dejavu.service.application),
    the global header_fields override must be looked up using the RENAMED
    (final) field name - matching real production param files - not the
    original database column name.
    """
    raw_field_name = "ALLELEFREQ"
    renamed_field_name = "dejavu.service.application"

    param = {
        "annotation": {
            "options": {
                "header_fields": {
                    renamed_field_name: {"description": "overridden by new name"}
                }
            }
        }
    }
    variants = _get_variants_object(tmp_path, param=param)
    source_info = _FakeParquetSourceInfo(num="1", type_="Float", desc="orig desc")

    number, type_, description = _parquet_natural_header_values(
        source_info, annotation_field=raw_field_name, parquet_type="variants"
    )

    # Override IS applied: lookup done with the renamed field name
    info_renamed = variants.build_info_with_header_override(
        field_name=renamed_field_name,
        field_number=number,
        field_type=type_,
        field_description=description,
        field_source="parquet_db",
        field_version="unknown",
        header_fields_override=variants.get_annotation_header_fields_override(),
    )
    assert info_renamed.desc == "overridden by new name"
    assert info_renamed.num == "1"
    assert info_renamed.type == "Float"

    # Override is NOT applied: lookup done with the original (raw) field name
    info_raw = variants.build_info_with_header_override(
        field_name=raw_field_name,
        field_number=number,
        field_type=type_,
        field_description=description,
        field_source="parquet_db",
        field_version="unknown",
        header_fields_override=variants.get_annotation_header_fields_override(),
    )
    assert info_raw.desc == "orig desc"


def test_parquet_reheader_real_diaggen_param_description_only(tmp_path):
    """
    Reproduces the real production param.DIAGGEN.GOMV1_GERMLINE.json use case:
    only 'description' is overridden for parquet-sourced fields; Number/Type
    natural values from the parquet DB's own header must be preserved.
    """
    param = {
        "annotation": {
            "options": {
                "header_fields": {
                    "dejavu.service.application": {
                        "description": "dejavu generated by vAnnot for GOMV1_GERMLINE"
                    }
                }
            }
        }
    }
    variants = _get_variants_object(tmp_path, param=param)

    # Simulate parquet DB's own header: ALLELEFREQ, Number=1, Type=Float
    source_info = _FakeParquetSourceInfo(
        num="1", type_="Float", desc="Allele frequency"
    )

    number, type_, description = _parquet_natural_header_values(
        source_info, annotation_field="ALLELEFREQ", parquet_type="variants"
    )
    info = variants.build_info_with_header_override(
        field_name="dejavu.service.application",
        field_number=number,
        field_type=type_,
        field_description=description,
        field_source="parquet_db",
        field_version="unknown",
        header_fields_override=variants.get_annotation_header_fields_override(),
    )

    assert info.id == "dejavu.service.application"
    assert info.num == "1"  # untouched, from parquet DB natural header
    assert info.type == "Float"  # untouched, from parquet DB natural header
    assert info.desc == "dejavu generated by vAnnot for GOMV1_GERMLINE"  # overridden
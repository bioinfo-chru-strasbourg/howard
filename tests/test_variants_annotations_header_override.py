"""
Tests for annotation header fields override feature:
- variants_annotation.get_annotation_header_fields_override(annotations)
- variants_annotation.build_info_with_header_override(...)

Design under test:
- get_annotation_header_fields_override(annotations) collects Number/Type/
  Description override dicts ONLY from each annotation's own "options"
  block (per-annotation level, sibling of "annotation_fields"):

      "annotation": {
          "annovar": {
              "annotations": {
                  "ALL.sites.2015_08": {
                      "annotation_fields": {
                          "ALL.sites.2015_08": "1000genomesALL"
                      },
                      "options": {
                          "header_fields": {
                              "1000genomesALL": {
                                  "number": ".",
                                  "type": "Float",
                                  "description": "1000genomesALL by Annovar"
                              }
                          }
                      }
                  }
              }
          }
      }

  There is NO tool-level/global header_fields config considered - only
  the per-annotation "options.header_fields" sub-key is read. If several
  annotations define an override for the same final field name, the LAST
  one encountered (dict iteration/insertion order) wins, since
  dict.update() is used internally.

- build_info_with_header_override(...) applies Number/Type/Description
  override (if any) on top of natural/source values, and validates 'type'
  against CODE_TYPE_MAP.

- Currently, only annotation_annovar() calls
  get_annotation_header_fields_override(annotations=annotations) and
  build_info_with_header_override(...), at the point where the final
  merged/multianno ANNOVAR VCF header is re-parsed.

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
# get_annotation_header_fields_override(annotations)
###############################################################################


def test_get_annotation_header_fields_override_empty_annotations(tmp_path):
    """
    Empty annotations dict -> empty override dict
    """
    variants = _get_variants_object(tmp_path)
    assert variants.get_annotation_header_fields_override(annotations={}) == {}


def test_get_annotation_header_fields_override_no_options_key(tmp_path):
    """
    Annotation entry with NO 'options' key at all (flat rename dict,
    e.g. "avsnp150": {"avsnp150": "dbSNP"}) -> no override collected,
    no error
    """
    variants = _get_variants_object(tmp_path)
    annotations = {"avsnp150": {"avsnp150": "dbSNP"}}
    assert (
        variants.get_annotation_header_fields_override(annotations=annotations) == {}
    )


def test_get_annotation_header_fields_override_empty_annotation_entry(tmp_path):
    """
    Annotation entry as an empty dict (e.g. "refGene": {}) -> no override
    collected, no error
    """
    variants = _get_variants_object(tmp_path)
    annotations = {"refGene": {}}
    assert (
        variants.get_annotation_header_fields_override(annotations=annotations) == {}
    )


def test_get_annotation_header_fields_override_options_without_header_fields(
    tmp_path,
):
    """
    Annotation entry has an 'options' block, but WITHOUT 'header_fields'
    key (e.g. only 'operation'/'parallelize') -> no override collected
    """
    variants = _get_variants_object(tmp_path)
    annotations = {
        "refGene": {"options": {"operation": "g"}},
    }
    assert (
        variants.get_annotation_header_fields_override(annotations=annotations) == {}
    )


def test_get_annotation_header_fields_override_single_annotation(tmp_path):
    """
    Single annotation defines 'options.header_fields' -> collected as-is
    """
    variants = _get_variants_object(tmp_path)
    annotations = {
        "ALL.sites.2015_08": {
            "annotation_fields": {"ALL.sites.2015_08": "1000genomesALL"},
            "options": {
                "header_fields": {
                    "1000genomesALL": {
                        "number": ".",
                        "type": "Float",
                        "description": "1000genomesALL by Annovar",
                    }
                }
            },
        }
    }
    override = variants.get_annotation_header_fields_override(annotations=annotations)
    assert override == {
        "1000genomesALL": {
            "number": ".",
            "type": "Float",
            "description": "1000genomesALL by Annovar",
        }
    }


def test_get_annotation_header_fields_override_multiple_annotations_merged(tmp_path):
    """
    Multiple annotations, each defining header_fields for DIFFERENT field
    names -> all merged into a single dict
    """
    variants = _get_variants_object(tmp_path)
    annotations = {
        "ALL.sites.2015_08": {
            "annotation_fields": {"ALL.sites.2015_08": "1000genomesALL"},
            "options": {"header_fields": {"1000genomesALL": {"type": "Float"}}},
        },
        "dbnsfp41a": {
            "annotation_fields": {"CADD_raw": "CADD_raw"},
            "options": {"header_fields": {"CADD_raw": {"type": "Float"}}},
        },
        # No options at all - must be ignored, not raise
        "refGene": {},
    }
    override = variants.get_annotation_header_fields_override(annotations=annotations)
    assert override == {
        "1000genomesALL": {"type": "Float"},
        "CADD_raw": {"type": "Float"},
    }


def test_get_annotation_header_fields_override_only_annotations_with_header_fields(
    tmp_path,
):
    """
    Mix of annotations with and without header_fields -> only those WITH
    header_fields contribute to the result
    """
    variants = _get_variants_object(tmp_path)
    annotations = {
        "avsnp150": {"avsnp150": "dbSNP"},
        "refGene": {},
        "ALL.sites.2015_08": {
            "annotation_fields": {"ALL.sites.2015_08": "1000genomesALL"},
            "options": {
                "header_fields": {"1000genomesALL": {"description": "custom desc"}}
            },
        },
        "clinvar_20210123": {"CLNSIG": "CLINVAR"},
    }
    override = variants.get_annotation_header_fields_override(annotations=annotations)
    assert override == {"1000genomesALL": {"description": "custom desc"}}


def test_get_annotation_header_fields_override_same_field_last_wins(tmp_path):
    """
    Two DIFFERENT annotation entries both define header_fields for the SAME
    final field name -> last one encountered (dict iteration/insertion
    order) wins, since dict.update() is used internally
    """
    variants = _get_variants_object(tmp_path)
    annotations = {
        "annotation_a": {
            "options": {"header_fields": {"shared_field": {"type": "Integer"}}}
        },
        "annotation_b": {
            "options": {"header_fields": {"shared_field": {"type": "Float"}}}
        },
    }
    override = variants.get_annotation_header_fields_override(annotations=annotations)
    assert override == {"shared_field": {"type": "Float"}}


def test_get_annotation_header_fields_override_missing_argument_raises(tmp_path):
    """
    'annotations' is a required argument (no default) -> calling without it
    raises TypeError
    """
    variants = _get_variants_object(tmp_path)
    with pytest.raises(TypeError):
        variants.get_annotation_header_fields_override()


def test_get_annotation_header_fields_override_real_annovar_param(tmp_path):
    """
    Reproduces the real production param structure
    (param.DIAGGEN.GOMV1_GERMLINE.json): most ANNOVAR annotations are flat
    rename dicts with no 'options', only 'ALL.sites.2015_08' defines a
    'header_fields' override.
    """
    annotations = {
        "avsnp150": {"avsnp150": "dbSNP"},
        "refGene": {},
        "refGeneWithVer": {},
        "ensGene": {},
        "iscaCuratedPathogenic": {},
        "iscaCuratedBenign": {},
        "gwasCatalog": {},
        "IARCTP53": {},
        "knownGene": {"Gene_knownGene": "knownGene"},
        "ALL.sites.2015_08": {
            "options": {
                "header_fields": {
                    "1000genomesALL": {
                        "number": ".",
                        "type": "Float",
                        "description": "1000genomesALL by Annovar",
                    }
                }
            }
        },
        "AFR.sites.2015_08": {},
        "EAS.sites.2015_08": {},
        "AMR.sites.2015_08": {},
        "EUR.sites.2015_08": {},
        "SAS.sites.2015_08": {},
        "snp138NonFlagged": {"snp138NonFlagged": "dbSNPNonFlagged"},
        "clinvar_20210123": {"CLNSIG": "CLINVAR"},
    }
    variants = _get_variants_object(tmp_path)
    override = variants.get_annotation_header_fields_override(annotations=annotations)
    assert override == {
        "1000genomesALL": {
            "number": ".",
            "type": "Float",
            "description": "1000genomesALL by Annovar",
        }
    }


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
    Matches the exact pattern used in production param files, e.g.
    "1000genomesALL": {"description": "1000genomesALL by Annovar"}
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
# Integration: get_annotation_header_fields_override + build_info_with_header_override
#
# Reproduces exactly what annotation_annovar() does at the point where the
# merged/multianno ANNOVAR VCF header is re-parsed:
#
#   annotations = param.get("annotation", {}).get("annovar", {}).get("annotations", {})
#   annotation_header_fields_override = self.get_annotation_header_fields_override(
#       annotations=annotations
#   )
#   ...
#   for ann in annovar_vcf_header.infos:
#       if ann not in self.get_header().infos:
#           ann_info = annovar_vcf_header.infos.get(ann)
#           vcf_reader.infos[ann] = self.build_info_with_header_override(
#               field_name=ann,
#               field_number=ann_info.num,
#               field_type=ann_info.type,
#               field_description=ann_info.desc,
#               field_source=ann_info.source,
#               field_version=ann_info.version,
#               header_fields_override=annotation_header_fields_override,
#           )
#
# ANNOVAR always writes Number=. / Type=String in its own multianno VCF
# output, regardless of the actual nature of the annotated value - this is
# exactly why the override is needed.
###############################################################################


class _FakeAnnovarInfo:
    """
    Minimal stand-in for a vcf.parser._Info object as produced by parsing
    ANNOVAR's own multianno VCF header (annovar_vcf_header.infos.get(ann)).
    ANNOVAR always writes Number='.' and Type='String' for custom fields.
    """

    def __init__(self, id_, num=".", type_="String", desc="", source="", version=""):
        self.id = id_
        self.num = num
        self.type = type_
        self.desc = desc
        self.source = source
        self.version = version


def test_annovar_integration_no_header_fields_keeps_annovar_defaults(tmp_path):
    """
    No 'header_fields' configured for any annotation -> ANNOVAR's own
    Number=./Type=String/Description are kept as-is.
    """
    variants = _get_variants_object(tmp_path)
    annotations = {"ALL.sites.2015_08": {"ALL.sites.2015_08": "1000genomesALL"}}
    override = variants.get_annotation_header_fields_override(annotations=annotations)

    ann_info = _FakeAnnovarInfo(
        id_="1000genomesALL",
        num=".",
        type_="String",
        desc="1000genomesALL generated by ANNOVAR",
    )
    info = variants.build_info_with_header_override(
        field_name=ann_info.id,
        field_number=ann_info.num,
        field_type=ann_info.type,
        field_description=ann_info.desc,
        field_source=ann_info.source,
        field_version=ann_info.version,
        header_fields_override=override,
    )

    assert info.type == "String"
    assert info.num == "."
    assert info.desc == "1000genomesALL generated by ANNOVAR"


def test_annovar_integration_real_param_overrides_type_and_description(tmp_path):
    """
    Real production param (param.DIAGGEN.GOMV1_GERMLINE.json): ANNOVAR
    writes '1000genomesALL' as Type=String, but the per-annotation
    'options.header_fields' overrides it to Type=Float with a custom
    description, Number staying '.'.
    """
    annotations = {
        "ALL.sites.2015_08": {
            "annotation_fields": {"ALL.sites.2015_08": "1000genomesALL"},
            "options": {
                "header_fields": {
                    "1000genomesALL": {
                        "number": ".",
                        "type": "Float",
                        "description": "1000genomesALL by Annovar",
                    }
                }
            },
        }
    }
    variants = _get_variants_object(tmp_path)
    override = variants.get_annotation_header_fields_override(annotations=annotations)

    # Simulate ANNOVAR's own multianno VCF header entry (Type=String always)
    ann_info = _FakeAnnovarInfo(
        id_="1000genomesALL",
        num=".",
        type_="String",
        desc="1000genomesALL generated by ANNOVAR",
        source="unknown",
        version="unknown",
    )
    info = variants.build_info_with_header_override(
        field_name=ann_info.id,
        field_number=ann_info.num,
        field_type=ann_info.type,
        field_description=ann_info.desc,
        field_source=ann_info.source,
        field_version=ann_info.version,
        header_fields_override=override,
    )

    assert info.id == "1000genomesALL"
    assert info.type == "Float"  # overridden
    assert info.num == "."  # overridden (same value, but explicitly set)
    assert info.desc == "1000genomesALL by Annovar"  # overridden
    assert info.type_code == variants.code_type_map["Float"]


def test_annovar_integration_other_field_not_overridden_untouched(tmp_path):
    """
    Override only targets '1000genomesALL' (via ALL.sites.2015_08's own
    options): a DIFFERENT ANNOVAR field (e.g. 'refGene', from another
    annotation with no header_fields) must remain untouched (Type=String,
    as produced by ANNOVAR).
    """
    annotations = {
        "ALL.sites.2015_08": {
            "annotation_fields": {"ALL.sites.2015_08": "1000genomesALL"},
            "options": {"header_fields": {"1000genomesALL": {"type": "Float"}}},
        },
        "refGene": {},
    }
    variants = _get_variants_object(tmp_path)
    override = variants.get_annotation_header_fields_override(annotations=annotations)

    ann_info = _FakeAnnovarInfo(
        id_="refGene", num=".", type_="String", desc="Gene refGene annotation"
    )
    info = variants.build_info_with_header_override(
        field_name=ann_info.id,
        field_number=ann_info.num,
        field_type=ann_info.type,
        field_description=ann_info.desc,
        field_source=ann_info.source,
        field_version=ann_info.version,
        header_fields_override=override,
    )

    assert info.type == "String"
    assert info.desc == "Gene refGene annotation"


def test_annovar_integration_via_real_param_file_structure(tmp_path):
    """
    Full pipeline test reproducing exactly what annotation_annovar() does:
    fetch 'annotations' from param.annotation.annovar.annotations (as set
    via variants.set_param(...)), then get_annotation_header_fields_override
    and build_info_with_header_override, in sequence.
    """
    param = {
        "annotation": {
            "annovar": {
                "annotations": {
                    "avsnp150": {"avsnp150": "dbSNP"},
                    "refGene": {},
                    "ALL.sites.2015_08": {
                        "annotation_fields": {
                            "ALL.sites.2015_08": "1000genomesALL"
                        },
                        "options": {
                            "header_fields": {
                                "1000genomesALL": {
                                    "number": ".",
                                    "type": "Float",
                                    "description": "1000genomesALL by Annovar",
                                }
                            }
                        },
                    },
                },
                "options": {"parallelize": "parallel"},
            }
        }
    }
    variants = _get_variants_object(tmp_path, param=param)

    # Reproduce annotation_annovar()'s own fetch of 'annotations'
    annotations = (
        variants.get_param()
        .get("annotation", {})
        .get("annovar", {})
        .get("annotations", {})
    )
    override = variants.get_annotation_header_fields_override(annotations=annotations)

    # Overridden field
    ann_info_overridden = _FakeAnnovarInfo(
        id_="1000genomesALL", num=".", type_="String", desc="ANNOVAR default desc"
    )
    info_overridden = variants.build_info_with_header_override(
        field_name=ann_info_overridden.id,
        field_number=ann_info_overridden.num,
        field_type=ann_info_overridden.type,
        field_description=ann_info_overridden.desc,
        field_source="unknown",
        field_version="unknown",
        header_fields_override=override,
    )
    assert info_overridden.type == "Float"
    assert info_overridden.desc == "1000genomesALL by Annovar"

    # Non-overridden field
    ann_info_untouched = _FakeAnnovarInfo(
        id_="dbSNP", num=".", type_="String", desc="ANNOVAR default desc for dbSNP"
    )
    info_untouched = variants.build_info_with_header_override(
        field_name=ann_info_untouched.id,
        field_number=ann_info_untouched.num,
        field_type=ann_info_untouched.type,
        field_description=ann_info_untouched.desc,
        field_source="unknown",
        field_version="unknown",
        header_fields_override=override,
    )
    assert info_untouched.type == "String"
    assert info_untouched.desc == "ANNOVAR default desc for dbSNP"
# -*- coding: utf-8 -*-
"""Tests for Docker-based variant annotation helpers."""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

import pytest

from howard.functions.commons import (
    get_assembly_mapping_config,
    normalize_assembly_mapping_source,
    resolve_assembly_mapping,
    inside_docker_container
)
from howard.objects.variants_mixin import annotation_docker as annotation_docker_module
from howard.objects.variants_mixin.annotation_docker import (
    _translate_mount_item_host_path,
    _resolve_entry_assembly,
    _resolve_entry_genome,
    _resolve_entry_spec,
    variants_annotation_docker,
)

from test_needed import tests_config


class DummyDocker(variants_annotation_docker):
    def __init__(self, config: dict, param: dict | None = None, tmp_dir: str | None = None):
        self.config = config
        self.param = param or {"assembly": "hg19"}
        self.tmp_dir = tmp_dir or "."
        self.updated = []
        self.export_calls = []

    def get_config(self):
        return self.config

    def get_param(self):
        return self.param

    def get_threads(self):
        return 4

    def get_memory(self, default=None):
        return default or "8G"

    def get_tmp_dir(self):
        return self.tmp_dir

    def _has_variants_to_annotate(self):
        return True

    def _has_variants_for_where_clause(self, where_clause=None):
        return True

    def export_variant_vcf(self, **kwargs):
        self.export_calls.append(kwargs)
        Path(kwargs["vcf_file"]).write_text("dummy\n")

    def update_from_vcf(self, *args, **kwargs):
        self.updated.append((args, kwargs))


@pytest.fixture()
def genome_folder(tmp_path: Path) -> Path:
    folder = tmp_path / "genomes"
    folder.mkdir()
    (folder / "hg19.fa").write_text(">chr1\nACGT\n")
    (folder / "hg38.fa").write_text(">chr1\nTGCA\n")
    return folder


@pytest.fixture()
def base_config(genome_folder: Path) -> dict:
    config = tests_config.copy()
    config = {
        **config,
        "assembly_mapping": {
            "default_source": "howard",
            "aliases": ["howard", "generic", "default", "HG"],
            "sources": {
                "GRCH": {"hg19": "GRCh37", "hg38": "GRCh38"},
                "CUSTOM": {"hg19": "CUSTOM37"},
            },
        },
        "folders": {
            **config["folders"],
            "databases": {
                **config["folders"]["databases"],
                "genomes": str(genome_folder),
            },
        },
        "tools": {
            "vep": {
                "docker": {
                    "image": "ensemblorg/vep:latest",
                    "command": "vep",
                    "annotation": {
                        "primary": {
                            "input": "--input_file",
                            "output": "--output_file",
                            "threads": "--fork",
                            "assembly": {
                                "flag": "--assembly",
                                "source": "GRCH",
                            },
                            "genome": {
                                "flag": "--fasta",
                                "source": "howard",
                            },
                        }
                    },
                }
            },
            "custom_tool": {
                "docker": {
                    "image": "custom/tool:latest",
                    "command": "custom-tool",
                    "annotation": {
                        "primary": {
                            "input": "--vcf",
                            "output": "--out",
                            "assembly": {
                                "flag": "--build",
                                "source": "CUSTOM",
                                "mapping": {"hg19": "OVERRIDE37"},
                            },
                            "genome": {
                                "flag": "--reference",
                                "source": "howard",
                            },
                        }
                    },
                }
            },
        },
    }
    return config


def test_assembly_mapping_helpers_default_identity():
    assembly_mapping_config = get_assembly_mapping_config({})
    assert normalize_assembly_mapping_source(None) == "howard"
    assert normalize_assembly_mapping_source("HG") == "howard"
    assert resolve_assembly_mapping(
        "hg19",
        source="howard",
        mapping=None,
        assembly_mapping_config=assembly_mapping_config,
    ) == "hg19"
    assert resolve_assembly_mapping(
        "hg19",
        source="GRCH",
        mapping=None,
        assembly_mapping_config=assembly_mapping_config,
    ) == "GRCh37"


def test_resolve_entry_spec_supports_custom_mapping_and_genome(base_config: dict):
    entry = {"tool": "custom_tool"}
    spec = _resolve_entry_spec(
        entry_name="custom",
        entry=entry,
        config=base_config,
        default_threads=8,
        default_memory="8G",
    )

    resolved_assembly = _resolve_entry_assembly(
        config=base_config,
        entry_name="custom",
        assembly="hg19",
        assembly_cfg=spec["assembly"],
    )
    assert resolved_assembly == "OVERRIDE37"

    resolved_genome, auto_mount = _resolve_entry_genome(
        config=base_config,
        entry_name="custom",
        assembly="hg19",
        genome_cfg=spec["genome"],
    )
    assert resolved_genome.endswith("hg19.fa")
    assert auto_mount is True


def test_empty_genome_block_disables_genome_resolution(base_config: dict):
    genome_path, auto_mount = _resolve_entry_genome(
        config=base_config,
        entry_name="no_genome",
        assembly="hg19",
        genome_cfg={},
    )

    assert genome_path is None
    assert auto_mount is False


def test_run_entry_injects_assembly_genome_and_auto_mount(monkeypatch, base_config: dict, tmp_path: Path):
    dummy = DummyDocker(config=base_config, tmp_dir=str(tmp_path))
    entry = {"tool": "vep"}
    spec = _resolve_entry_spec(
        entry_name="vep_entry",
        entry=entry,
        config=base_config,
        default_threads=8,
        default_memory="8G",
    )

    captured = {}

    def fake_prepare(self, config, tool, spec, run_dir, assembly, entry_name, command_in_container):
        captured["spec"] = dict(spec)
        captured["command_in_container"] = command_in_container
        captured["run_dir"] = run_dir
        captured["assembly"] = assembly
        return "docker run fake"

    class FakeResult:
        stdout = ""
        stderr = ""

        def check_returncode(self):
            return None

    monkeypatch.setattr(DummyDocker, "_ensure_docker_image", lambda *args, **kwargs: None)
    monkeypatch.setattr(annotation_docker_module.subprocess, "run", lambda *args, **kwargs: FakeResult())
    monkeypatch.setattr(DummyDocker, "_prepare_docker_command", fake_prepare)

    dummy._run_annotation_docker_entry(
        config=base_config,
        entry_name="vep_entry",
        spec=spec,
        assembly="hg19",
    )

    assert "--assembly GRCh37" in captured["command_in_container"]
    assert f"--fasta {base_config['folders']['databases']['genomes']}/hg19.fa" in captured["command_in_container"]
    assert captured["spec"]["resolved_assembly"] == "GRCh37"
    assert captured["spec"]["resolved_genome_path"].endswith("hg19.fa")
    assert captured["spec"]["genome_mounts"] == [
        {
            "host_path": captured["spec"]["resolved_genome_path"],
            "container_path": captured["spec"]["resolved_genome_path"],
            "mode": "ro",
        }
    ]


def test_run_entry_adds_rw_genome_directory_mount(monkeypatch, base_config: dict, tmp_path: Path):
    config = {
        **base_config,
        "tools": {
            **base_config["tools"],
            "vep": {
                "docker": {
                    **base_config["tools"]["vep"]["docker"],
                    "runtime": {
                        "databases": [
                            {
                                "name": "genomes",
                                "container_path": "/opt/vep/.genomes",
                                "mode": "rw",
                            }
                        ],
                    },
                    "annotation": {
                        **base_config["tools"]["vep"]["docker"]["annotation"],
                        "primary": {
                            **base_config["tools"]["vep"]["docker"]["annotation"]["primary"],
                            "genome": {
                                "flag": "--fasta",
                                "source": "howard",
                                "mode": "rw",
                            },
                        },
                    },
                },
            },
        },
    }
    dummy = DummyDocker(config=config, tmp_dir=str(tmp_path))
    entry = {"tool": "vep"}
    spec = _resolve_entry_spec(
        entry_name="vep_entry",
        entry=entry,
        config=config,
        default_threads=8,
        default_memory="8G",
    )

    captured = {}

    def fake_prepare(self, config, tool, spec, run_dir, assembly, entry_name, command_in_container):
        captured["spec"] = dict(spec)
        captured["command_in_container"] = command_in_container
        return "docker run fake"

    class FakeResult:
        stdout = ""
        stderr = ""

        def check_returncode(self):
            return None

    monkeypatch.setattr(DummyDocker, "_ensure_docker_image", lambda *args, **kwargs: None)
    monkeypatch.setattr(annotation_docker_module.subprocess, "run", lambda *args, **kwargs: FakeResult())
    monkeypatch.setattr(DummyDocker, "_prepare_docker_command", fake_prepare)

    dummy._run_annotation_docker_entry(
        config=config,
        entry_name="vep_entry",
        spec=spec,
        assembly="hg19",
    )

    assert "--fasta /opt/vep/.genomes/hg19.fa" in captured["command_in_container"]
    assert captured["spec"]["resolved_genome_path"] == "/opt/vep/.genomes/hg19.fa"
    assert captured["spec"]["genome_mounts"] == [
        {
            "host_path": captured["spec"]["resolved_genome_host_path"],
            "container_path": "/opt/vep/.genomes/hg19.fa",
            "mode": "ro",
        },
        {
            "host_path": str(Path(captured["spec"]["resolved_genome_host_path"]).parent),
            "container_path": "/opt/vep/.genomes",
            "mode": "rw",
        },
    ]


def test_prepare_docker_command_deduplicates_genome_mount(monkeypatch, base_config: dict, tmp_path: Path):
    dummy = DummyDocker(config=base_config, tmp_dir=str(tmp_path))
    spec = {
        "tool": "vep",
        "threads": 4,
        "memory": "8G",
        "mounts": [base_config["folders"]["databases"]["genomes"]],
        "databases": [],
        "paths": [],
        "genome_mounts": [
            {
                "host_path": f"{base_config['folders']['databases']['genomes']}/hg19.fa",
                "container_path": f"{base_config['folders']['databases']['genomes']}/hg19.fa",
                "mode": "ro",
            }
        ],
    }

    monkeypatch.setattr(
        annotation_docker_module,
        "get_bin_command",
        lambda **kwargs: kwargs["add_options"],
    )

    docker_cmd = dummy._prepare_docker_command(
        config=base_config,
        tool="vep",
        spec=spec,
        run_dir=str(tmp_path),
        assembly="hg19",
        entry_name="vep_entry",
        command_in_container="vep --input_file input.vcf --output_file output.vcf.gz",
    )

    genome_mount = f"-v {base_config['folders']['databases']['genomes']}/hg19.fa:{base_config['folders']['databases']['genomes']}/hg19.fa:ro"
    assert docker_cmd.count(genome_mount) == 1
    assert f"-v {tmp_path}:{tmp_path}:rw" in docker_cmd


def test_translate_mount_item_host_path_from_parent_mounts():
    parent_mounts = [
        {
            "host_path": "/home/howard/databases",
            "container_path": "/root/howard/databases",
        }
    ]

    translated = _translate_mount_item_host_path(
        {
            "host_path": "/root/howard/databases/vep/current/hg19",
            "container_path": "/opt/vep/.vep",
            "mode": "ro",
        },
        parent_mounts,
    )

    assert translated == {
        "host_path": "/home/howard/databases/vep/current/hg19",
        "container_path": "/opt/vep/.vep",
        "mode": "ro",
    }


def test_prepare_docker_command_adds_volumes_from_inside_container(monkeypatch, base_config: dict, tmp_path: Path):
    dummy = DummyDocker(config=base_config, tmp_dir=str(tmp_path))
    spec = {
        "tool": "vep",
        "threads": 4,
        "memory": "8G",
        "mounts": [],
        "databases": [],
        "paths": [],
        "genome_mounts": [],
    }

    monkeypatch.setattr(annotation_docker_module, "inside_docker_container", lambda: True)
    monkeypatch.setattr(annotation_docker_module, "get_container_id", lambda: "howard-cli")
    monkeypatch.setattr(
        annotation_docker_module,
        "get_bin_command",
        lambda **kwargs: kwargs["add_options"],
    )

    docker_cmd = dummy._prepare_docker_command(
        config=base_config,
        tool="vep",
        spec=spec,
        run_dir=str(tmp_path),
        assembly="hg19",
        entry_name="vep_entry",
        command_in_container="vep --input_file input.vcf --output_file output.vcf.gz",
    )

    assert "--volumes-from howard-cli" in docker_cmd


def test_prepare_docker_command_translates_specific_db_mount_to_host(monkeypatch, base_config: dict, tmp_path: Path):

    if inside_docker_container():

        config = {
            **base_config,
            "tools": {
                **base_config["tools"],
                "vep": {
                    "docker": {
                        **base_config["tools"]["vep"]["docker"],
                        "options": "-v /home/databases:/root/howard/databases:rw",
                    }
                },
            },
        }

        dummy = DummyDocker(config=config, tmp_dir=str(tmp_path))
        spec = {
            "tool": "vep",
            "threads": 4,
            "memory": "8G",
            "mounts": [],
            "databases": [
                {
                    "path": "/root/howard/databases/vep/current",
                    "container_path": "/opt/vep/.vep",
                }
            ],
            "paths": [],
            "genome_mounts": [],
        }

        monkeypatch.setattr(
            annotation_docker_module,
            "get_bin_command",
            lambda **kwargs: kwargs["add_options"],
        )

        docker_cmd = dummy._prepare_docker_command(
            config=config,
            tool="vep",
            spec=spec,
            run_dir=str(tmp_path),
            assembly="hg19",
            entry_name="vep_online",
            command_in_container="vep --input_file input.vcf --output_file output.vcf.gz",
        )

        assert "-v /home/databases/vep/current/hg19:/opt/vep/.vep:ro" in docker_cmd
        assert "-v /root/howard/databases/vep/current/hg19:/opt/vep/.vep:ro" not in docker_cmd

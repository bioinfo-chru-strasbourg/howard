import copy
import glob
import json
import os
import shlex
import subprocess
from tempfile import TemporaryDirectory
from typing import Any

import logging as log

from howard.functions.commons import (
	DEFAULT_TOOLS_FOLDER,
	check_docker_image_exists,
	command,
	extract_memory_in_go,
	get_bin_command,
	get_random,
)


def _as_list(value: Any) -> list:
	if value is None:
		return []
	if isinstance(value, list):
		return value
	return [value]


def _shell_join(parts: list[str]) -> str:
	return " ".join(shlex.quote(str(part)) for part in parts if part is not None)


def _format_cli_param(parameter: Any) -> list[str]:
	"""
	Convert one parameter description to CLI tokens.

	Supported forms:
	- "--flag"
	- ["--opt", "value"]
	- {"key": "--opt", "value": "x"}
	- {"key": "--flag", "value": true}
	"""

	if parameter is None:
		return []

	if isinstance(parameter, str):
		return [parameter]

	if isinstance(parameter, list):
		return [str(item) for item in parameter]

	if isinstance(parameter, dict):
		key = parameter.get("key", None)
		value = parameter.get("value", None)
		if key is None:
			return []
		if value is None:
			return [str(key)]
		if isinstance(value, bool):
			return [str(key)] if value else []
		return [str(key), str(value)]

	return [str(parameter)]


def _collect_cli_params(parameters: Any) -> list[str]:
	cli_tokens: list[str] = []
	for parameter in _as_list(parameters):
		cli_tokens.extend(_format_cli_param(parameter))
	return cli_tokens


def _cli_block_override_key(block: Any) -> str | None:
	"""
	Return an override key for one CLI block when possible.

	Examples:
	- "--cache" -> "--cache"
	- ["--compress_output", "bgzip"] -> "--compress_output"
	- {"key": "--compress_output", "value": "bgzip"} -> "--compress_output"
	"""

	if isinstance(block, dict):
		key = block.get("key", None)
		return str(key) if key else None

	if isinstance(block, list) and len(block) > 0:
		first = str(block[0])
		if first.startswith("-"):
			return first
		return None

	if isinstance(block, str):
		first = block.strip().split()[0] if block.strip() else ""
		if first.startswith("-"):
			return first
		return None

	return None


def _merge_cli_param_blocks(
	default_parameters: Any,
	default_options: Any,
	entry_parameters: Any,
	entry_options: Any,
) -> list:
	"""
	Merge parameter blocks from defaults and entry values.
	
	Order and precedence:
	- defaults first (parameters + options)
	- entry values next (parameters + options)

	If a later block has the same CLI key (e.g. --compress_output), it overrides
	the previous one. Non-keyed blocks are de-duplicated by exact representation.
	"""
	merged: list[Any] = []
	seen_nonkey: set[str] = set()
	index_by_key: dict[str, int] = {}

	for block in (
		_as_list(default_parameters)
		+ _as_list(default_options)
		+ _as_list(entry_parameters)
		+ _as_list(entry_options)
	):
		override_key = _cli_block_override_key(block)
		if override_key:
			if override_key in index_by_key:
				merged[index_by_key[override_key]] = block
			else:
				index_by_key[override_key] = len(merged)
				merged.append(block)
			continue

		block_repr = json.dumps(block, sort_keys=True, ensure_ascii=True, default=str)
		if block_repr in seen_nonkey:
			continue
		seen_nonkey.add(block_repr)
		merged.append(block)

	return merged


def _expand_path(path: str) -> str:
	return os.path.abspath(os.path.expanduser(path))


def _path_contains_assembly(path: str, assembly: str) -> bool:
	if not assembly:
		return True
	path_norm = os.path.normpath(path).lower()
	assembly_norm = str(assembly).lower()
	return any(assembly_norm in part for part in path_norm.split(os.sep))


def _extract_db_paths_from_item(db_item: Any, databases_config: dict) -> list[str]:
	if isinstance(db_item, str):
		if db_item in databases_config:
			return _as_list(databases_config.get(db_item))
		return [db_item]

	if isinstance(db_item, dict):
		if db_item.get("path", None):
			return [db_item.get("path")]
		db_key = db_item.get("name", db_item.get("key", None))
		if db_key:
			return _as_list(databases_config.get(db_key))

	return []


def _extract_db_mount_options(db_item: Any, db_mount_defaults: dict) -> tuple[str | None, str]:
	container_path = None
	mode = "ro"

	if isinstance(db_item, dict):
		container_path = db_item.get("container_path", None)
		mode = str(db_item.get("mode", mode))

	if db_mount_defaults:
		if container_path is None:
			container_path = db_mount_defaults.get("container_path", None)
		if mode == "ro":
			mode = str(db_mount_defaults.get("mode", mode))

	return container_path, mode


def _extract_db_key(db_item: Any, databases_config: dict) -> str | None:
	if isinstance(db_item, str) and db_item in databases_config:
		return db_item
	if isinstance(db_item, dict):
		return db_item.get("name", db_item.get("key", None))
	return None


def _resolve_assembly_path(
	entry_name: str,
	resolved_path: str,
	assembly: str | None,
) -> str:
	if not assembly or _path_contains_assembly(resolved_path, assembly):
		return resolved_path

	assembly_subpath = _expand_path(os.path.join(resolved_path, assembly))
	if os.path.exists(assembly_subpath):
		return assembly_subpath

	raise ValueError(
		f"annotation_docker entry '{entry_name}': database path '{resolved_path}' does not include assembly '{assembly}' and '{assembly_subpath}' does not exist"
	)


def _resolve_database_paths(
	entry_name: str,
	config: dict,
	databases_value: Any,
	assembly: str | None,
) -> list[dict]:
	resolved_paths: list[dict] = []
	databases_config = config.get("folders", {}).get("databases", {})
	databases_mounts_config = config.get("folders", {}).get("databases_mounts", {})

	for db_item in _as_list(databases_value):
		db_paths = _extract_db_paths_from_item(db_item=db_item, databases_config=databases_config)
		db_key = _extract_db_key(db_item=db_item, databases_config=databases_config)
		db_mount_defaults = {}
		if db_key:
			db_mount_defaults = databases_mounts_config.get(db_key, {})

		container_path, mode = _extract_db_mount_options(
			db_item=db_item,
			db_mount_defaults=db_mount_defaults,
		)
		for db_path in db_paths:
			if not db_path:
				continue

			resolved_path = _expand_path(str(db_path))
			resolved_path = _resolve_assembly_path(
				entry_name=entry_name,
				resolved_path=resolved_path,
				assembly=assembly,
			)

			entry_container_path = container_path if container_path else resolved_path
			resolved_paths.append(
				{
					"host_path": resolved_path,
					"container_path": entry_container_path,
					"mode": mode,
				}
			)

	return resolved_paths


def _format_mount(path_item: Any, default_mode: str = "ro") -> str | None:
	if path_item is None:
		return None

	if isinstance(path_item, str):
		path_str = path_item.strip()
		if not path_str:
			return None
		if path_str.startswith("-v "):
			return path_str
		host_path = _expand_path(path_str)
		return f"-v {host_path}:{host_path}:{default_mode}"

	if isinstance(path_item, dict):
		host_value = path_item.get("host_path", path_item.get("path", None))
		if not host_value:
			return None
		host_path = _expand_path(str(host_value))
		container_path = str(path_item.get("container_path", host_path))
		mode = str(path_item.get("mode", default_mode))
		return f"-v {host_path}:{container_path}:{mode}"

	path_str = str(path_item)
	if not path_str:
		return None
	host_path = _expand_path(path_str)
	return f"-v {host_path}:{host_path}:{default_mode}"


def _warn_entry_forbidden_overrides(entry_name: str, entry: dict, tool: str) -> None:
	location_map = {
		"primary": f"config.tools.{tool}.docker.annotation.primary",
		"databases": f"config.tools.{tool}.docker.runtime.databases",
		"mounts": f"config.tools.{tool}.docker.runtime.mounts",
		"mount": f"config.tools.{tool}.docker.runtime.mounts",
		"paths": f"config.tools.{tool}.docker.runtime.paths",
		"path": f"config.tools.{tool}.docker.runtime.paths",
		"remove_info": f"config.tools.{tool}.docker.annotation.vcf_update.remove_info",
		"add_samples": f"config.tools.{tool}.docker.annotation.vcf_update.add_samples",
		"update_header": f"config.tools.{tool}.docker.annotation.vcf_update.update_header",
		"command": f"config.tools.{tool}.docker.command",
	}

	for forbidden_key in [
		"primary",
		"databases",
		"mounts",
		"mount",
		"paths",
		"path",
		"remove_info",
		"add_samples",
		"update_header",
		"command",
	]:
		if forbidden_key not in entry:
			continue
		location = location_map.get(forbidden_key, f"config.tools.{tool}.docker")

		log.warning(
			f"annotation_docker entry '{entry_name}': key '{forbidden_key}' should be configured in {location} (entry value ignored)"
		)


def _normalize_runtime_mount_settings(runtime_cfg: dict, annotation_cfg: dict) -> tuple[list, list, list]:
	"""
	Resolve runtime mount settings with backward compatibility and aliases.

	Preferred location:
	- docker.runtime.databases/mounts/paths

	Backward-compatible fallbacks:
	- docker.annotation.databases/mounts/paths
	- docker.runtime.mount (alias mounts)
	- docker.runtime.path (alias paths)
	"""

	databases = _as_list(
		runtime_cfg.get("databases", annotation_cfg.get("databases", []))
	)

	runtime_mounts = runtime_cfg.get("mounts", None)
	if runtime_mounts is None:
		runtime_mounts = runtime_cfg.get("mount", None)
	if runtime_mounts is None:
		runtime_mounts = annotation_cfg.get("mounts", [])
	mounts = _as_list(runtime_mounts)

	runtime_paths = runtime_cfg.get("paths", None)
	if runtime_paths is None:
		runtime_paths = runtime_cfg.get("path", None)
	if runtime_paths is None:
		runtime_paths = annotation_cfg.get("paths", annotation_cfg.get("path", []))
	paths = _as_list(runtime_paths)

	return databases, mounts, paths


def _get_effective_threads_limit(default_threads: int | None) -> int:
	"""Compute the effective thread limit for one entry."""

	system_threads = os.cpu_count() or 1
	if not default_threads or int(default_threads) <= 0:
		return int(system_threads)

	return int(min(int(default_threads), int(system_threads)))


def _get_effective_memory_limit_gb(default_memory: Any) -> int:
	"""Compute the effective memory limit (GB) for one entry."""

	default_memory_gb = extract_memory_in_go(default_memory)
	if default_memory_gb < 1:
		default_memory_gb = 1

	try:
		import psutil  # type: ignore

		available_memory_gb = int(psutil.virtual_memory().available / 1024 / 1024 / 1024)
		if available_memory_gb < 1:
			available_memory_gb = 1
		return int(min(default_memory_gb, available_memory_gb))
	except Exception:
		return int(default_memory_gb)


def _normalize_entry_resources(
	entry_name: str,
	entry_threads: Any,
	entry_memory: Any,
	default_memory: Any,
	threads_limit: int,
	memory_limit_gb: int,
) -> tuple[int, str]:
	"""Normalize and cap entry resources to effective limits."""

	resolved_threads = int(entry_threads) if entry_threads is not None else int(threads_limit)
	if resolved_threads <= 0:
		resolved_threads = int(threads_limit)
	if resolved_threads > threads_limit:
		log.warning(
			f"annotation_docker entry '{entry_name}': requested threads={resolved_threads} exceeds available threads={threads_limit}, using {threads_limit}"
		)
		resolved_threads = int(threads_limit)

	if entry_memory is None:
		entry_memory = default_memory
	requested_memory_gb = extract_memory_in_go(entry_memory)
	if requested_memory_gb < 1:
		requested_memory_gb = 1
	if requested_memory_gb > memory_limit_gb:
		log.warning(
			f"annotation_docker entry '{entry_name}': requested memory={requested_memory_gb}G exceeds available memory={memory_limit_gb}G, using {memory_limit_gb}G"
		)
		requested_memory_gb = int(memory_limit_gb)

	return int(resolved_threads), f"{int(requested_memory_gb)}G"


def _resolve_entry_spec(
	entry_name: str,
	entry: dict,
	config: dict,
	default_threads: int,
	default_memory: str,
) -> dict:
	"""
	Normalize one annotation_docker entry.

	Expected minimal fields:
	- tool: tool name in config.tools

	Configuration-centric behavior:
	- Structural runtime settings are read from config.tools.<tool>.docker.annotation
	  (primary, databases, mounts, paths, remove_info, add_samples, update_header).
	- Param entry keeps user-specific run settings (parameters, options, where_clause),
	  with allowed overrides for resources and output_pattern.
	"""

	tool = entry.get("tool", None)
	if not tool:
		raise ValueError(f"annotation_docker entry '{entry_name}': missing 'tool'")

	tool_config = config.get("tools", {}).get(tool, {})
	if not tool_config:
		raise ValueError(
			f"annotation_docker entry '{entry_name}': tool '{tool}' not configured in config.tools"
		)

	docker_cfg = tool_config.get("docker", {})
	annotation_cfg = docker_cfg.get("annotation", {})
	runtime_cfg = docker_cfg.get("runtime", {})
	vcf_update_cfg = annotation_cfg.get("vcf_update", {})
	defaults_cfg = annotation_cfg.get("defaults", {})

	if not isinstance(runtime_cfg, dict):
		runtime_cfg = {}
	if not isinstance(vcf_update_cfg, dict):
		vcf_update_cfg = {}
	if not isinstance(defaults_cfg, dict):
		defaults_cfg = {}

	_warn_entry_forbidden_overrides(entry_name=entry_name, entry=entry, tool=tool)

	command_name = docker_cfg.get("command", None)

	if not command_name:
		log.warning(
			f"annotation_docker entry '{entry_name}': empty tool docker command (will rely on docker entrypoint if configured)"
		)

	primary = annotation_cfg.get("primary", {})
	if not isinstance(primary, dict):
		raise ValueError(
			f"annotation_docker entry '{entry_name}': config.tools.{tool}.docker.annotation.primary must be a dict"
		)

	input_flag = primary.get("input", None)
	output_flag = primary.get("output", None)
	threads_flag = primary.get("threads", None)
	memory_flag = primary.get("memory", None)

	if not input_flag:
		raise ValueError(
			f"annotation_docker entry '{entry_name}': missing input flag in config.tools.{tool}.docker.annotation.primary.input"
		)
	if not output_flag:
		raise ValueError(
			f"annotation_docker entry '{entry_name}': missing output flag in config.tools.{tool}.docker.annotation.primary.output"
		)

	cfg_resources = defaults_cfg.get("resources", annotation_cfg.get("resources", {}))
	entry_resources = entry.get("resources", {})
	if not isinstance(cfg_resources, dict):
		cfg_resources = {}
	if not isinstance(entry_resources, dict):
		entry_resources = {}

	entry_threads = entry_resources.get(
		"threads",
		cfg_resources.get("threads", default_threads),
	)
	entry_memory = entry_resources.get(
		"memory",
		cfg_resources.get("memory", default_memory),
	)

	if entry_threads is None:
		entry_threads = default_threads
	if entry_memory is None:
		entry_memory = default_memory

	threads_limit = _get_effective_threads_limit(default_threads=default_threads)
	memory_limit_gb = _get_effective_memory_limit_gb(default_memory=default_memory)
	entry_threads, entry_memory = _normalize_entry_resources(
		entry_name=entry_name,
		entry_threads=entry_threads,
		entry_memory=entry_memory,
		default_memory=default_memory,
		threads_limit=threads_limit,
		memory_limit_gb=memory_limit_gb,
	)

	default_parameters = defaults_cfg.get(
		"parameters",
		annotation_cfg.get("parameters", []),
	)
	default_options = defaults_cfg.get(
		"options",
		annotation_cfg.get("options", []),
	)
	entry_options = entry.get("options", [])
	entry_parameters = entry.get("parameters", [])

	if entry_parameters and entry_options:
		log.warning(
			f"annotation_docker entry '{entry_name}': both 'parameters' and 'options' are provided; values are merged and exact duplicates are removed"
		)

	merged_param_blocks = _merge_cli_param_blocks(
		default_parameters=default_parameters,
		default_options=default_options,
		entry_parameters=entry_parameters,
		entry_options=entry_options,
	)
	databases, mounts, paths = _normalize_runtime_mount_settings(
		runtime_cfg=runtime_cfg,
		annotation_cfg=annotation_cfg,
	)

	remove_info = vcf_update_cfg.get("remove_info", annotation_cfg.get("remove_info", True))
	add_samples = vcf_update_cfg.get("add_samples", annotation_cfg.get("add_samples", False))
	update_header = vcf_update_cfg.get("update_header", annotation_cfg.get("update_header", True))
	update_existing_fields = vcf_update_cfg.get(
		"update_existing_fields",
		annotation_cfg.get("update_existing_fields", False),
	)
	update_existing_fields = bool(
		entry.get("update_existing_fields", update_existing_fields)
	)

	return {
		"tool": tool,
		"command": command_name,
		"input_flag": input_flag,
		"output_flag": output_flag,
		"threads_flag": threads_flag,
		"memory_flag": memory_flag,
		"threads": entry_threads,
		"memory": entry_memory,
		"parameters": _collect_cli_params(merged_param_blocks),
		"where_clause": entry.get("where_clause", None),
		"remove_info": remove_info,
		"add_samples": add_samples,
		"update_header": update_header,
		"update_existing_fields": update_existing_fields,
		"output_pattern": entry.get(
			"output_pattern",
			defaults_cfg.get("output_pattern", annotation_cfg.get("output_pattern", None)),
		),
		"mounts": mounts,
		"databases": databases,
		"paths": paths,
	}


class variants_annotation_docker:
	"""
	Generic docker-based annotation mixin.

	Parameter schema:
	section -> annotation_docker -> entries -> <entry_name>
	"""

	def _find_output_vcf(
		self,
		output_vcf_expected: str,
		run_dir: str,
		output_pattern: str | None,
	) -> str | None:
		if os.path.exists(output_vcf_expected):
			return output_vcf_expected

		candidates = []
		if output_pattern:
			pattern = (
				output_pattern
				if os.path.isabs(output_pattern)
				else os.path.join(run_dir, output_pattern)
			)
			candidates.extend(sorted(glob.glob(pattern)))

		if not candidates:
			candidates.extend(sorted(glob.glob(os.path.join(run_dir, "*.vcf"))))
			candidates.extend(sorted(glob.glob(os.path.join(run_dir, "*.vcf.gz"))))

		for candidate in candidates:
			if os.path.exists(candidate):
				return candidate

		return None

	def _get_annotation_docker_entries(self, section: str) -> dict:
		param = self.get_param()
		entries = (
			param.get(section, {})
			.get("docker", {})
			.get("entries", {})
		)
		if not entries:
			log.info("No annotation_docker entries configured")
		return entries

	def _has_variants_to_annotate(self) -> bool:
		table_variants = self.get_table_variants()
		sql_query_chromosomes = (
			f"""SELECT count(*) as count FROM {table_variants} as table_variants"""
		)
		return bool(self.get_query_to_df(sql_query_chromosomes)["count"][0])

	def _has_variants_for_where_clause(self, where_clause: str | None = None) -> bool:
		table_variants = self.get_table_variants()
		where_sql = where_clause if where_clause else ""
		query = f"SELECT 1 as has_variant FROM {table_variants}{where_sql} LIMIT 1"
		query_df = self.get_query_to_df(query)
		return len(query_df) > 0

	def _ensure_docker_image(self, entry_name: str, docker_image: str) -> None:
		# check_docker_image_exists expects image in "name:tag" format.
		# If no tag is provided, Docker defaults to "latest".
		docker_image_for_check = docker_image
		image_last_part = docker_image.rsplit("/", 1)[-1]
		if ":" not in image_last_part and "@" not in image_last_part:
			docker_image_for_check = f"{docker_image}:latest"

		try:
			if check_docker_image_exists(docker_image_for_check):
				return
		except ValueError:
			# Fallback for uncommon docker references (e.g. digests or custom formats).
			log.debug(
				f"annotation_docker entry '{entry_name}': skip local image check for '{docker_image}'"
			)

		log.warning(
			f"Annotation: docker image {docker_image} not found locally, trying to pull"
		)
		try:
			command(f"docker pull {docker_image}")
		except subprocess.CalledProcessError:
			msg_err = f"Unable to pull docker image {docker_image}"
			log.error(msg_err)
			raise ValueError(
				f"annotation_docker entry '{entry_name}': {msg_err}"
			)

	def _prepare_docker_command(
		self,
		config: dict,
		tool: str,
		spec: dict,
		run_dir: str,
		assembly: str | None,
		entry_name: str,
		command_in_container: str | None,
	) -> str:
		config_runtime = copy.deepcopy(config)
		if command_in_container is not None:
			config_runtime["tools"][tool]["docker"]["command"] = command_in_container

		proxy = [
			f"-e {var}={os.getenv(var)}"
			for var in ["https_proxy", "http_proxy", "ftp_proxy"]
			if os.getenv(var) is not None
		]
		entry_mounts = [f"-v {run_dir}:{run_dir}:rw"]

		for mount in spec["mounts"]:
			formatted_mount = _format_mount(mount, default_mode="rw")
			if formatted_mount:
				entry_mounts.append(formatted_mount)

		for db_path in _resolve_database_paths(
			entry_name=entry_name,
			config=config,
			databases_value=spec.get("databases", []),
			assembly=assembly,
		):
			formatted_mount = _format_mount(db_path, default_mode="ro")
			if formatted_mount:
				entry_mounts.append(formatted_mount)

		for path_item in spec.get("paths", []):
			formatted_mount = _format_mount(path_item, default_mode="ro")
			if formatted_mount:
				entry_mounts.append(formatted_mount)

		run_name = f"HOWARD-ANNOT-{entry_name}-{get_random()}"
		docker_add_options = (
			f"--name {run_name} {' '.join(entry_mounts)} {' '.join(proxy)}"
		)

		return get_bin_command(
			tool=tool,
			bin_type="docker",
			config=config_runtime,
			param={"threads": spec["threads"], "memory": spec["memory"]},
			default_folder=f"{DEFAULT_TOOLS_FOLDER}/docker",
			add_options=docker_add_options,
		)

	def _run_annotation_docker_entry(
		self,
		config: dict,
		entry_name: str,
		spec: dict,
		assembly: str | None,
	) -> None:
		tool = spec["tool"]

		tool_config = config.get("tools", {}).get(tool, {})
		if not tool_config:
			raise ValueError(
				f"annotation_docker entry '{entry_name}': tool '{tool}' not configured in config.tools"
			)

		docker_cfg = tool_config.get("docker", {})
		docker_image = docker_cfg.get("image", None)
		docker_entrypoint = docker_cfg.get("entrypoint", None)
		docker_config_command = docker_cfg.get("command", None)
		# command is treated as an executable path (entrypoint-like), not as a full command line.
		effective_command = docker_config_command
		if not docker_image:
			raise ValueError(
				f"annotation_docker entry '{entry_name}': missing docker image in config.tools.{tool}.docker.image"
			)
		if (
			not effective_command
			and not docker_entrypoint
		):
			raise ValueError(
				f"annotation_docker entry '{entry_name}': missing executable 'command' and no docker entrypoint configured in config.tools.{tool}.docker"
			)

		self._ensure_docker_image(entry_name=entry_name, docker_image=docker_image)

		with TemporaryDirectory(
			prefix=f"annotation_docker-{entry_name}-", dir=self.get_tmp_dir()
		) as run_dir:
			if not self._has_variants_for_where_clause(
				where_clause=spec.get("where_clause", None)
			):
				log.debug(
					f"annotation_docker entry '{entry_name}': no variants to annotate after where_clause"
				)
				return None

			input_vcf = os.path.join(run_dir, "input.vcf")
			output_vcf_expected = os.path.join(run_dir, "output.vcf.gz")

			self.export_variant_vcf(
				vcf_file=input_vcf,
				remove_info=spec["remove_info"],
				add_samples=spec["add_samples"],
				index=False,
				where_clause=spec["where_clause"],
			)

			command_in_container = None
			cmd_tokens = []

			# Build command from configured executable path + generated args.
			# If executable is missing but docker entrypoint exists, pass only args.
			if effective_command:
				cmd_tokens.append(str(effective_command))
				cmd_tokens.extend([spec["input_flag"], input_vcf])
				cmd_tokens.extend(spec["parameters"])
				cmd_tokens.extend([spec["output_flag"], output_vcf_expected])
			else:
				cmd_tokens.extend([spec["input_flag"], input_vcf])
				cmd_tokens.extend(spec["parameters"])
				cmd_tokens.extend([spec["output_flag"], output_vcf_expected])

			if spec["threads_flag"] and spec["threads"] is not None:
				cmd_tokens.extend([spec["threads_flag"], str(spec["threads"])])

			if spec["memory_flag"] and spec["memory"] is not None:
				cmd_tokens.extend([spec["memory_flag"], str(spec["memory"])])

			command_in_container = _shell_join(cmd_tokens)
			log.debug(
				f"annotation_docker entry '{entry_name}' command: {command_in_container}"
			)

			docker_cmd = self._prepare_docker_command(
				config=config,
				tool=tool,
				spec=spec,
				run_dir=run_dir,
				assembly=assembly,
				entry_name=entry_name,
				command_in_container=command_in_container,
			)

			log.debug(f"annotation_docker entry '{entry_name}' docker cmd: {docker_cmd}")
			res = subprocess.run(
				docker_cmd, shell=True, capture_output=True, text=True
			)
			log.debug(res.stdout)
			if res.stderr:
				log.error(res.stderr)
			res.check_returncode()

			output_vcf = self._find_output_vcf(
				output_vcf_expected=output_vcf_expected,
				run_dir=run_dir,
				output_pattern=spec["output_pattern"],
			)
			if not output_vcf:
				raise ValueError(
					f"annotation_docker entry '{entry_name}': output VCF not found in '{run_dir}'"
				)

			self.update_from_vcf(
				output_vcf,
				update_header=spec["update_header"],
				update_existing_fields=spec["update_existing_fields"],
			)

	def annotation_docker(self, section: str = "annotation", threads: int = None) -> None:
		"""
		Generic Docker annotation runner.

		It exports variants to VCF, executes one or multiple dockerized tools,
		then imports annotations back with update_from_vcf.
		"""

		log.debug("Start annotation with generic docker tools")

		config = self.get_config()

		if not threads:
			threads = self.get_threads()
		log.debug(f"Threads: {threads}")

		default_memory = self.get_memory("8G")
		log.debug(f"Default memory: {default_memory}")

		param = self.get_param()
		assembly = param.get("assembly", config.get("assembly", None))
		log.debug(f"Assembly for annotation_docker: {assembly}")

		entries = self._get_annotation_docker_entries(section=section)
		if not entries:
			return None

		if not self._has_variants_to_annotate():
			log.info("VCF empty")
			return None

		for entry_name, entry in entries.items():
			log.info(f"Annotation 'annotation_docker/{entry_name}' ...")

			spec = _resolve_entry_spec(
				entry_name=entry_name,
				entry=entry,
				config=config,
				default_threads=threads,
				default_memory=default_memory,
			)

			self._run_annotation_docker_entry(
				config=config,
				entry_name=entry_name,
				spec=spec,
				assembly=assembly,
			)

		return None

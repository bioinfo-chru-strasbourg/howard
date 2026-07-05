import copy
import glob
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


def _resolve_entry_spec(entry_name: str, entry: dict, default_threads: int, default_memory: str) -> dict:
	"""
	Normalize one annotation_docker entry.

	Expected minimal fields:
	- tool: tool name in config.tools
	- command: main command in container
	- primary.input: input flag (e.g. --input)
	- primary.output: output flag (e.g. --output)
	"""

	tool = entry.get("tool", None)
	command_name = entry.get("command", None)

	if not tool:
		raise ValueError(f"annotation_docker entry '{entry_name}': missing 'tool'")
	if not command_name:
		log.warning(
			f"annotation_docker entry '{entry_name}': empty 'command' (will rely on docker entrypoint if configured)"
		)

	primary = entry.get("primary", {})
	if not isinstance(primary, dict):
		raise ValueError(
			f"annotation_docker entry '{entry_name}': 'primary' must be a dict"
		)

	input_flag = primary.get("input", entry.get("input", None))
	output_flag = primary.get("output", entry.get("output", None))
	threads_flag = primary.get("threads", entry.get("threads", None))
	memory_flag = primary.get("memory", entry.get("memory", None))

	if not input_flag:
		raise ValueError(
			f"annotation_docker entry '{entry_name}': missing input flag (primary.input)"
		)
	if not output_flag:
		raise ValueError(
			f"annotation_docker entry '{entry_name}': missing output flag (primary.output)"
		)

	entry_threads = entry.get("resources", {}).get("threads", default_threads)
	entry_memory = entry.get("resources", {}).get("memory", default_memory)

	if entry_threads is None:
		entry_threads = default_threads
	if entry_memory is None:
		entry_memory = default_memory

	options = entry.get("options", [])
	parameters = entry.get("parameters", [])

	return {
		"tool": tool,
		"command": command_name,
		"input_flag": input_flag,
		"output_flag": output_flag,
		"threads_flag": threads_flag,
		"memory_flag": memory_flag,
		"threads": entry_threads,
		"memory": entry_memory,
		"parameters": _collect_cli_params(parameters) + _collect_cli_params(options),
		"where_clause": entry.get("where_clause", None),
		"remove_info": entry.get("remove_info", True),
		"add_samples": entry.get("add_samples", True),
		"update_header": entry.get("update_header", True),
		"output_pattern": entry.get("output_pattern", None),
		"mounts": _as_list(entry.get("mounts", [])),
		"databases": _as_list(entry.get("databases", [])),
		"paths": _as_list(entry.get("paths", [])),
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
		entry_command = spec.get("command", None)
		# command is treated as an executable path (entrypoint-like), not as a full command line.
		effective_command = entry_command if entry_command else docker_config_command
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

			# Build command from executable path + generated args.
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

			self.update_from_vcf(output_vcf, update_header=spec["update_header"])

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

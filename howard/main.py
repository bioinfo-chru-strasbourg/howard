#!/usr/bin/env python

import os
import argparse
import logging as log
import psutil  # type: ignore
import json
import yaml  # type: ignore

from howard.tools.tools import (
    commands_arguments,
    arguments,
    shared_arguments,
    tool_gui_enable,
)
from howard.functions.plugins import plugins_infos, plugins_list, plugins_to_load
from howard.functions.commons import (
    folder_plugins,
    full_path,
    get_log_level,
    load_args,
    subfolder_plugins,
    help_header,
    set_log_level,
    help_generation,
    DEFAULT_CHUNK_SIZE,
)


# Usage
# python -m pip install -e .
# howard --help
# howard query --help
# howard process --input=my.vcf.gz --output=my.output.vcf --annotations=my.annotations.vcf.gz
# howard gui
# python -m howard.main --input=my.vcf.gz --output=my.output.vcf --annotations=my.annotations.vcf.gz


msg_gui_disable = "HOWARD GUI disabled"

# Load command submodule
for command in list(commands_arguments.keys()):
    try:
        exec(
            "from {module} import {submodule}".format(
                module="howard.tools.tools", submodule=command
            )
        )
    except Exception as e:
        print(e)


# Load plugns
if os.path.exists(folder_plugins):

    # Load plugins infos
    plugins = plugins_infos(
        plugins_dir=folder_plugins, subfolder_plugins=subfolder_plugins
    )
    list_plugins_dict = plugins_list(
        plugins, plugins_dir=folder_plugins, subfolder_plugins=subfolder_plugins
    )
    plugins_to_load_dict = plugins_to_load(plugins_list_dict=list_plugins_dict)

    # Import plugins
    for plugin_name in plugins_to_load_dict:

        # Plugin infos
        plugin_infos = plugins_to_load_dict[plugin_name]
        plugin_main_file = plugin_infos.get("__main_file__", "__main__")
        plugin_main_function = plugin_infos.get("__main_function__", "main")
        plugin_arguments = plugin_infos.get("__arguments__", "arguments")
        plugins_commands_arguments = plugin_infos.get(
            "__commands_arguments__", "commands_arguments"
        )

        try:
            # Load plugin as module
            exec(
                "from {module}.{plugin_name}.{main_file} import {main_function} as {plugin_name}, {plugin_arguments} as plugin_arguments, {plugins_commands_arguments} as plugins_commands_arguments".format(
                    module=subfolder_plugins,
                    main_function=plugin_main_function,
                    plugin_name=plugin_name,
                    main_file=plugin_main_file,
                    plugin_arguments=plugin_arguments,
                    plugins_commands_arguments=plugins_commands_arguments,
                )
            )

            # Identify as plugin
            plugins_commands_arguments[plugin_name]["help"] = "(plugin) " + str(
                plugins_commands_arguments.get(plugin_name, {}).get("help")
            )

            # Add command arguments on argparse and help
            for plugin_command_arguments in plugins_commands_arguments:
                if plugin_command_arguments not in commands_arguments:
                    commands_arguments[plugin_command_arguments] = (
                        plugins_commands_arguments[plugin_command_arguments]
                    )

            # Add arguments need for plugin
            for plugin_argument in plugin_arguments:
                if plugin_argument not in arguments:
                    arguments[plugin_argument] = plugin_arguments[plugin_argument]

        # Warning if plugon not loaded
        except Exception as inst:
            log.warning(f"{inst}")
            msg_warning = f"Plugin '{plugin_name}' NOT loaded"
            log.warning(msg_warning)

main_folder = os.path.dirname(__file__)


# Main function
def main() -> None:
    """
    It loads a VCF file in multiple format (VCF, parquet, DB), and process, query, export data
    """

    # Config infos
    setup_cfg = f"{main_folder}/../setup.cfg"
    arguments_dict = {
        "arguments": arguments,
        "commands_arguments": commands_arguments,
        "shared_arguments": shared_arguments,
    }

    # Generate parser
    parser = argparse.ArgumentParser()
    parser = help_generation(
        arguments_dict=arguments_dict,
        parser=parser,
        output_type="parser",
        setup=setup_cfg,
    )

    # Parse args
    args = parser.parse_args()

    # Quiet
    if "quiet" in args and args.quiet:
        args.verbosity = "warning"
    # Verbose
    if "verbose" in args and args.verbose:
        args.verbosity = "info"
    # Debug
    if "debug" in args and args.debug:
        args.verbosity = "debug"
    # verbosity
    if "verbosity" not in args:
        args.verbosity = "info"

    # log
    if "log" not in args:
        args.log = None
    elif args.log and not isinstance(args.log, str):
        args.log = args.log.name

    # Config infos
    args.arguments_dict = arguments_dict
    args.setup_cfg = setup_cfg

    # Logging
    set_log_level(args.verbosity, args.log)

    # Header
    if get_log_level() <= 20:
        print(help_header(setup=setup_cfg))

    # Threads
    nb_threads = os.cpu_count()
    if "threads" in args and args.threads:
        threads = args.threads
    else:
        threads = nb_threads
    if threads == -1:
        threads = nb_threads

    # Memory
    mem = psutil.virtual_memory()
    mem_total = mem.total / 1024 / 1024 / 1024
    mem_default = int(mem_total * 0.8)
    if mem_default < 1:
        mem_default = 1
    if "memory" in args and args.memory:
        memory = args.memory
    else:
        memory = f"{mem_default}G"

    # chunk_size
    chunk_size = DEFAULT_CHUNK_SIZE
    if "chunk_size" in args and args.chunk_size:
        chunk_size = int(args.chunk_size)

    # Temporary folder
    tmp = None
    if "tmp" in args and args.tmp:
        tmp = args.tmp

    # duckdb settings
    duckdb_settings = None
    if "duckdb_settings" in args and args.duckdb_settings:
        if isinstance(args.duckdb_settings, str) and os.path.exists(
            full_path(args.duckdb_settings)
        ):
            with open(full_path(args.duckdb_settings)) as config_file:
                duckdb_settings = yaml.safe_load(config_file)
        else:
            duckdb_settings = json.loads(args.duckdb_settings)

    # Assembly
    if "assembly" in args and args.assembly:
        assembly = args.assembly
    else:
        assembly = "hg19"

    # Load configuration in JSON format
    if "config" in args:
        if isinstance(args.config, str) and os.path.exists(full_path(args.config)):
            with open(full_path(args.config)) as config_file:
                config = yaml.safe_load(config_file)
        else:
            config = json.loads(args.config)
    else:
        config = {}

    # add to config

    # Verbosity
    config["verbosity"] = args.verbosity

    # Threads
    if "threads" not in config or not config.get("threads", None):
        config["threads"] = threads

    # Memory
    if "memory" not in config or not config.get("memory", None):
        config["memory"] = memory

    # Chunk size
    if "chunk_size" not in config or not config.get("chunk_size", None):
        config["chunk_size"] = chunk_size

    # Tmp
    if "tmp" not in config or not config.get("tmp", None):
        config["tmp"] = tmp

    # duckDB settings
    if "duckdb_settings" not in config or not config.get("duckdb_settings", None):
        config["duckdb_settings"] = duckdb_settings

    # Assembly
    if "assembly" not in config or not config.get("assembly", None):
        config["assembly"] = assembly

    # Change config
    args.config = config
    log.debug(f"config: {config}")

    # Variants object
    vcfdata_obj = None

    # Command eval
    command_function = None
    if not args.command:
        parser.print_help()
    else:
        if args.command == "gui" and not tool_gui_enable:
            log.error(msg_gui_disable)
            log.error("""ModuleNotFoundError: No module named 'gooey'""")
            log.error("""Install module 'gooey': pip install gooey""")
            log.error(
                """Or install requirements: pip install -r requirements-gui.txt"""
            )
            raise ValueError(msg_gui_disable)
        command_function = commands_arguments[args.command]["function"]

        log.debug(f"Command/Tool: {command_function}")
        vcfdata_obj = eval(f"{command_function}(args)")

    # Interactive option
    interactive = False

    # Load args into param
    if vcfdata_obj is not None:
        param = vcfdata_obj.get_param()
    else:
        param = {}

    param = load_args(
        param=param,
        args=args,
        arguments_dict=arguments_dict,
        command="query",
        strict=False,
    )

    # Query
    query = param.get("query", {}).get("query", None)

    # Specific "query" tool interactivity
    # if command_function == "query" and ("query" not in args or not args.query):
    if command_function == "query" and not query:
        interactive = True
        log.debug("Interactivity terminal activated")

    # Launch interactive terminal if --interactive is specified
    if interactive or ("interactive" in args and args.interactive):

        # Check variants object
        if vcfdata_obj is None:
            msg_err = f"Command/Tool '{command_function}' does not support interactive terminal"
            log.warning(msg_err)
            return None

        # Import
        from howard.tools.interactive import launch_interactive_terminal

        # Launch interactive terminal
        log.info("Start interative terminal")
        launch_interactive_terminal(args=args, variants=vcfdata_obj)
        log.info("End interative terminal")

    # Close Variants object connexion
    if vcfdata_obj is not None:
        log.debug("Close connexion")
        vcfdata_obj.close_connexion()
        log.debug("Connexion closed")


if __name__ == "__main__":
    main()

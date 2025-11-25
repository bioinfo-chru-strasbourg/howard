import argparse
import glob
import logging as log
import os
import tempfile
from pathlib import Path
import shutil
from tabulate import tabulate  # type: ignore

from howard.functions.commons import (
    DEFAULT_CHUNK_SIZE,
    load_args,
    load_config_args,
    get_tmp,
    get_random,
    remove_if_exists,
)
from howard.objects.variants import Variants


def process(args: argparse) -> None:
    """
    The "process" function processes input arguments, loads parameters in JSON format, creates a VCF
    object, performs annotations, calculations, prioritizations, and queries, exports output, and
    closes the connection.

    This enhanced version supports processing large files by:
    1. Optionally chunking the input file into smaller temporary files
    2. Processing each chunk independently
    3. Merging the results at the end

    :param args: args is a variable that contains the arguments passed to the function "process". It is
    assumed to be an object with several attributes, including "config", "param", "input", "output",
    "annotations", "calculations", "prioritizations", and "query". These attributes are used to
    :type args: argparse
    """

    log.info("Start")

    # Load config args
    arguments_dict, _, config, param = load_config_args(args)

    # Load args into param
    param = load_args(
        param=param,
        args=args,
        arguments_dict=arguments_dict,
        command="process",
        strict=False,
    )

    # Get chunking parameters from the dedicated "chunking" section in param
    chunking_config_config = config.get("chunking", {})
    chunking_config_param = param.get("chunking", {})
    if isinstance(chunking_config_config, dict) and isinstance(
        chunking_config_param, dict
    ):
        chunking_config = {**chunking_config_config, **chunking_config_param}
    elif isinstance(chunking_config_config, dict):
        chunking_config = chunking_config_config
    elif isinstance(chunking_config_param, dict):
        chunking_config = chunking_config_param
    else:
        chunking_config = {}

    # Check if chunking is enabled via the 'enable' flag in the chunking section
    chunking_enable = chunking_config.get("chunking_enable", False)

    # Check chunking mode from args (overrides config/param)
    chunking_mode = chunking_config.get("chunking_mode", "parquet")

    # If chunking is not enabled, run the standard processing flow
    if not chunking_enable:
        return _process_standard(args, config, param)
    else:
        # Log chunking settings
        log.info("Chunking enabled")
        for key, value in chunking_config.items():
            if key != "chunking_enable":
                log.info(f"  {key}: {value}")

        # Call chunked processing in parquet mode
        if chunking_mode.lower() == "parquet" or not chunking_mode:
            return _process_chunked(args, config, param, chunking_config)
        elif chunking_mode.lower() == "duckdb":
            log.debug(f"config1: {config}")
            return _process_duckdb(args, config, param, chunking_config)
        else:
            log.error(f"Unknown chunking mode: {chunking_mode}")
            raise ValueError(f"Unknown chunking mode: {chunking_mode}")


def _process_duckdb(args, config, param, chunking_config):
    """
    Process a large input file by chunking it into smaller parts using DuckDB's capabilities.

    :param args: Command line arguments
    :param config: Configuration dictionary
    :param param: Parameters dictionary
    :param chunking_config: Chunking configuration dictionary
    :return: Processed Variants object
    """

    # Currently not implemented
    log.warning("Chunking mode 'duckdb' is an experimental feature.")

    config_chunking_mode_duckdb = config.copy()
    config_chunking_mode_duckdb["connexion_type"] = "tmpfile"

    # Get chunking parameters from chunking configuration
    # chunk_size (default to 1,000,000 if not found, see Variants class for details)
    chunk_size = (
        chunking_config.get("chunking_size")
        or (
            args.chunk_size
            if "chunk_size" in args and args.chunk_size != DEFAULT_CHUNK_SIZE
            else None
        )
        or param.get("chunk_size")
        or config.get("chunk_size")
        or DEFAULT_CHUNK_SIZE
    )
    config_chunking_mode_duckdb["chunk_size"] = chunk_size

    return _process_standard(args, config_chunking_mode_duckdb, param)


def _process_standard(args, config, param):
    """
    Standard processing function without chunking.

    :param args: Command line arguments
    :param config: Configuration dictionary
    :param param: Parameters dictionary
    :return: Processed Variants object
    """

    # Create variants object
    vcfdata_obj = Variants(
        input=args.input, output=args.output, config=config, param=param
    )

    # Get Config and Params
    config = vcfdata_obj.get_config()
    param = vcfdata_obj.get_param()

    # Re-Load Config and Params
    vcfdata_obj.set_param(param)
    vcfdata_obj.set_config(config)

    # Load data
    vcfdata_obj.load_data()

    # Annotation HGVS
    vcfdata_obj.annotation_hgvs()

    # Annotations
    vcfdata_obj.annotation()

    # Calculations
    vcfdata_obj.calculation()

    # Prioritization
    vcfdata_obj.prioritization()

    # Explode infos
    if param.get("explode", {}).get("explode_infos", False):
        vcfdata_obj.explode_infos()

    # Query
    if param.get("query", {}).get("query", None):
        _handle_query(vcfdata_obj, param)

    # Export
    vcfdata_obj.export_output(query=param.get("query", {}).get("query", None))

    # Close connexion
    # vcfdata_obj.close_connexion()

    # Log
    log.info("End")

    # Return Variants object
    return vcfdata_obj


def _process_chunked(args, config, param, chunking_config):
    """
    Process a large input file by chunking it into smaller parts, using HOWARD's built-in
    partitioning functionality.

    This implementation follows the approach from the bash script:
    1. Convert input to partitioned parquet format using howard convert
    2. Process each partition independently
    3. Merge results with the convert tool

    :param args: Command line arguments
    :param config: Configuration dictionary
    :param param: Parameters dictionary
    :param chunking_config: Chunking configuration dictionary
    :return: Final merged Variants object
    """

    # Store original input and output
    original_input = args.input
    original_output = args.output

    # Chunking config

    # Get chunking parameters from chunking configuration
    # chunk_size (default to 1,000,000 if not found, see Variants class for details)
    chunk_size = (
        chunking_config.get("chunking_size")
        or (
            args.chunk_size
            if "chunk_size" in args and args.chunk_size != DEFAULT_CHUNK_SIZE
            else None
        )
        or param.get("chunk_size")
        or config.get("chunk_size")
        or DEFAULT_CHUNK_SIZE
    )
    log.debug(f"Chunking - Using chunk size: {chunk_size}")

    # Get chunking partitioning strategy (default to 'None' if not specified)
    # Can be 'None', '#CHROM', '#CHROM,FILTER', etc.
    # Prevents if chunking_partitions is empty string or None
    chunking_partitions = chunking_config.get("chunking_partitions", "None")
    if chunking_partitions is None or str(chunking_partitions).strip() == "":
        chunking_partitions = "None"
    log.debug(f"Chunking - Using partitioning strategy: {chunking_partitions}")

    # Get base temporary directory
    base_temp_dir = get_tmp(config=config, param=param)

    # Create a dedicated subdirectory for chunking within the temp directory
    chunk_temp_dir = os.path.join(base_temp_dir, f"howard_chunks_{get_random()}")
    os.makedirs(chunk_temp_dir, exist_ok=True)

    # Create temp directory for chunks
    temp_dir = tempfile.mkdtemp(dir=chunk_temp_dir)
    log.debug(f"Chunking - Created temporary directory for chunks: {temp_dir}")

    # Define paths for partitioned data
    partition_dir = os.path.join(temp_dir, "partitioned.parquet")

    try:
        # Create an object to convert input to partitioned parquet format
        log.debug(
            f"Chunking - Converting input file to partitioned parquet format with chunk size {chunk_size}"
        )

        # Simplified approach that mimics the bash script
        # Create a copy of param and config to modify the access mode to RO for chunking
        chunk_param = param.copy() if isinstance(param, dict) else {}
        chunk_config = config.copy() if isinstance(config, dict) else {}

        # Transcripts export handling
        if param.get("transcripts", {}).get("export", None):
            transcripts_export = True
        else:
            transcripts_export = False

        # Transcripts export
        if transcripts_export:
            log.debug("Chunking - Transcripts export chunk processing enabled")
            original_transcripts_export = (
                chunk_param.copy().get("transcripts", {}).get("export", {})
            )
            original_transcripts_export_output = original_transcripts_export.get(
                "output", None
            )

        # Set access mode to RO to reduce memory usage during chunking
        if isinstance(chunk_param, dict):
            chunk_param["access"] = "RO"

        # Set access mode in config as well, which is required for the RO mode to work
        if isinstance(chunk_config, dict):
            chunk_config["access"] = "RO"

        # Force connexion type in memory to allow RO access mode
        chunk_config["connexion_type"] = "memory"

        # Create Variants object for conversion into Parquet with chunking
        convert_obj = Variants(
            input=original_input,
            output=partition_dir,
            config=chunk_config,
            param=chunk_param,
        )

        # Load the data in read-only mode
        convert_obj.load_data()

        # Check number of variants to decide if chunking is needed
        log.debug("Chunking - Checking number of variants in input file")
        check_nb_of_variants = convert_obj.get_query_to_df(
            f"SELECT count(1) AS count FROM variants LIMIT ({chunk_size+1})"
        )
        if check_nb_of_variants["count"][0] <= chunk_size:
            log.debug(
                f"Chunking - Number of variants is less than or equal to chunk size ({chunk_size}), processing without chunking"
            )
            return _process_standard(args, config, param)
        else:
            log.debug(
                f"Chunking - Number of variants exceeds chunk size ({chunk_size}), proceeding with chunked processing"
            )

        # Export to partitioned parquet format using parquet_partitions='None'
        # This ensures the file is split into multiple parquet files in the output directory
        convert_obj.export_output(
            output_file=partition_dir,
            parquet_partitions=chunking_partitions,
            chunk_size=chunk_size,
            export_header=True,
        )

        # Close connexion for conversion object
        convert_obj.close_connexion()

        # Find all parquet chunks
        # Use glob to find all ".parquet" file in a folder, recursively
        glob_pattern = os.path.join(partition_dir, "**", "*.parquet")
        chunk_files = sorted(
            glob.glob(glob_pattern, recursive=True),
            key=lambda p: (os.path.getmtime(p)),
        )
        log.debug(f"chunk_files: {chunk_files}")

        # If no chunk files found, raise error
        if len(chunk_files) <= 1:
            # Suppose no variants were found (empty input file)
            log.debug(
                f"Chunking - Only {len(chunk_files)} parquet chunks created, processing without chunking"
            )
            # Process with standard flow to create empty output
            return _process_standard(args, config, param)
        else:
            log.debug(f"Chunking - Created {len(chunk_files)} parquet chunks")

        # Create directory to hold processed chunks
        processed_dir = os.path.join(temp_dir, "processed_chunks.parquet")
        os.makedirs(processed_dir, exist_ok=True)
        processed_chunks = []
        processed_chunks_header = []
        processed_transcripts_dir = os.path.join(
            temp_dir, "processed_transcripts_chunks.parquet"
        )
        os.makedirs(processed_transcripts_dir, exist_ok=True)
        processed_transcripts_chunks = []

        for i, chunk_file in enumerate(chunk_files):
            log.info(
                f"Chunking - Processing chunk {i+1}/{len(chunk_files)} - {(i+1)/len(chunk_files)*100:.2f}%"
            )
            log.debug(f"Chunking - Processing chunk file {chunk_file}")

            # Copy the header of the converted input partitionned parquet
            # for the chunk parquet file
            if os.path.exists(f"{partition_dir}.hdr"):
                shutil.copy(f"{partition_dir}.hdr", f"{chunk_file}.hdr")

            # Process the chunk
            chunk_output = os.path.join(processed_dir, f"chunk_{i+1}.parquet")

            # Create a new args object to avoid modifying the original
            chunk_args = argparse.Namespace()
            for attr, value in vars(args).items():
                setattr(chunk_args, attr, value)

            # Set input/output for this chunk
            chunk_args.input = chunk_file
            chunk_args.output = chunk_output

            # Transcripts export
            if transcripts_export:
                log.debug("Chunking - Process Transcripts chunk")
                chunk_transcripts_output = os.path.join(
                    processed_transcripts_dir, f"chunk_{i+1}.parquet"
                )
                chunk_param["transcripts"]["export"][
                    "output"
                ] = chunk_transcripts_output
            else:
                chunk_transcripts_output = None

            # Create a copy of parameters for processing individual chunks
            # NOTE: We don't set access=RO here, as we want normal processing for individual chunks
            chunk_proc_param = param.copy() if isinstance(param, dict) else {}
            chunk_proc_config = config.copy() if isinstance(config, dict) else {}

            # Remove query in param to avoid querying each chunk
            if "query" in chunk_proc_param:
                del chunk_proc_param["query"]

            try:
                # Process the chunk with standard processing (not in read-only mode)
                _process_standard(chunk_args, chunk_proc_config, chunk_proc_param)
                processed_chunks.append(chunk_output)
                processed_chunks_header.append(f"{chunk_output}.hdr")
                if chunk_transcripts_output and os.path.exists(
                    chunk_transcripts_output
                ):
                    processed_transcripts_chunks.append(chunk_transcripts_output)

                # Remove original chunk to save space
                remove_if_exists([chunk_file, f"{chunk_file}.hdr"])

            except Exception as e:
                log.error(f"Error processing chunk {i+1}: {str(e)}")
                raise ValueError(f"Error processing chunk {i+1}: {str(e)}")

        # Check if any chunks were processed
        if not processed_chunks:
            log.error("Chunking - Process failed")
            raise RuntimeError("Chunking - Process failed")
        else:
            log.info(f"Chunking - Merging {len(processed_chunks)} processed chunks")

        # Create a merged Variants object for final output

        # Create read-only parameters and config for final merging
        merge_param = param.copy() if isinstance(param, dict) else {}
        merge_config = config.copy() if isinstance(config, dict) else {}

        if isinstance(merge_param, dict):
            merge_param["access"] = "RO"  # Use read-only access for merging

        if isinstance(merge_config, dict):
            merge_config["access"] = "RO"  # Use read-only access in config as well

        # Force connexion type in memory to allow RO access mode
        merge_config["connexion_type"] = "memory"

        # Copy output header files for processed chunks
        if os.path.exists(processed_chunks_header[0]):
            shutil.copy(processed_chunks_header[0], f"{processed_dir}.hdr")

        # Create final Variants object for merging
        final_obj = Variants(
            input=processed_dir,
            output=original_output,
            config=merge_config,
            param=merge_param,
            load=True,
        )

        # Sort if VCF format
        # Needed because chunk and merging parquet files does not guarantee sorted order
        if final_obj.get_output_format() in [
            "vcf",
            "vcf.gz",
            "bcf",
        ] and chunk_param.get("chunking_sort", False):
            log.debug("Output can be sorted as VCF format")
            sort = True
        else:
            log.debug("Output sorting skipped for non-VCF format")
            sort = False

        # Export to final output format
        final_obj.export_output(
            sort=sort, query=param.get("query", {}).get("query", None)
        )

        # Check if final export created the output file
        if os.path.exists(original_output):
            log.debug(f"Chunking - Final output written to {original_output}")
        else:
            log.error(f"Chunking - Final output file {original_output} not created")

        # Transcript export
        if transcripts_export:
            if processed_transcripts_chunks:

                # Create final Variants object for merging transcripts
                final_transcripts_obj = Variants(
                    input=processed_transcripts_dir,
                    output=original_transcripts_export_output,
                    config=merge_config,
                    load=True,
                )
                log.debug(
                    f"Chunking - Exporting final transcripts to {original_transcripts_export_output}"
                )

                # Export transcripts
                final_transcripts_obj.export_output()

                # Close connexion for final object
                final_transcripts_obj.close_connexion()

            elif original_transcripts_export_output:

                log.debug("Chunking - No transcripts to export")
                log.debug(
                    f"Chunking - Creating empty transcripts export file at {original_transcripts_export_output}"
                )
                remove_if_exists([original_transcripts_export_output])
                Path(original_transcripts_export_output).touch()

        # Handle query if needed
        if param.get("query", {}).get("query", None):
            _handle_query(final_obj, param)

        log.info("Chunking - Chunked processing completed successfully")

        # Log
        log.info("End")

        # Return None because chunking do not support returning Variants object (for interactive mode)
        return None

    except Exception as e:
        import traceback

        log.error(f"Chunking - Error during chunked processing: {str(e)}")
        log.error(f"Chunking - Error traceback: {traceback.format_exc()}")
        raise RuntimeError(f"Chunking - Error during chunked processing: {str(e)}")

    finally:
        # Clean up and restore original arguments
        try:
            log.debug(f"Chunking - Cleaning up temporary directory: {temp_dir}")
            shutil.rmtree(temp_dir)
        except Exception as e:
            log.warning(
                f"Chunking - Failed to clean up temporary directory {temp_dir}: {str(e)}"
            )

        try:
            log.debug(f"Chunking - Cleaning up chunking directory: {chunk_temp_dir}")
            shutil.rmtree(chunk_temp_dir)
        except Exception as e:
            log.warning(
                f"Chunking - Failed to clean up chunking directory {chunk_temp_dir}: {str(e)}"
            )

        args.input = original_input
        args.output = original_output


def _handle_query(vcfdata_obj, param):
    """
    Handle query operations on a Variants object.

    :param vcfdata_obj: Variants object
    :param param: Parameters dictionary
    """
    log.info("Querying...")

    # Parameters
    query = param.get("query", {}).get("query", None)
    query_limit = param.get("query", {}).get("query_limit", None)
    query_print_mode = param.get("query", {}).get("query_print_mode", None)

    # Print query
    if query_print_mode in ["markdown"]:
        print(vcfdata_obj.get_query_to_df(query, limit=query_limit).to_markdown())
    elif query_print_mode in ["tabulate"]:
        print(
            tabulate(
                vcfdata_obj.get_query_to_df(query, limit=query_limit),
                headers="keys",
                tablefmt="psql",
            )
        )
    else:
        print(vcfdata_obj.get_query_to_df(query, limit=query_limit))

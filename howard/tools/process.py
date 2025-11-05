import argparse
import logging as log
import os
import tempfile
from pathlib import Path
import shutil
import time
from tabulate import tabulate  # type: ignore

from howard.functions.commons import (
    load_args,
    load_config_args,
    get_tmp,
    get_threads,
    get_random,
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
    chunking_config = param.get("chunking", {})

    # Check if chunking is enabled via the 'enable' flag in the chunking section
    use_chunking = chunking_config.get("chunking_enable", False)

    # If chunking is not enabled, run the standard processing flow
    if not use_chunking:
        log.info(f"Chunking disabled")
        return _process_standard(args, config, param)
    else:
        # Get chunking parameters from chunking configuration
        # chunk_size (default to 100,000 if not found)
        if "chunk_size" in args:
            chunk_size = args.chunk_size
        else:
            chunk_size = 100000
        chunk_size = chunking_config.get(
            "chunk_size",
            param.get("chunk_size", config.get("chunk_size", chunk_size)),
        )
        chunk_size = chunking_config.get(
            "chunk_size",
            param.get("chunk_size", config.get("chunk_size", chunk_size)),
        )

        # Get base temporary directory
        base_temp_dir = get_tmp(config=config, param=param)

        # Create a dedicated subdirectory for chunking within the temp directory
        import tempfile
        import os

        chunk_temp_dir = os.path.join(base_temp_dir, f"howard_chunks_{get_random()}")
        os.makedirs(chunk_temp_dir, exist_ok=True)

        # Get threads using the existing get_threads utility function
        # threads = get_threads(config=config, param=param)

        # Get chunking partitioning strategy (default to 'None' if not specified)
        # Can be 'None', '#CHROM', '#CHROM,FILTER', etc.
        chunking_partitions = chunking_config.get("partitions", "None")

        # Log chunking settings
        log.info(f"Chunking enabled")
        log.debug(
            f"Chunking enabled: chunk_size={chunk_size}, tmp_dir={chunk_temp_dir}, partitioning={chunking_partitions}"
        )

        # param_chunk = param.copy()

        # # copy param dict to avoid modifying original
        # param_chunk = param.copy() if isinstance(param, dict) else {}

        # # remove query in param
        # if "query" in param_chunk:
        #     del param_chunk["query"]

        # print(param)
        # print(param_chunk)

        # Call chunked processing
        return _process_chunked(
            args,
            config,
            param,
            chunk_size,
            chunk_temp_dir,
            chunking_partitions,
        )


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
    if param.get("explode_infos", {}).get("explode_infos", False):
        vcfdata_obj.explode_infos()

    # Query
    if param.get("query", {}).get("query", None):
        _handle_query(vcfdata_obj, param)

    # Export
    vcfdata_obj.export_output(query=param.get("query", {}).get("query", None))

    # Log
    log.info("End")

    # Return Variants object
    return vcfdata_obj


def _process_chunked(
    args, config, param, chunk_size, chunk_temp_dir, chunking_partitions="None"
):
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
    :param chunk_size: Number of variants per chunk
    :param chunk_temp_dir: Directory to store temporary chunks
    :param chunking_partitions: Strategy for partitioning chunks ('None', '#CHROM', etc.)
    :return: Final merged Variants object
    """
    # Store original input and output
    original_input = args.input
    original_output = args.output

    # Create temp directory for chunks
    temp_dir = tempfile.mkdtemp(dir=chunk_temp_dir)
    log.debug(f"Chunking - Created temporary directory for chunks: {temp_dir}")

    # Log the chunking partition strategy passed from the main function
    log.debug(f"Chunking - Using partitioning strategy: {chunking_partitions}")

    # Define paths for partitioned data
    partition_dir = os.path.join(temp_dir, "partitioned.parquet")

    try:
        # Step 1: Create an object to convert input to partitioned parquet format
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

        if transcripts_export:
            log.debug(f"Chunking - Transcripts export chunk processing enabled")
            original_transcripts_export = (
                chunk_param.copy().get("transcripts", {}).get("export", {})
            )
            original_transcripts_export_output = original_transcripts_export.get(
                "output", None
            )
            # log.debug(
            #     f"original_transcripts_export_output: {original_transcripts_export_output}"
            # )
            # transcripts_partition_dir = os.path.join(
            #     temp_dir, "transcripts.partitioned.parquet"
            # )

        # log.debug(f"original_transcripts_export : {original_transcripts_export}")

        # Set access mode to RO to reduce memory usage during chunking
        if isinstance(chunk_param, dict):
            # Force access mode to RO for chunking to reduce memory usage
            chunk_param["access"] = "RO"

        # Set access mode in config as well, which is required for the RO mode to work
        if isinstance(chunk_config, dict):
            # Force access mode to RO in config for chunking to reduce memory usage
            chunk_config["access"] = "RO"

        convert_obj = Variants(
            input=original_input,
            output=partition_dir,
            config=chunk_config,
            param=chunk_param,
        )

        # Load the data in read-only mode
        convert_obj.load_data()

        # Export to partitioned parquet format using parquet_partitions='None'
        # This ensures the file is split into multiple parquet files in the output directory
        convert_obj.export_output(
            output_file=partition_dir,
            parquet_partitions=chunking_partitions,
            chunk_size=chunk_size,
            export_header=True,
        )

        # Step 2: Find all parquet chunks
        chunk_files = []
        for root, _, files in os.walk(partition_dir):
            for file in files:
                if file.endswith(".parquet"):
                    chunk_files.append(os.path.join(root, file))

        if not chunk_files:
            log.error(
                "Chunking - No parquet chunks created in {}".format(partition_dir)
            )
            raise RuntimeError("Chunking - Failed to create parquet chunks")

        log.debug("Chunking - Created {} parquet chunks".format(len(chunk_files)))

        # Step 3: Process each chunk
        processed_dir = os.path.join(temp_dir, "processed_chunks.parquet")
        os.makedirs(processed_dir, exist_ok=True)
        processed_chunks = []
        processed_transcripts_dir = os.path.join(
            temp_dir, "processed_transcripts_chunks.parquet"
        )
        os.makedirs(processed_transcripts_dir, exist_ok=True)
        processed_transcripts_chunks = []

        for i, chunk_file in enumerate(chunk_files):
            log.info(f"Chunking - Processing chunk {i+1}/{len(chunk_files)}")

            # Copy the header for the chunk
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

            # # Transcripts export handling
            chunk_transcripts_output = None
            if transcripts_export:
                log.debug(f"Chunking - Transcripts export chunk processing - ")
                chunk_transcripts_output = os.path.join(
                    processed_transcripts_dir, f"chunk_{i+1}.parquet"
                )
                chunk_param["transcripts"]["export"][
                    "output"
                ] = chunk_transcripts_output
                log.debug(
                    f"Chunking - Transcripts export path: {chunk_transcripts_output}"
                )
                log.debug(f"chunk_param: {chunk_param}")
            else:
                chunk_transcripts_output = None

            # # Transcript export
            # if param.get("transcripts", {}).get("export", None):
            #     log.debug(
            #         f"Chunking - Enabling transcripts export for chunk {i+1} processing"
            #     )
            #     if "transcripts" not in param:
            #         param["transcripts"] = {}
            #     param["transcripts"]["export"] = original_transcripts_export

            # Create a copy of parameters for processing individual chunks
            # NOTE: We don't set access=RO here, as we want normal processing for individual chunks
            chunk_proc_param = param.copy() if isinstance(param, dict) else {}
            chunk_proc_config = config.copy() if isinstance(config, dict) else {}

            if "query" in chunk_proc_param:
                del chunk_proc_param["query"]

            try:
                # Process the chunk with standard processing (not in read-only mode)
                _process_standard(chunk_args, chunk_proc_config, chunk_proc_param)
                processed_chunks.append(chunk_output)
                if chunk_transcripts_output:
                    processed_transcripts_chunks.append(chunk_transcripts_output)

                # Copy the header from the processed file
                if os.path.exists(f"{chunk_output}.hdr"):
                    shutil.copy(f"{chunk_output}.hdr", f"{processed_dir}.hdr")

                # Optional: Remove original chunk to save space
                # if os.path.exists(chunk_file):
                #     os.remove(chunk_file)
            except Exception as e:
                log.error(f"Error processing chunk {i+1}: {str(e)}")
                continue

        if not processed_chunks:
            log.error("Chunking - No chunks were successfully processed")
            raise RuntimeError("Chunking - All chunks failed to process")

        # Step 4: Create a merged Variants object for final output
        log.info(f"Chunking - Merging {len(processed_chunks)} processed chunks")

        # Create read-only parameters and config for final merging
        merge_param = param.copy() if isinstance(param, dict) else {}
        merge_config = config.copy() if isinstance(config, dict) else {}

        if isinstance(merge_param, dict):
            merge_param["access"] = "RO"  # Use read-only access for merging

        if isinstance(merge_config, dict):
            merge_config["access"] = "RO"  # Use read-only access in config as well

        if len(processed_chunks) == 1 and True:
            # If only one chunk was processed, just copy it
            # Use read-only mode for the final merge operation to reduce memory usage
            # single_chunk_param = param.copy() if isinstance(param, dict) else {}
            # single_chunk_config = config.copy() if isinstance(config, dict) else {}

            # if isinstance(single_chunk_param, dict):
            #     single_chunk_param["access"] = (
            #         "RO"  # Use read-only access for final merge
            #     )

            # if isinstance(single_chunk_config, dict):
            #     single_chunk_config["access"] = (
            #         "RO"  # Use read-only access in config for final merge
            #     )

            final_obj = Variants(
                input=processed_chunks[0],
                output=original_output,
                config=merge_config,
                param=merge_param,
                load=True,
            )

            # Transcript export
            if transcripts_export:
                # log.debug(
                #     f"original_transcripts_export : {original_transcripts_export}"
                # )
                # log.debug(
                #     f"""original_transcripts_export : {original_transcripts_export.get("output", None)}"""
                # )
                # log.debug(
                #     f"original_transcripts_export_output: {original_transcripts_export_output}"
                # )

                # shutil.copy(chunk_file, dest_file)
                final_transcripts_obj = Variants(
                    # input=processed_transcripts_dir,
                    input=processed_transcripts_chunks[0],
                    output=original_transcripts_export_output,
                    config=merge_config,
                    load=True,
                )
                # df = final_transcripts_obj.get_query_to_df(
                #     """ SELECT * FROM variants """
                # )
                # log.debug(f"Final transcripts dataframe  {df}")

        else:
            # Create a merge directory
            merge_dir = os.path.join(temp_dir, "merged")
            os.makedirs(merge_dir, exist_ok=True)

            # Copy all processed files to merge directory
            for i, chunk_file in enumerate(processed_chunks):
                dest_file = os.path.join(merge_dir, f"part_{i}.parquet")
                shutil.copy(chunk_file, dest_file)

                # Copy header from the last processed chunk
                if i == len(processed_chunks) - 1 and os.path.exists(
                    f"{chunk_file}.hdr"
                ):
                    shutil.copy(f"{chunk_file}.hdr", f"{merge_dir}.hdr")

            # Create a final object to export
            # Explicitly mark the directory as a parquet format by naming it with .parquet extension
            parquet_merge_dir = merge_dir + ".parquet"
            os.rename(merge_dir, parquet_merge_dir)

            # Copy header if it exists
            if os.path.exists(f"{merge_dir}.hdr"):
                shutil.copy(f"{merge_dir}.hdr", f"{parquet_merge_dir}.hdr")

            # # Create read-only parameters and config for final merging
            # merge_param = param.copy() if isinstance(param, dict) else {}
            # merge_config = config.copy() if isinstance(config, dict) else {}

            # if isinstance(merge_param, dict):
            #     merge_param["access"] = "RO"  # Use read-only access for merging

            # if isinstance(merge_config, dict):
            #     merge_config["access"] = "RO"  # Use read-only access in config as well

            final_obj = Variants(
                input=parquet_merge_dir,
                output=original_output,
                config=merge_config,
                param=merge_param,
                load=True,
            )

            # Transcript export
            if transcripts_export:
                final_transcripts_obj = Variants(
                    input=processed_transcripts_dir,
                    output=original_transcripts_export_output,
                    config=merge_config,
                    load=True,
                )
                # df = final_transcripts_obj.get_query_to_df(
                #     """ SELECT * FROM variants """
                # )
                # log.debug(f"Final transcripts dataframe  {df}")

        # Sort if VCF format
        if final_obj.get_output_format() in ["vcf", "vcf.gz", "bcf"]:
            log.debug("Output can be sorted as VCF format")
            sort = True
        else:
            log.debug("Output sorting skipped for non-VCF format")
            sort = False

        # Export to final output format
        # log.debug(f"original_output: {original_output}")
        # final_obj.export_output(output_file=original_output, sort=True, index=True)
        # df = final_obj.get_query_to_df(""" SELECT * FROM variants """)
        # log.debug(f"Final dataframe has {len(df)} rows")
        # log.debug(f"Final dataframe columns: {df.to_string()}")
        final_obj.export_output(sort=sort)

        # Transcript export
        if transcripts_export:
            log.debug(
                f"Chunking - Exporting final transcripts to {original_transcripts_export_output}"
            )
            final_transcripts_obj.export_output()

        if os.path.exists(original_output):
            log.debug(f"Chunking - Final output written to {original_output}")
        else:
            log.error(f"Chunking - Final output file {original_output} not created")

        # # Transcripts export
        # if param.get("transcripts", {}).get("export", None):
        #     final_transcripts_obj.export_output()

        # Handle query if needed
        if param.get("query", {}).get("query", None):
            _handle_query(final_obj, param)

        log.info("Chunking - Chunked processing completed successfully")

        # Log
        log.info("End")

        return final_obj

    except Exception as e:
        import traceback

        log.error(f"Chunking - Error during chunked processing: {str(e)}")
        log.error(f"Chunking - Error traceback: {traceback.format_exc()}")
        raise RuntimeError(f"Chunking - Error during chunked processing: {str(e)}")

        # # If chunking fails, fall back to standard processing
        # log.info("Chunking - Falling back to standard processing...")
        # args.input = original_input
        # args.output = original_output
        # return _process_standard(args, config, param)

    finally:
        # Clean up and restore original arguments
        try:
            log.debug(f"Chunking - Cleaning up temporary directory: {temp_dir}")
            # shutil.rmtree(temp_dir)
        except Exception as e:
            log.warning(
                f"Chunking - Failed to clean up temporary directory {temp_dir}: {str(e)}"
            )

        try:
            log.debug(f"Chunking - Cleaning up chunking directory: {chunk_temp_dir}")
            # shutil.rmtree(chunk_temp_dir)
        except Exception as e:
            log.warning(
                f"Chunking - Failed to clean up chunking directory {chunk_temp_dir}: {str(e)}"
            )

        args.input = original_input
        args.output = original_output


# Function removed as it's now integrated directly into _process_chunked


# def _merge_processed_chunks(chunk_files, output_file, config, param):
#     """
#     Merge multiple processed chunk files into a single output file.

#     Note: This function is kept for compatibility but is no longer used in the
#     simplified implementation, as merging is now handled directly in _process_chunked
#     using HOWARD's built-in functionality for handling partitioned parquet files.

#     :param chunk_files: List of processed chunk files
#     :param output_file: Final output file path
#     :param config: Configuration dictionary
#     :param param: Parameters dictionary
#     """
#     log.warning(
#         "_merge_processed_chunks is deprecated, merging now handled in _process_chunked"
#     )

#     if not chunk_files:
#         log.error("No chunks to merge")
#         return

#     log.info(f"Merging {len(chunk_files)} chunks into {output_file}")

#     # Create a temporary directory to hold all processed chunks
#     temp_dir = os.path.dirname(chunk_files[0])
#     partition_dir = os.path.join(temp_dir, "merged_partition")
#     os.makedirs(partition_dir, exist_ok=True)

#     # Copy all chunk files to the partition directory with proper naming
#     for i, chunk_file in enumerate(chunk_files):
#         dest_file = os.path.join(partition_dir, f"part_{i}.parquet")
#         shutil.copy(chunk_file, dest_file)

#         # Copy header file if it exists
#         if os.path.exists(f"{chunk_file}.hdr"):
#             shutil.copy(f"{chunk_file}.hdr", f"{partition_dir}.hdr")

#     # Create a Variants object for the partition directory
#     merged_obj = Variants(
#         input=partition_dir, output=output_file, config=config, param=param, load=True
#     )

#     # Export to final format
#     merged_obj.export_output(
#         output_file=output_file,
#         sort=True,  # Sort by genomic coordinates
#         index=True,  # Create index if applicable
#     )


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

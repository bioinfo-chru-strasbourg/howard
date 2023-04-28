




def parquet_info_explode(input_file:str, output_file:str, threads:int = 1, tmp:str = "/tmp") -> None:
    """
    > It takes a parquet file, splits it by chromosome, explodes the INFO column, and then merges the
    exploded files back together
    
    :param input_file: The input file to be exploded
    :param output_file: The name of the output file
    """
    
    # Config
    config_ro = {
        "access": "RO",
        "threads": threads
    }
    config_threads = {
        "threads": threads
    }
    param_explode = {
        "explode_infos": True,
        "export_extra_infos": True
    }

    # List of exploded parquet files
    list_of_exploded_files = []

    # Create 
    parquet_input = Variants(input=input_file, config=config_ro)
    parquet_chromosomes = parquet_input.get_query_to_df(query=f""" SELECT "#CHROM" FROM '{input_file}' WHERE "#CHROM" NOT NULL GROUP BY "#CHROM" """)
    for chrom in parquet_chromosomes["#CHROM"]:
        
        # Split Chrom
        query_chrom = f""" SELECT * FROM '{input_file}' WHERE "#CHROM"='{chrom}'  """
        query_chrom_output = f"{output_file}.{chrom}.parquet"
        query_chrom_explode_output = f"{output_file}.{chrom}.explode.parquet"

        # Log
        log.info(f"Explode infos chromosome {chrom}")


        # Extract
        log.debug(f"Extract chromosome {chrom}: {query_chrom_explode_output}")
        remove_if_exists([query_chrom_output,query_chrom_output+".hdr"])
        parquet_input.export_output(output_file=query_chrom_output, query=query_chrom, export_header=True)
        
        # Explode
        log.debug(f"Explode infos for chromosome {chrom}: {query_chrom_explode_output}")
        list_of_exploded_files.append(query_chrom_explode_output)
        parquet_explode = Variants(input=query_chrom_output, output=query_chrom_explode_output, config=config_threads, param=param_explode)
        parquet_explode.load_data()
        remove_if_exists([query_chrom_explode_output,query_chrom_explode_output+".hdr"])
        parquet_explode.export_output()
        remove_if_exists([query_chrom_output,query_chrom_output+".hdr"])


    # list_of_exploded_files
    log.info(f"Merge explode infos")
    log.debug(f"Merge explode infos files: {list_of_exploded_files}")
    query_explode = f""" SELECT * FROM read_parquet({list_of_exploded_files}) """
    query_explode_output = output_file
    remove_if_exists([query_explode_output,query_explode_output+".hdr"])

    parquet_explode.export_output(output_file=query_explode_output, query=query_explode, export_header=True)

    # Clear
    for file in list_of_exploded_files:
        remove_if_exists([file,file+".hdr"])

    df = parquet_input.get_query_to_df(f"SELECT * FROM '{query_explode_output}' LIMIT 10 ")
    log.debug(df)



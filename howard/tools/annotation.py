import argparse
import logging as log

from howard.functions.commons import load_args, load_config_args, launch_pipeline, load_param_and_config
from howard.objects.variants import Variants
from howard.tools.process import process


def annotation(args: argparse) -> None:
    """
    The `annotation` function performs annotation on a VCF file based on specified parameters and
    exports the annotated data.

    :param args: The `args` parameter is likely an object or dictionary containing various arguments
    passed to the `annotation` function. It is not clear from the code snippet what specific arguments
    are expected or required
    :type args: argparse
    """

    # Load args, param, config and vcfdata_obj
    #arguments_dict, config, param, vcfdata_obj = load_param_and_config(args=args, command="annotation", strict=False, load_data=False)

    vcfdata_obj = process(args=args, tools=["annotation"])

    return vcfdata_obj

    # log.info("Start")

    # # Load config args
    # arguments_dict, _, config, param = load_config_args(args)

    # # Create variants object
    # vcfdata_obj = Variants(
    #     input=args.input, output=args.output, config=config, param=param
    # )

    # # Get Config and Params
    # config = vcfdata_obj.get_config()
    # param = vcfdata_obj.get_param()

    # # Load args into param
    # param = load_args(
    #     param=param,
    #     args=args,
    #     arguments_dict=arguments_dict,
    #     command="annotation",
    #     strict=False,
    # )

    # # Re-Load Config and Params
    # vcfdata_obj.set_param(param)
    # vcfdata_obj.set_config(config)

    # # Load data
    # vcfdata_obj.load_data()


    # # # Pipeline
    # # pipeline = param.get("pipeline", default_pipeline)
    # # steps = pipeline.get("steps", [])
    # # log.debug(f"{param=}")


    # # Start pipeline
    # launch_pipeline(vcfdata_obj=vcfdata_obj, param=param, allowed_tools=["annotation"])
    # # step_i = 0
    # # for step in steps:
    # #     step_i += 1
    # #     log.debug(f"Processing step: {step} [{step_i}/{len(steps)}]")
    # #     for step_name in step:
    # #         step_tool = step.get(step_name, "annotation")
    # #         if step_tool in ["annotation"]:
    # #             log.info(f"Processing pipeline [{step_i}/{len(steps)}] - '{step_name}' [{step_tool}]...")
    # #             if step_name not in param :
    # #                 log.warning(f"Processing pipeline [{step_i}/{len(steps)}] - '{step_name}' [{step_tool}] - Not found in parameters. Try without...")
    # #             try:
    # #                 test = eval(f"vcfdata_obj.{step_tool}(section='{step_name}')")
    # #                 log.debug(f"Processing pipeline [{step_i}/{len(steps)}] - '{step_name}' [{step_tool}] completed successfully.")
    # #                 log.debug(f"Result of step '{step_name}' with tool '{step_tool}': {test}")
    # #             except Exception as e:
    # #                 log.error(f"Error processing step '{step_name}' with tool '{step_tool}': {str(e)}")
    # #                 raise ValueError(f"Error processing step '{step_name}' with tool '{step_tool}': {str(e)}")
    # #             # else:
    # #             #     log.warning(f"Processing pipeline [{step_i}/{len(steps)}] - '{step_name}' [{step_tool}] not found in parameters, skipping.")
    # #         else:
    # #             log.warning(f"Processing pipeline [{step_i}/{len(steps)}] - '{step_name}' [{step_tool}] not in ['annotation'], skipping.")


    # # # Annotation
    # # vcfdata_obj.annotation()

    # # Export
    # vcfdata_obj.export_output()

    # # Log
    # log.info("End")

    # # Return Variants object
    # return vcfdata_obj

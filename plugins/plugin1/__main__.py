import argparse


arguments = {}
commands_arguments = {
    "plugin1": {
        "function": "plugin1",
        "description": """plugin1 description""",
        "help": """Help for plugin1""",
        "epilog": """Usage examples:\n"""
        """   howard plugin1 --input=tests/data/example.vcf.gz \n"""
        """    \n""",
        "groups": {"main": {"input": True, "output": False, "param": False}},
    }
}


def main(args: argparse) -> None:
    print("START Plugin1")
    # print(f"args={args}")
    print(f"Input file: {args.input}")

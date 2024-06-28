# HOWARD Help

How to generate doc

## Generate MarkDown Help from help and JSON

```
howard help --help_md=docs/help.md
howard help --help_json_input=docs/json/help.config.json --help_json_input_title='HOWARD Configuration' --help_md=docs/help.config.md
howard help --help_json_input=docs/json/help.param.json --help_json_input_title='HOWARD Parameters' --help_md=docs/help.param.md
howard help --help_json_input=docs/json/help.param.databases.json --help_json_input_title='HOWARD Parameters Databases' --help_md=docs/help.param.databases.md
```

## Generate HTML and PDF

Use VSCode `Markdown PDF` plugin to generate HTML and PDF from MarkDown files

## Generate pydoc html

```bash
./docs/pdoc.sh
```

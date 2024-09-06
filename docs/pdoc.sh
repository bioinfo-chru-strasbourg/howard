
folder=docs/pdoc
mkdir -p $folder
pip install pdoc==14.5.1
export PDOC_DISPLAY_ENV_VARS=1
pdoc -o $folder howard

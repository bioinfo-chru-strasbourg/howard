# Generate all docs and help
# Usage: ./docs/docs.sh

# Script docs_folder
# Absolute path to this script. /home/user/bin/foo.sh
script=$(readlink -f $0)
# Absolute path this script is in. /home/user/bin
script_path=`dirname $script`

# Init
docs_folder=$script_path
# html_folder=$docs_folder/html
# pdf_folder=$docs_folder/pdf
# html_folder=$docs_folder
# pdf_folder=$docs_folder
list=""

# Folders
# mkdir -p $docs_folder $html_folder $pdf_folder

# Generate HOWARD Main Help
howard help --help_md=$docs_folder/help.md --help_html=$docs_folder/help.html --help_pdf=$docs_folder/help.pdf --code_type="json"
list=$list" $docs_folder/help.md"

# # Generate HOWARD Help in JSON
# for json_file in $docs_folder/json/*.json; do
#     # Init files names
#     md_file=$docs_folder"/"$(basename $json_file | sed s/json\$/md/gi)
#     html_file=$docs_folder"/"$(basename $json_file | sed s/json\$/html/gi)
#     pdf_file=$docs_folder"/"$(basename $json_file | sed s/json\$/pdf/gi)
#     # Init Help title
#     help_title=""
#     for var in $(basename $json_file | sed s/.json\$//gi | tr "." " " | tr "_" " "); do
#         help_title=$help_title" ${var^}"
#     done;
#     # Generate
#     howard help --help_json_input=$json_file --help_md=$md_file --help_html=$html_file --help_pdf=$pdf_file --help_json_input_title="HOWARD$help_title" --code_type="json"
#     list=$list" $md_file"
# done;

# # Generate HOWARD Help in JSON
# for md_file in $docs_folder/*.md; do
# # Check if file already processed
#     if echo "$list" | grep -q "$md_file"; then
#         echo $md_file" already exists";
#     else
#         # Init files names
#         html_file=$(dirname $md_file)"/"$(basename $md_file | sed s/md\$/html/gi)
#         pdf_file=$(dirname $md_file)"/"$(basename $md_file | sed s/md\$/pdf/gi)
#         # Init Help title
#         help_title=""
#         for var in $(basename $md_file | sed s/.md\$//gi | tr "." " " | tr "_" " "); do
#             help_title=$help_title" ${var^}"
#         done;
#         # Generate
#         howard help --help_md_input=$md_file --help_html=$html_file --help_pdf=$pdf_file --help_json_input_title="HOWARD$help_title" --code_type="json"
#     fi;
# done;

# # Generate HOWARD Main README
# for md_file in $script_path/../*.md; do
# # Check if file already processed
#     if echo "$list" | grep -q "$md_file"; then
#         echo $md_file" already exists";
#     else
#         # Init files names
#         html_file=$(dirname $md_file)"/"$(basename $md_file | sed s/md\$/html/gi)
#         pdf_file=$(dirname $md_file)"/"$(basename $md_file | sed s/md\$/pdf/gi)
#         # Init Help title
#         help_title=""
#         for var in $(basename $md_file | sed s/.md\$//gi | tr "." " " | tr "_" " "); do
#             help_title=$help_title" ${var^}"
#         done;
#         # Generate
#         howard help --help_md_input=$md_file --help_html=$html_file --help_pdf=$pdf_file --help_json_input_title="HOWARD$help_title" --code_type="json"
#     fi;
# done;

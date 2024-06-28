#!/bin/bash -ue
mkdir -p /workspaces/MetaChlam/output/test
python /workspaces/MetaChlam/bin/report-lin.py --lin_file /workspaces/MetaChlam/bin/lingroups.txt         --data_file /workspaces/MetaChlam/databases/LINtax_db/taxonomy/data.txt         --in_file_report test.kreport         --in_file_output test.koutput         --output test.LINreport.txt

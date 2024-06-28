#!/bin/bash -ue
python /workspaces/MetaChlam/bin/global_report.py -st test.strains.tsv         -lr test.LINreport.txt -sc strainscan_output         -sm test.gather -s test         -o test_summary.csv

#!/bin/bash -ue
mkdir -p /workspaces/MetaChlam/output/test
kraken2 --db /workspaces/MetaChlam/databases/LINtax_db --paired test_R1.fastq test_R2.fastq         --minimum-hit-groups 4         --confidence 0.45         --output test.koutput         --report test.kreport

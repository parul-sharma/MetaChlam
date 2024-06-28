#!/bin/bash -ue
mkdir -p /workspaces/MetaChlam/output/test
strainscan -i test_R1.fastq -j test_R2.fastq -d /workspaces/MetaChlam/databases/strainscan_db -o strainscan_output

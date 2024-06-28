#!/bin/bash -ue
mkdir -p /workspaces/MetaChlam/output/test
straingst kmerize -k 23 -o test.hdf5 test_R1.fastq test_R2.fastq

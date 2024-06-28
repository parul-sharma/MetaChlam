#!/bin/bash -ue
mkdir -p /workspaces/MetaChlam/output/test
sourmash sketch dna -p scaled=1000,k=31 test_R1.fastq test_R2.fastq -o test.zip --name test

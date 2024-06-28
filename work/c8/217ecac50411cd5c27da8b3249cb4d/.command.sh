#!/bin/bash -ue
mkdir -p /workspaces/MetaChlam/output/test
straingst run -O -o test /workspaces/MetaChlam/databases/pan-genome-db_99.hdf5 test.hdf5

#!/bin/bash -ue
mkdir -p /workspaces/MetaChlam/output/test
sourmash gather test.zip /workspaces/MetaChlam/databases/sourmash_db.sbt.zip > test.gather

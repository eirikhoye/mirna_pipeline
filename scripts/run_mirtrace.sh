#!/bin/bash

# Usage : bash scripts/run_mirtracesh <path/to/config_file>

dt=`date '+%d%m%Y_%H%M%S'`

mirtrace qc --species hsa --custom-db-folder custom_databases/ --config "$1" --write-fasta --uncollapse-fasta -o data/mirtrace_out/"$dt"

